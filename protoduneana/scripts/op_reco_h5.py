#!/usr/bin/env python3
"""Readers for HDF5 files written by ``compare_op_reco.py``.

The producer script stores optical detector data in channel-indexed HDF5
groups:

* ``/opdetwf/<channel>/timestamps`` and ``waveforms`` for raw
  ``raw::OpDetWaveform`` objects.
* ``/opwf/<channel>/timestamps`` and ``waveforms`` for deconvolved
  ``recob::OpWaveform`` objects.
* ``/ophit/<channel>/hits`` for ``recob::OpHit`` objects.

Time convention
---------------
All timestamp-like quantities are in microseconds. The raw waveform timestamp
origin is not the Unix epoch: the original 16 ns tick counter was trimmed to
its 40 least significant bits, then converted to microseconds and stored as a
double in ``OpDetWaveform``. ``OpWaveform`` and ``OpHit`` times use the same
trimmed-reference convention copied through from ``OpDetWaveform``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Literal, Union

import awkward as ak
import h5py
import numpy as np

WaveformCollection = Literal["opdetwf", "opwf"]
ArrayLibrary = Literal["numpy", "awkward"]
ChannelMapName = Literal["same", "mod40_plus120"]
ChannelMap = Union[ChannelMapName, Mapping[int, int], Callable[[int], int]]

HIT_FIELDS: tuple[str, ...] = (
    "peak_time",
    "start_time",
    "width",
    "amplitude",
    "area",
    "pe",
)
ASSOCIATION_FIELDS: tuple[str, ...] = ("waveform_index", "sample_index")
ASSOCIATED_HIT_FIELDS: tuple[str, ...] = HIT_FIELDS + ASSOCIATION_FIELDS
SAMPLE_PERIOD_US = 0.016
STREAMING_APA0_TIMESTAMP_OFFSET_US = -208.0

HIT_UNITS: Mapping[str, str] = {
    "peak_time": "microseconds, trimmed 40-bit reference",
    "start_time": "microseconds, trimmed 40-bit reference",
    "width": "microseconds",
    "amplitude": "producer units from recob::OpHit::Amplitude()",
    "area": "producer units from recob::OpHit::Area()",
    "pe": "photoelectrons",
}

TIME_CONVENTION = (
    "Times are stored in microseconds using the trimmed optical-detector "
    "reference: the original uint64 counter in 16 ns ticks was reduced to its "
    "40 least significant bits before conversion to microseconds. These times "
    "are therefore not measured from the Unix epoch."
)


def open_h5(file_or_path: str | Path | h5py.File):
    """Return a context manager yielding an ``h5py.File``.

    Existing ``h5py.File`` objects are not closed by this helper. Paths are
    opened read-only and closed when the context exits.
    """
    if isinstance(file_or_path, h5py.File):
        return nullcontext(file_or_path)
    return h5py.File(file_or_path, "r")


def channels(file_or_path: str | Path | h5py.File, group: str) -> list[int]:
    """Return sorted integer channel IDs present under an HDF5 group."""
    with open_h5(file_or_path) as h5:
        if group not in h5:
            return []
        return sorted(int(key) for key in h5[group].keys())


def available_groups(file_or_path: str | Path | h5py.File) -> list[str]:
    """Return top-level HDF5 groups available in the file."""
    with open_h5(file_or_path) as h5:
        return sorted(key for key, value in h5.items() if isinstance(value, h5py.Group))


def read_waveforms_numpy(
    file_or_path: str | Path | h5py.File,
    collection: WaveformCollection,
    channel: int,
) -> dict[str, np.ndarray]:
    """Read one channel of waveform data as NumPy arrays.

    Parameters
    ----------
    file_or_path:
        HDF5 path or open ``h5py.File``.
    collection:
        ``"opdetwf"`` for raw ADC waveforms, or ``"opwf"`` for deconvolved
        floating-point waveforms.
    channel:
        Optical channel number.

    Returns
    -------
    dict
        ``{"timestamps": ..., "waveforms": ...}``. Timestamps are in
        microseconds using the trimmed-reference convention described in the
        module docstring.
    """
    with open_h5(file_or_path) as h5:
        group = h5[f"/{collection}/{int(channel)}"]
        return {
            "timestamps": group["timestamps"][()],
            "waveforms": group["waveforms"][()],
        }


def read_waveforms_awkward(
    file_or_path: str | Path | h5py.File,
    collection: WaveformCollection,
    channel: int,
) -> ak.Array:
    """Read one channel of waveform data as an Awkward record array."""
    arrays = read_waveforms_numpy(file_or_path, collection, channel)
    return ak.Array(
        {
            "timestamp": arrays["timestamps"],
            "waveform": arrays["waveforms"],
        }
    )


def read_all_waveforms_numpy(
    file_or_path: str | Path | h5py.File,
    collection: WaveformCollection,
    selected_channels: Iterable[int] | None = None,
) -> dict[int, dict[str, np.ndarray]]:
    """Read waveform data for many channels as a dict of NumPy arrays.

    A dict is used instead of one stacked NumPy array because the number of
    waveforms can differ by channel, and raw channels can mix self-triggered
    and streaming waveform lengths.
    """
    with open_h5(file_or_path) as h5:
        channel_ids = (
            sorted(int(chan) for chan in selected_channels)
            if selected_channels is not None
            else channels(h5, collection)
        )
        return {
            chan: {
                "timestamps": h5[f"/{collection}/{chan}/timestamps"][()],
                "waveforms": h5[f"/{collection}/{chan}/waveforms"][()],
            }
            for chan in channel_ids
        }


def read_all_waveforms_awkward(
    file_or_path: str | Path | h5py.File,
    collection: WaveformCollection,
    selected_channels: Iterable[int] | None = None,
) -> ak.Array:
    """Read waveform data for many channels as a ragged Awkward array.

    The result has one top-level record per channel:
    ``channel``, ``timestamps``, and ``waveforms``.
    """
    by_channel = read_all_waveforms_numpy(file_or_path, collection, selected_channels)
    channels = []
    timestamps = []
    waveforms = []
    for chan, arrays in sorted(by_channel.items()):
        channels.append(chan)
        timestamps.append(arrays["timestamps"])
        waveforms.append(arrays["waveforms"])

    return ak.zip(
        {
            "channel": np.asarray(channels, dtype=np.int64),
            "timestamps": _ragged_1d(timestamps),
            "waveforms": _ragged_2d(waveforms),
        },
        depth_limit=1,
    )


def read_hits_numpy(
    file_or_path: str | Path | h5py.File,
    channel: int,
    *,
    structured: bool = False,
) -> np.ndarray:
    """Read ``/ophit/<channel>/hits`` as a NumPy array.

    The unstructured producer layout has columns in ``HIT_FIELDS`` order:
    ``peak_time``, ``start_time``, ``width``, ``amplitude``, ``area``, ``pe``.
    Set ``structured=True`` to return a structured array with those field names.
    """
    with open_h5(file_or_path) as h5:
        hits = h5[f"/ophit/{int(channel)}/hits"][()]

    if not structured:
        return hits

    out = np.empty(hits.shape[0], dtype=[(name, hits.dtype) for name in HIT_FIELDS])
    for index, name in enumerate(HIT_FIELDS):
        out[name] = hits[:, index]
    return out


def read_hits_awkward(file_or_path: str | Path | h5py.File, channel: int) -> ak.Array:
    """Read one channel of optical hits as an Awkward record array."""
    hits = read_hits_numpy(file_or_path, channel)
    return ak.zip({name: hits[:, index] for index, name in enumerate(HIT_FIELDS)})


def read_all_hits_numpy(
    file_or_path: str | Path | h5py.File,
    selected_channels: Iterable[int] | None = None,
    *,
    structured: bool = False,
) -> dict[int, np.ndarray]:
    """Read optical hits for many channels as ``{channel: numpy_array}``."""
    with open_h5(file_or_path) as h5:
        channel_ids = (
            sorted(int(chan) for chan in selected_channels)
            if selected_channels is not None
            else channels(h5, "ophit")
        )
        out = {}
        for chan in channel_ids:
            hits = h5[f"/ophit/{chan}/hits"][()]
            out[chan] = _hits_to_structured(hits) if structured else hits
        return out


def read_all_hits_awkward(
    file_or_path: str | Path | h5py.File,
    selected_channels: Iterable[int] | None = None,
) -> ak.Array:
    """Read optical hits for many channels as a ragged Awkward array.

    The result has one top-level record per channel. Its ``hits`` field is a
    variable-length array of hit records, not a record of variable-length
    arrays.
    """
    by_channel = read_all_hits_numpy(file_or_path, selected_channels)
    channels = []
    hit_records = []
    for chan, hits in sorted(by_channel.items()):
        channels.append(chan)
        hit_records.append(_hit_columns(hits))

    return ak.zip(
        {
            "channel": np.asarray(channels, dtype=np.int64),
            "hits": _ragged_records(hit_records),
        },
        depth_limit=1,
    )


def associate_hits_to_waveforms_numpy(
    file_or_path: str | Path | h5py.File,
    channel: int,
    collection: WaveformCollection = "opdetwf",
    *,
    sample_period_us: float = SAMPLE_PERIOD_US,
    channel_map: ChannelMap = "mod40_plus120",
    timestamp_offset_us: float = STREAMING_APA0_TIMESTAMP_OFFSET_US,
    structured: bool = True,
) -> np.ndarray:
    """Associate one channel's optical hits with waveforms by timestamp.

    By default this uses the long APA0 streaming raw waveforms seen in
    `compare_op_reco.py` outputs: hit channel ``ch`` is matched to waveform
    channel ``(ch % 40) + 120``, and the waveform time window starts 208 us
    before the stored streaming timestamp.

    A hit is matched to waveform ``i`` when ``waveform_start[i] <=
    hit.peak_time < waveform_start[i] + n_samples * sample_period_us``.
    The returned association fields are:

    * ``waveform_index``: index into that channel's waveform array, or ``-1``.
    * ``sample_index``: sample containing the hit peak, or ``-1``.

    Set ``structured=False`` to return a plain 2D array with columns in
    ``ASSOCIATED_HIT_FIELDS`` order.
    """
    waveform_channel = _resolve_waveform_channel(channel, channel_map)
    waveforms = read_waveforms_numpy(file_or_path, collection, waveform_channel)
    hits = read_hits_numpy(file_or_path, channel)
    waveform_index, sample_index = _associate_hit_times_to_waveforms(
        hits[:, HIT_FIELDS.index("peak_time")],
        waveforms["timestamps"] + timestamp_offset_us,
        waveforms["waveforms"].shape[1],
        sample_period_us,
    )
    associated = np.column_stack((hits, waveform_index, sample_index))

    if not structured:
        return associated

    out = np.empty(
        associated.shape[0],
        dtype=[
            *((name, hits.dtype) for name in HIT_FIELDS),
            ("waveform_index", np.int64),
            ("sample_index", np.int64),
        ],
    )
    for index, name in enumerate(HIT_FIELDS):
        out[name] = hits[:, index]
    out["waveform_index"] = waveform_index
    out["sample_index"] = sample_index
    return out


def associate_hits_to_waveforms_awkward(
    file_or_path: str | Path | h5py.File,
    channel: int,
    collection: WaveformCollection = "opdetwf",
    *,
    sample_period_us: float = SAMPLE_PERIOD_US,
    channel_map: ChannelMap = "mod40_plus120",
    timestamp_offset_us: float = STREAMING_APA0_TIMESTAMP_OFFSET_US,
) -> ak.Array:
    """Associate one channel's hits with waveforms as hit records."""
    associated = associate_hits_to_waveforms_numpy(
        file_or_path,
        channel,
        collection,
        sample_period_us=sample_period_us,
        channel_map=channel_map,
        timestamp_offset_us=timestamp_offset_us,
        structured=True,
    )
    return ak.zip({name: associated[name] for name in ASSOCIATED_HIT_FIELDS})


def associate_all_hits_to_waveforms_numpy(
    file_or_path: str | Path | h5py.File,
    collection: WaveformCollection = "opdetwf",
    selected_channels: Iterable[int] | None = None,
    *,
    sample_period_us: float = SAMPLE_PERIOD_US,
    channel_map: ChannelMap = "mod40_plus120",
    timestamp_offset_us: float = STREAMING_APA0_TIMESTAMP_OFFSET_US,
    structured: bool = True,
) -> dict[int, np.ndarray]:
    """Associate hits with waveforms for many channels."""
    hit_channels = set(channels(file_or_path, "ophit"))
    waveform_channels = set(channels(file_or_path, collection))
    if selected_channels is None:
        channel_ids = sorted(
            chan
            for chan in hit_channels
            if _resolve_waveform_channel(chan, channel_map) in waveform_channels
        )
    else:
        channel_ids = sorted(
            int(chan)
            for chan in selected_channels
            if int(chan) in hit_channels
            and _resolve_waveform_channel(int(chan), channel_map) in waveform_channels
        )

    return {
        chan: associate_hits_to_waveforms_numpy(
            file_or_path,
            chan,
            collection,
            sample_period_us=sample_period_us,
            channel_map=channel_map,
            timestamp_offset_us=timestamp_offset_us,
            structured=structured,
        )
        for chan in channel_ids
    }


def associate_all_hits_to_waveforms_awkward(
    file_or_path: str | Path | h5py.File,
    collection: WaveformCollection = "opdetwf",
    selected_channels: Iterable[int] | None = None,
    *,
    sample_period_us: float = SAMPLE_PERIOD_US,
    channel_map: ChannelMap = "mod40_plus120",
    timestamp_offset_us: float = STREAMING_APA0_TIMESTAMP_OFFSET_US,
) -> ak.Array:
    """Associate hits with waveforms for many channels as a ragged array.

    The result has one top-level record per channel. Its ``hits`` field is an
    array of hit records with the normal hit fields plus ``waveform_index`` and
    ``sample_index``.
    """
    by_channel = associate_all_hits_to_waveforms_numpy(
        file_or_path,
        collection,
        selected_channels,
        sample_period_us=sample_period_us,
        channel_map=channel_map,
        timestamp_offset_us=timestamp_offset_us,
        structured=True,
    )
    channels_out = []
    hit_records = []
    for chan, hits in sorted(by_channel.items()):
        channels_out.append(chan)
        hit_records.append({name: hits[name] for name in ASSOCIATED_HIT_FIELDS})

    return ak.zip(
        {
            "channel": np.asarray(channels_out, dtype=np.int64),
            "hits": _ragged_records(hit_records),
        },
        depth_limit=1,
    )


def read(
    file_or_path: str | Path | h5py.File,
    *,
    library: ArrayLibrary = "awkward",
    include_raw: bool = True,
    include_deco: bool = True,
    include_hits: bool = True,
) -> dict[str, Any]:
    """Read all known optical reco collections from a file.

    Parameters select the output library and which top-level collections to
    include. NumPy output is a nested dict keyed by channel. Awkward output is a
    dict whose values are ragged Awkward arrays.
    """
    if library not in ("numpy", "awkward"):
        raise ValueError("library must be 'numpy' or 'awkward'")

    readers = {
        "numpy": {
            "opdetwf": read_all_waveforms_numpy,
            "opwf": read_all_waveforms_numpy,
            "ophit": read_all_hits_numpy,
        },
        "awkward": {
            "opdetwf": read_all_waveforms_awkward,
            "opwf": read_all_waveforms_awkward,
            "ophit": read_all_hits_awkward,
        },
    }[library]

    out: dict[str, Any] = {}
    groups = set(available_groups(file_or_path))
    if include_raw and "opdetwf" in groups:
        out["opdetwf"] = readers["opdetwf"](file_or_path, "opdetwf")
    if include_deco and "opwf" in groups:
        out["opwf"] = readers["opwf"](file_or_path, "opwf")
    if include_hits and "ophit" in groups:
        out["ophit"] = readers["ophit"](file_or_path)
    return out


def _hits_to_structured(hits: np.ndarray) -> np.ndarray:
    out = np.empty(hits.shape[0], dtype=[(name, hits.dtype) for name in HIT_FIELDS])
    for index, name in enumerate(HIT_FIELDS):
        out[name] = hits[:, index]
    return out


def _associate_hit_times_to_waveforms(
    hit_times_us: np.ndarray,
    waveform_starts_us: np.ndarray,
    waveform_length: int | np.ndarray,
    sample_period_us: float,
) -> tuple[np.ndarray, np.ndarray]:
    hit_times_us = np.asarray(hit_times_us)
    waveform_starts_us = np.asarray(waveform_starts_us)
    lengths = np.broadcast_to(np.asarray(waveform_length), waveform_starts_us.shape)
    waveform_stops_us = waveform_starts_us + lengths * sample_period_us

    waveform_index = np.full(hit_times_us.shape, -1, dtype=np.int64)
    sample_index = np.full(hit_times_us.shape, -1, dtype=np.int64)

    for hit_index, hit_time in enumerate(hit_times_us):
        matches = np.flatnonzero(
            (waveform_starts_us <= hit_time) & (hit_time < waveform_stops_us)
        )
        if matches.size == 0:
            continue

        matched_waveform = int(matches[0])
        sample = int(
            np.floor(
                (hit_time - waveform_starts_us[matched_waveform]) / sample_period_us
            )
        )
        if 0 <= sample < lengths[matched_waveform]:
            waveform_index[hit_index] = matched_waveform
            sample_index[hit_index] = sample

    return waveform_index, sample_index


def _resolve_waveform_channel(channel: int, channel_map: ChannelMap) -> int:
    if channel_map == "same":
        return int(channel)
    if channel_map == "mod40_plus120":
        return int(channel) % 40 + 120
    if isinstance(channel_map, Mapping):
        return int(channel_map[int(channel)])
    return int(channel_map(int(channel)))


def _hit_columns(hits: np.ndarray) -> dict[str, np.ndarray]:
    return {name: hits[:, index] for index, name in enumerate(HIT_FIELDS)}


def _ragged_1d(arrays: Iterable[np.ndarray]) -> ak.Array:
    arrays = [np.asarray(array) for array in arrays]
    offsets = np.empty(len(arrays) + 1, dtype=np.int64)
    offsets[0] = 0
    for index, array in enumerate(arrays, start=1):
        offsets[index] = offsets[index - 1] + len(array)

    if arrays:
        content = np.concatenate(arrays)
    else:
        content = np.asarray([], dtype=np.float64)

    layout = ak.contents.ListOffsetArray(
        ak.index.Index64(offsets),
        ak.contents.NumpyArray(content),
    )
    return ak.Array(layout)


def _ragged_2d(arrays: Iterable[np.ndarray]) -> ak.Array:
    arrays = [np.asarray(array) for array in arrays]
    outer_offsets = np.empty(len(arrays) + 1, dtype=np.int64)
    outer_offsets[0] = 0
    inner_lengths = []
    flat_arrays = []

    for index, array in enumerate(arrays, start=1):
        if array.ndim != 2:
            raise ValueError(f"expected 2D waveform array, got shape {array.shape}")
        outer_offsets[index] = outer_offsets[index - 1] + array.shape[0]
        inner_lengths.extend([array.shape[1]] * array.shape[0])
        flat_arrays.append(array.reshape(-1))

    inner_offsets = np.empty(len(inner_lengths) + 1, dtype=np.int64)
    inner_offsets[0] = 0
    for index, length in enumerate(inner_lengths, start=1):
        inner_offsets[index] = inner_offsets[index - 1] + length

    if flat_arrays:
        content = np.concatenate(flat_arrays)
    else:
        content = np.asarray([], dtype=np.float64)

    inner = ak.contents.ListOffsetArray(
        ak.index.Index64(inner_offsets),
        ak.contents.NumpyArray(content),
    )
    outer = ak.contents.ListOffsetArray(ak.index.Index64(outer_offsets), inner)
    return ak.Array(outer)


def _ragged_records(arrays: Iterable[Mapping[str, np.ndarray]]) -> ak.Array:
    arrays = [{name: np.asarray(values) for name, values in array.items()} for array in arrays]
    offsets = np.empty(len(arrays) + 1, dtype=np.int64)
    offsets[0] = 0
    for index, array in enumerate(arrays, start=1):
        lengths = {len(values) for values in array.values()}
        if len(lengths) != 1:
            raise ValueError("record fields must have equal lengths")
        length = lengths.pop() if lengths else 0
        offsets[index] = offsets[index - 1] + length

    field_names = list(arrays[0].keys()) if arrays else []
    flat = {}
    for name in field_names:
        values = [array[name] for array in arrays]
        flat[name] = np.concatenate(values) if values else np.asarray([])

    records = ak.zip(flat).layout
    layout = ak.contents.ListOffsetArray(ak.index.Index64(offsets), records)
    return ak.Array(layout)
