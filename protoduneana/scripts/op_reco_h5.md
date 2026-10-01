# Optical reco HDF5 reader notes

`op_reco_h5.py` reads HDF5 files produced by `compare_op_reco.py`.

## Layout

The file is organized by top-level product type and optical channel:

```text
/opdetwf/<channel>/timestamps
/opdetwf/<channel>/waveforms
/opwf/<channel>/timestamps
/opwf/<channel>/waveforms
/ophit/<channel>/hits
```

`opdetwf` contains raw `raw::OpDetWaveform` ADC samples. `opwf` contains
deconvolved `recob::OpWaveform` samples. `ophit` contains `recob::OpHit` values.

The `ophit` dataset is a plain floating-point matrix. Its columns are:

```python
["peak_time", "start_time", "width", "amplitude", "area", "pe"]
```

## Time convention

All timestamp-like quantities are in microseconds.

For raw waveforms, the reference is not the Unix epoch. The original timestamp
was a `uint64` counter in 16 ns ticks since the epoch. Before storage in
`OpDetWaveform` as a double, it was trimmed to the 40 least significant bits
and then converted to microseconds. Therefore these timestamps are relative to
that trimmed 40-bit reference.

The deconvolved waveform (`OpWaveform`) times use the same convention, copied
from `OpDetWaveform`. `OpHit` timing values use the same convention.

Waveform lengths are counts of optical-detector samples. The sample period is
16 ns, or `0.016` microseconds. A waveform starting at timestamp `t0` with
`n` samples covers the half-open interval:

```python
t0 <= time < t0 + n * 0.016
```

## Examples

Read one deconvolved waveform channel into NumPy:

```python
from op_reco_h5 import read_waveforms_numpy

ch0 = read_waveforms_numpy("test.h5", "opwf", 0)
timestamps = ch0["timestamps"]
waveforms = ch0["waveforms"]
```

Read one hit channel with named NumPy fields:

```python
from op_reco_h5 import read_hits_numpy

hits = read_hits_numpy("test.h5", 0, structured=True)
peak_time = hits["peak_time"]
pe = hits["pe"]
```

Read all hit channels as Awkward. The outer array is per channel, and each
channel's `hits` field is an array of hit records:

```python
from op_reco_h5 import read_all_hits_awkward

ophit = read_all_hits_awkward("test.h5")
print(ophit.type)
print(ophit[0].hits[0].peak_time)
```

Associate hits to deconvolved waveforms by comparing hit peak times to waveform
time ranges. For this file's streaming raw waveforms, the matching waveform
channel is `(hit_channel % 40) + 120`. The stored streaming timestamp is about
208 microseconds after the first stored sample, so the reader applies
`timestamp_offset_us=-208.0` by default:

```python
from op_reco_h5 import associate_all_hits_to_waveforms_awkward

associated = associate_all_hits_to_waveforms_awkward("test.h5")
print(associated.type)

first_hit = associated[0].hits[0]
print(first_hit.waveform_index)  # index in this channel's waveform array, or -1
print(first_hit.sample_index)    # sample containing the hit peak, or -1
```

To use same-numbered self-triggered waveform snippets instead, override the
defaults:

```python
associated = associate_all_hits_to_waveforms_awkward(
    "test.h5",
    collection="opwf",
    channel_map="same",
    timestamp_offset_us=0.0,
)
```

Read all deconvolved waveforms as a ragged Awkward array:

```python
from op_reco_h5 import read_all_waveforms_awkward

opwf = read_all_waveforms_awkward("test.h5", "opwf")
print(opwf.channel)
print(opwf[0].waveforms)
```

Read every available collection:

```python
from op_reco_h5 import read

arrays = read("test.h5", library="awkward")
raw = arrays["opdetwf"]
deco = arrays["opwf"]
hits = arrays["ophit"]
```
