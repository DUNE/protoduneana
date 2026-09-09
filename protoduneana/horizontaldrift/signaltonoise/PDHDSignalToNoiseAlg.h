////////////////////////////////////////////////////////////////////////
// Class:       PDHDSignalToNoiseAlg
// File:        PDHDSignalToNoiseAlg.h
//
// Signal and noise extraction for the ProtoDUNE-HD S/N measurement.
//
// Deliberately free of art::Event so the same code can be driven either by
// a single job reading a reco file with the raw file as a secondary input,
// or by a two-pass arrangement.
//
// The noise definition follows section 4.6 of the ProtoDUNE-SP performance
// paper: the standard deviation of a Gaussian fit to the distribution of raw
// ADC values in signal-free regions of a channel's waveform.
////////////////////////////////////////////////////////////////////////

#ifndef PDHDSIGNALTONOISEALG_H
#define PDHDSIGNALTONOISEALG_H

#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Wire.h"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace pdhd {

  /// Per-channel noise, from a Gaussian fit to signal-free raw ADC values.
  struct NoiseResult {
    float ped = -999.f;       ///< Gaussian mean, i.e. the pedestal
    float sigma = -999.f;     ///< Gaussian sigma, i.e. the noise
    float rms = -999.f;       ///< plain RMS over the same ticks (QA cross-check)
    float chi2ndf = -999.f;
    int nfree = 0;            ///< ticks that went into the fit
    int fitstatus = -1;       ///< 0 == converged
    bool ok = false;
  };

  /// Signal / signal-free structure of one deconvolved wire.
  struct WireMask {
    std::vector<bool> isSignal; ///< inside an ROI, with a non-zero value
    std::vector<bool> isFree;   ///< signal-free after short-run removal and erosion
    int nfree = 0;
  };

  /// Maximum of a raw waveform over some tick window.
  struct MaxResult {
    float amp = -999.f; ///< maximum, pedestal subtracted
    float raw = -999.f; ///< maximum, as read from the digit
    int tick = -1;
  };

  class PDHDSignalToNoiseAlg {
  public:
    explicit PDHDSignalToNoiseAlg(fhicl::ParameterSet const& p);

    /// Signal / signal-free masks for one deconvolved wire.
    ///
    /// Works for both dense and sparse WireCell output. PDHD's wclsdatahd is
    /// configured dense (signal_output_form: "dense"), so SignalROI() holds a
    /// single range spanning the whole readout and it is the exact zeros
    /// within it that delimit the regions of interest.
    WireMask MakeWireMask(recob::Wire const& wire, std::size_t nticks) const;

    /// Half-open [lo,hi) extent of the ROI containing `tick`, or {-1,-1}.
    std::pair<int, int> RoiContaining(WireMask const& mask, int tick) const;

    /// Gaussian fit to the raw ADC distribution on the mask's free ticks.
    ///
    /// The mask is indexed in the deconvolved-wire tick frame, which is not
    /// necessarily the digit's own. `scale` is nDigitTicks / nWireTicks, so
    /// digit tick t corresponds to wire tick t / scale. For PDHD this is 1 for
    /// wclsdatahdfilter:raw (5999 ticks, 500 ns) and 5859/5999 = 500/512 for
    /// tpcrawdecoder:daq (5859 ticks, 512 ns), which are the same 3 ms window
    /// sampled differently.
    NoiseResult FitNoise(std::vector<short> const& adcs,
                         WireMask const& mask,
                         double scale = 1.0) const;

    /// Wire-frame tick -> digit-frame tick.
    static int ToDigitTick(double wireTick, double scale)
    {
      return static_cast<int>(std::lround(wireTick * scale));
    }

    /// Maximum raw ADC over the half-open window [lo,hi), pedestal subtracted.
    MaxResult MaxInWindow(std::vector<short> const& adcs, int lo, int hi, float ped) const;

    int HitHalfWindow() const { return fHitHalfWindow; }
    unsigned MinFreeTicks() const { return fMinFreeTicks; }

  private:
    unsigned fMinZeroRun;    ///< shortest run of zeros that counts as a real gap
    unsigned fGuardBand;     ///< ticks trimmed from each end of a free run
    unsigned fMinFreeTicks;  ///< below this the channel's noise is unusable
    int fHitHalfWindow;      ///< +/- ticks searched around a hit peak time
    float fFitRangeSigma;    ///< fit range, in robust sigmas about the median
  };

} // namespace pdhd

#endif
