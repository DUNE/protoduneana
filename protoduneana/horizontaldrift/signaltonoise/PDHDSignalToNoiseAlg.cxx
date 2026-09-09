////////////////////////////////////////////////////////////////////////
// PDHDSignalToNoiseAlg.cxx
////////////////////////////////////////////////////////////////////////

#include "protoduneana/horizontaldrift/signaltonoise/PDHDSignalToNoiseAlg.h"

#include "TF1.h"
#include "TH1F.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
  /// Scale factor turning a median absolute deviation into a Gaussian sigma.
  constexpr float kMADtoSigma = 1.4826f;
}

pdhd::PDHDSignalToNoiseAlg::PDHDSignalToNoiseAlg(fhicl::ParameterSet const& p)
  : fMinZeroRun(p.get<unsigned>("MinZeroRun", 20))
  , fGuardBand(p.get<unsigned>("GuardBand", 20))
  , fMinFreeTicks(p.get<unsigned>("MinFreeTicks", 1000))
  , fHitHalfWindow(p.get<int>("HitHalfWindow", 40))
  , fFitRangeSigma(p.get<float>("FitRangeSigma", 4.0))
{}

pdhd::WireMask
pdhd::PDHDSignalToNoiseAlg::MakeWireMask(recob::Wire const& wire, std::size_t nticks) const
{
  WireMask m;
  if (nticks == 0) return m;

  m.isSignal.assign(nticks, false);

  // A tick carries signal only if it lies inside an ROI *and* its deconvolved
  // value is non-zero. Under dense output the single range covers everything,
  // so the zeros do the delimiting; under sparse output the ranges themselves
  // do it and the non-zero test is a no-op. Both cases are handled here.
  for (auto const& range : wire.SignalROI().get_ranges()) {
    std::size_t tick = range.begin_index();
    for (float v : range) {
      if (tick >= nticks) break;
      if (v != 0.0f) m.isSignal[tick] = true;
      ++tick;
    }
  }

  // Free ticks are the non-signal ones, but only where the gap is long enough
  // to be a genuine quiet region -- an isolated zero inside a real ROI must not
  // punch a hole in the mask -- and with a guard band trimmed off each end so
  // that induction tails leaking past the ROI edge do not inflate the noise.
  m.isFree.assign(nticks, false);
  std::size_t i = 0;
  while (i < nticks) {
    if (m.isSignal[i]) {
      ++i;
      continue;
    }
    std::size_t j = i;
    while (j < nticks && !m.isSignal[j]) ++j; // [i,j) is a run of non-signal ticks

    const std::size_t len = j - i;
    if (len >= fMinZeroRun && len > 2u * fGuardBand) {
      for (std::size_t k = i + fGuardBand; k + fGuardBand < j; ++k) m.isFree[k] = true;
    }
    i = j;
  }

  m.nfree = static_cast<int>(std::count(m.isFree.begin(), m.isFree.end(), true));
  return m;
}

std::pair<int, int>
pdhd::PDHDSignalToNoiseAlg::RoiContaining(WireMask const& mask, int tick) const
{
  const int n = static_cast<int>(mask.isSignal.size());
  if (tick < 0 || tick >= n || !mask.isSignal[tick]) return {-1, -1};

  int lo = tick;
  while (lo > 0 && mask.isSignal[lo - 1]) --lo;
  int hi = tick;
  while (hi + 1 < n && mask.isSignal[hi + 1]) ++hi;

  return {lo, hi + 1};
}

pdhd::NoiseResult
pdhd::PDHDSignalToNoiseAlg::FitNoise(std::vector<short> const& adcs,
                                     WireMask const& mask,
                                     double scale) const
{
  NoiseResult r;
  if (!(scale > 0.)) return r;

  const std::size_t nmask = mask.isFree.size();
  std::vector<short> vals;
  vals.reserve(adcs.size());

  // Walk the digit's own ticks and ask the mask, which lives in the wire tick
  // frame, whether the corresponding moment in time is signal-free. With
  // scale == 1 this is a straight index lookup.
  for (std::size_t t = 0; t < adcs.size(); ++t) {
    const long w = std::lround(t / scale);
    if (w < 0 || static_cast<std::size_t>(w) >= nmask) continue;
    if (mask.isFree[w]) vals.push_back(adcs[t]);
  }

  r.nfree = static_cast<int>(vals.size());
  if (vals.size() < fMinFreeTicks) return r;

  // Robust centre and width first, so the fit range is not dragged around by
  // whatever signal did leak through the mask.
  std::vector<short> work = vals;
  const std::size_t mid = work.size() / 2;
  std::nth_element(work.begin(), work.begin() + mid, work.end());
  const float median = work[mid];

  std::vector<float> dev;
  dev.reserve(vals.size());
  for (short v : vals) dev.push_back(std::abs(v - median));
  std::nth_element(dev.begin(), dev.begin() + mid, dev.end());
  const float robustSigma = kMADtoSigma * dev[mid];

  if (!(robustSigma > 0.f)) return r; // dead or stuck channel

  // Plain RMS over the same ticks, as an independent cross-check on the fit.
  double sum = 0., sum2 = 0.;
  for (short v : vals) {
    sum += v;
    sum2 += static_cast<double>(v) * v;
  }
  const double mean = sum / vals.size();
  r.rms = static_cast<float>(std::sqrt(std::max(0., sum2 / vals.size() - mean * mean)));

  // Histogram at 1 ADC binning over a robust window, then fit a Gaussian.
  const float half = std::max(5.0f, fFitRangeSigma * robustSigma);
  const int lo = static_cast<int>(std::floor(median - half));
  const int hi = static_cast<int>(std::ceil(median + half));
  const int nbins = std::max(1, hi - lo);

  TH1F h("h_noise", "", nbins, lo - 0.5, hi + 0.5);
  h.SetDirectory(nullptr);
  for (short v : vals) h.Fill(v);

  TF1 g("g_noise", "gaus", lo - 0.5, hi + 0.5);
  g.SetParameters(h.GetMaximum(), median, robustSigma);

  r.fitstatus = h.Fit(&g, "QNR");
  r.ped = g.GetParameter(1);
  r.sigma = std::abs(g.GetParameter(2));
  if (g.GetNDF() > 0) r.chi2ndf = g.GetChisquare() / g.GetNDF();
  r.ok = (r.fitstatus == 0) && (r.sigma > 0.f);

  return r;
}

pdhd::MaxResult
pdhd::PDHDSignalToNoiseAlg::MaxInWindow(std::vector<short> const& adcs,
                                        int lo,
                                        int hi,
                                        float ped) const
{
  MaxResult r;
  const int n = static_cast<int>(adcs.size());
  lo = std::max(lo, 0);
  hi = std::min(hi, n);

  float best = std::numeric_limits<float>::lowest();
  for (int t = lo; t < hi; ++t) {
    if (adcs[t] > best) {
      best = adcs[t];
      r.tick = t;
    }
  }

  if (r.tick >= 0) {
    r.raw = best;
    r.amp = best - ped;
  }
  return r;
}
