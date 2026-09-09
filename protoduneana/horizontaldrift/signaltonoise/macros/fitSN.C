////////////////////////////////////////////////////////////////////////
// fitSN.C
//
// Angle correction and MPV fit for the ProtoDUNE-HD S/N measurement.
// Consumes the tree written by PDHDSignalToNoise_module.cc.
//
//   root -l -b -q 'fitSN.C("sn_*.root","filt")'
//
// `files` may be a single file, a wildcard, or a .txt file listing one path
// per line, so a grid campaign can be analysed without hadd'ing first.
//
// The final angle cuts live here rather than in the module. Because collection
// and induction wires differ by 35.7 degrees, no track is perpendicular to all
// three planes at once, so each plane gets its own sample: a common cut on the
// angle out of the anode plane (theta_drift), plus a per-plane cut on the
// in-plane angle to that plane's perpendicular-to-wire axis (phi_wire).
//
// The angle correction MUST be derived on a loose sample. Under the final
// cuts the pitch spans only ~0.475-0.49 cm, which is no lever arm at all: a
// free-slope fit there is degenerate, returns an unphysical negative slope,
// and has a denominator that crosses zero just outside the range.
////////////////////////////////////////////////////////////////////////

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {
  /// Resolve LanGausFit.C from the mrb source tree, so the macro runs from any
  /// working directory. Override with SN_LANGAUSFIT.
  TString LanGausFitPath()
  {
    if (const char* e = gSystem->Getenv("SN_LANGAUSFIT")) return TString(e);
    if (const char* s = gSystem->Getenv("MRB_SOURCE"))
      return TString(s) +
             "/protoduneana/protoduneana/singlephase/michelremoving/scripts/LanGausFit.C";
    return TString("LanGausFit.C");
  }
}

TF1* langaufit(TH1*, Double_t*, Double_t*, Double_t*, Double_t*, Double_t*, Double_t*, Double_t*,
               Int_t*, Int_t*);
Int_t langaupro(Double_t*, Double_t&, Double_t&);

namespace {

  double Percentile(std::vector<double> v, double frac)
  {
    if (v.empty()) return -1.;
    std::sort(v.begin(), v.end());
    std::size_t i = static_cast<std::size_t>(frac * v.size());
    if (i >= v.size()) i = v.size() - 1;
    return v[i];
  }

  struct FitOut {
    double mpv = -1.;
    double mean = -1.;
    double fwhm = -1.;
    long n = 0;
    bool ok = false;
  };

  /// Landau (x) Gaussian MPV. Returns ok=false rather than a sentinel number
  /// dressed up as a measurement.
  FitOut FitOne(TH1F* h, long minEntries)
  {
    FitOut o;
    o.n = static_cast<long>(h->GetEntries());
    o.mean = h->GetMean();
    if (o.n < minEntries || h->GetRMS() <= 0.) return o;

    Double_t fr[2] = {h->GetMean() - 2.0 * h->GetRMS(), h->GetMean() + 4.0 * h->GetRMS()};
    if (fr[0] < 0.) fr[0] = 0.;

    // par: [0] Landau width, [1] MPV, [2] area, [3] Gaussian sigma
    Double_t sv[4] = {0.2 * h->GetRMS(), h->GetMean(), h->GetEntries() * h->GetRMS(),
                      0.4 * h->GetRMS()};
    Double_t pllo[4] = {0.001, 0.1 * h->GetMean(), 1., 0.001};
    Double_t plhi[4] = {5.0 * h->GetRMS(), 3.0 * h->GetMean(), 1e9, 5.0 * h->GetRMS()};
    Double_t fp[4], fe[4], chisqr;
    Int_t ndf, status;

    TF1* f = langaufit(h, fr, sv, pllo, plhi, fp, fe, &chisqr, &ndf, &status);
    if (!f) return o;

    Double_t maxx, fwhm;
    langaupro(fp, maxx, fwhm);
    if (maxx <= 0. || maxx > 10. * h->GetMean()) return o;
    o.mpv = maxx;
    o.fwhm = fwhm;
    o.ok = true;
    return o;
  }

  void PrintFit(const char* label, int apa, FitOut const& o)
  {
    if (apa < 0)
      printf("%-12s %8ld ", label, o.n);
    else
      printf("%-12s %5d %8ld ", label, apa, o.n);
    if (o.ok)
      printf("%8.2f %8.2f\n", o.mpv, o.mean);
    else
      printf("%8s %8.2f   (fit failed / too few)\n", "--", o.mean);
  }

} // namespace

void fitSN(const char* files = "pdhd_signaltonoise.root",
           const char* nick = "filt",
           const char* dir = "sigtonoise",
           double maxThetaDrift = 20.,
           double maxPhiWire = 20.,
           double minLen = 100.,
           bool useCaloPitch = true,
           int runMin = -1,
           int runMax = -1,
           bool empiricalCorr = false,
           long minEntries = 200)
{
  const TString lg = LanGausFitPath();
  if (gSystem->AccessPathName(lg)) {
    printf("cannot find LanGausFit.C at %s -- set SN_LANGAUSFIT\n", lg.Data());
    return;
  }
  gROOT->ProcessLine(Form(".L %s+", lg.Data()));

  // ---- input: file, wildcard, or list ---------------------------------------
  TChain ch(Form("%s/hits", dir));
  TString f(files);
  if (f.EndsWith(".txt")) {
    std::ifstream in(f.Data());
    std::string line;
    int n = 0;
    while (std::getline(in, line))
      if (!line.empty() && line[0] != '#') {
        ch.Add(line.c_str());
        ++n;
      }
    printf("added %d files from list %s\n", n, f.Data());
  }
  else {
    ch.Add(f);
  }
  const Long64_t n = ch.GetEntries();
  if (n <= 0) {
    printf("no entries found for %s/hits in %s\n", dir, files);
    return;
  }
  printf("chained %lld hits from %d file(s)\n", n, ch.GetNtrees());

  int plane, apa, status, run;
  float thetaDrift, phiWire, pitchGeo, pitchCalo, trklen, amp, sigma;
  ch.SetBranchAddress("run", &run);
  ch.SetBranchAddress("plane", &plane);
  ch.SetBranchAddress("apa", &apa);
  ch.SetBranchAddress("theta_drift", &thetaDrift);
  ch.SetBranchAddress("phi_wire", &phiWire);
  ch.SetBranchAddress("pitch_geo", &pitchGeo);
  ch.SetBranchAddress("pitch_calo", &pitchCalo);
  ch.SetBranchAddress("trklen", &trklen);
  ch.SetBranchAddress(Form("amp_hit_%s", nick), &amp);
  ch.SetBranchAddress(Form("sigma_%s", nick), &sigma);
  ch.SetBranchAddress(Form("status_%s", nick), &status);

  // pitch_calo carries the space-charge correction and so is the truer path
  // length; it agrees with pitch_geo to 1.0000 against the nosce variant and
  // differs by ~1.3% with SCE on. Falls back where calorimetry is absent.
  auto pitchOf = [&]() -> double {
    return (useCaloPitch && pitchCalo > 0.f) ? pitchCalo : pitchGeo;
  };
  auto basicOK = [&]() -> bool {
    if (status != 0 || sigma <= 0.f || amp <= 0.f) return false;
    if (trklen < minLen) return false;
    if (plane < 0 || plane > 2) return false;
    if (runMin >= 0 && run < runMin) return false;
    if (runMax >= 0 && run > runMax) return false;
    return true;
  };

  // ---- pass 1: angle correction, derived on a LOOSE sample ------------------
  // The final cuts leave pitch spanning ~0.475-0.49 cm. Fitting there is
  // degenerate. Use everything the module wrote, where pitch runs out to ~1 cm.
  std::vector<std::vector<double>> loosePitch(3), tightPitch(3);
  TProfile* prof[3];
  for (int p = 0; p < 3; ++p)
    prof[p] = new TProfile(Form("prof%d", p),
                           Form("plane %d;pitch [cm];amplitude [ADC]", p), 60, 0.4, 1.6);

  for (Long64_t i = 0; i < n; ++i) {
    ch.GetEntry(i);
    if (!basicOK()) continue;
    const double pitch = pitchOf();
    if (pitch <= 0.) continue;
    loosePitch[plane].push_back(pitch);
    prof[plane]->Fill(pitch, amp);
    if (thetaDrift >= 0.f && thetaDrift <= maxThetaDrift && phiWire >= 0.f &&
        phiWire <= maxPhiWire)
      tightPitch[plane].push_back(pitch);
  }

  double pref[3], kprop[3], c0[3], c1[3];
  printf("\n=== angle correction (derived on the loose sample) ===\n");
  printf("%-6s %10s %10s %22s %22s\n", "plane", "N loose", "N tight", "loose pitch range",
         "tight pitch range");
  for (int p = 0; p < 3; ++p) {
    pref[p] = Percentile(tightPitch[p], 0.01);
    if (pref[p] <= 0.) pref[p] = Percentile(loosePitch[p], 0.01);
    printf("%-6s %10zu %10zu   %6.3f - %-10.3f   %6.3f - %-10.3f\n", p == 0 ? "U" : (p == 1 ? "V" : "Coll"),
           loosePitch[p].size(), tightPitch[p].size(), Percentile(loosePitch[p], 0.01),
           Percentile(loosePitch[p], 0.99), Percentile(tightPitch[p], 0.01),
           Percentile(tightPitch[p], 0.99));

    kprop[p] = 0.;
    c0[p] = 0.;
    c1[p] = 0.;
    if (prof[p]->GetEntries() > 100) {
      // Physically motivated: charge on a wire scales with path length through
      // it, so amplitude should be proportional to pitch through the origin.
      TF1 prop("prop", "[0]*x", 0.4, 1.6);
      prof[p]->Fit(&prop, "QNR");
      kprop[p] = prop.GetParameter(0);
      // Free-slope alternative, reported so the assumption can be checked.
      TF1 lin("lin", "pol1", 0.4, 1.6);
      prof[p]->Fit(&lin, "QNR");
      c0[p] = lin.GetParameter(0);
      c1[p] = lin.GetParameter(1);
    }
  }
  printf("\n%-6s %14s %26s %10s\n", "plane", "prop: amp=k*p", "linear: amp=c0+c1*p", "ref pitch");
  for (int p = 0; p < 3; ++p)
    printf("%-6s %14.1f %13.1f + %10.1f*p %10.4f\n", p == 0 ? "U" : (p == 1 ? "V" : "Coll"),
           kprop[p], c0[p], c1[p], pref[p]);
  printf("\nusing the %s correction\n", empiricalCorr ? "linear (empirical)" : "proportional");
  for (int p = 0; p < 3; ++p)
    if (c1[p] <= 0.)
      printf("  NOTE plane %d: linear slope is %.1f (<=0, unphysical) -- "
             "do not trust empiricalCorr here\n",
             p, c1[p]);

  // ---- pass 2: angle-corrected S/N ------------------------------------------
  // APA0 is fitted separately from APAs 1-3. Its V and collection planes have
  // exchanged responses (the known APA0 V/W swap), which puts its collection
  // S/N about an order of magnitude below the others. Pooling all four APAs
  // makes the collection distribution bimodal, and a single Landau (x) Gaussian
  // then describes neither peak -- the symptom is a fitted MPV above the
  // histogram mean, which a Landau tail should make impossible.
  TH1F* hA123[3]; // APAs 1,2,3 combined -- the number to quote
  TH1F* hA0[3];   // APA0 alone
  TH1F* hApa[3][4];
  for (int p = 0; p < 3; ++p) {
    hA123[p] = new TH1F(Form("sn_p%d_apa123", p), Form("plane %d APA1-3;S/N;hits", p), 120, 0., 120.);
    hA0[p] = new TH1F(Form("sn_p%d_apa0", p), Form("plane %d APA0;S/N;hits", p), 120, 0., 120.);
    for (int a = 0; a < 4; ++a)
      hApa[p][a] = new TH1F(Form("sn_p%d_apa%d", p, a),
                            Form("plane %d APA %d;S/N;hits", p, a), 120, 0., 120.);
  }

  for (Long64_t i = 0; i < n; ++i) {
    ch.GetEntry(i);
    if (!basicOK()) continue;
    if (thetaDrift < 0.f || thetaDrift > maxThetaDrift) continue;
    if (phiWire < 0.f || phiWire > maxPhiWire) continue;
    const double pitch = pitchOf();
    if (pitch <= 0.) continue;

    double corr = 1.;
    if (empiricalCorr) {
      const double den = c0[plane] + c1[plane] * pitch;
      const double num = c0[plane] + c1[plane] * pref[plane];
      if (den > 0. && num > 0.) corr = num / den;
    }
    else {
      corr = pref[plane] / pitch; // amp proportional to path length
    }

    const double sn = amp * corr / sigma;
    if (apa == 0)
      hA0[plane]->Fill(sn);
    else if (apa >= 1 && apa <= 3)
      hA123[plane]->Fill(sn);
    if (apa >= 0 && apa < 4) hApa[plane][apa]->Fill(sn);
  }

  // ---- results --------------------------------------------------------------
  const char* pname[3] = {"U", "V", "Collection"};

  printf("\n=== angle-corrected S/N, %s waveforms ===\n", nick);
  printf("APAs 1-3 combined -- the numbers to quote\n");
  printf("%-12s %8s %8s %8s\n", "plane", "hits", "MPV", "mean");
  for (int p = 0; p < 3; ++p)
    PrintFit(pname[p], -1, FitOne(hA123[p], minEntries));

  printf("\nAPA0 alone (V/W responses exchanged -- not comparable to the above)\n");
  printf("%-12s %8s %8s %8s\n", "plane", "hits", "MPV", "mean");
  for (int p = 0; p < 3; ++p)
    PrintFit(pname[p], -1, FitOne(hA0[p], minEntries));

  printf("\nper-APA breakdown (diagnostic)\n");
  printf("%-12s %5s %8s %8s %8s\n", "plane", "APA", "hits", "MPV", "mean");
  for (int p = 0; p < 3; ++p)
    for (int a = 0; a < 4; ++a) {
      FitOut o = FitOne(hApa[p][a], minEntries);
      if (o.n == 0) continue;
      PrintFit(pname[p], a, o);
    }

  TCanvas* c = new TCanvas("c", "S/N", 1500, 1200);
  c->Divide(3, 3);
  for (int p = 0; p < 3; ++p) {
    c->cd(p + 1);
    hA123[p]->Draw();
    c->cd(p + 4);
    hA0[p]->Draw();
    c->cd(p + 7);
    prof[p]->Draw();
  }
  c->SaveAs(Form("sn_%s.png", nick));
}
