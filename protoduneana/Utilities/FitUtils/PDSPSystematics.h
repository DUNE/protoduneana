#ifndef PDSPSYSTEMATICS_hh
#define PDSPSYSTEMATICS_hh

#include <map>
#include "ThinSliceDriver.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TSpline.h"

namespace protoana {

class PDSPSystematics {
 public:
  PDSPSystematics(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<std::string, ThinSliceSystematic> & pars,
      const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
      TFile & output_file, int upstream_ID, int no_track_ID, int decay_ID,
      int past_FV_ID, int beam_cut_ID, int past_FV_sel_ID);
  virtual ~PDSPSystematics(){};

  double GetEventWeight(
      const ThinSliceEvent & event,
      int signal_index,
      //int selection_bin,
      const std::map<std::string, ThinSliceSystematic> & pars);

  void SetupSyst_G4RWCoeff(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_TiedG4RWCoeff(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSignal_G4RWCoeff(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSignal_TiedG4RWCoeff(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_G4RWCoeff(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par);
      //const std::map<std::string, ThinSliceSystematic> & pars);

  void SetupSyst_CosThetaEff(
      const std::map<std::string, ThinSliceSystematic> & pars);
  static double GetSystWeight_CosThetaEff(
    const ThinSliceEvent & event, int signal_idnex);

  static double GetSystWeight_G4RWCoeffNoPar(
      const ThinSliceEvent & event, int signal_index);
  static double GetSystWeight_TiedG4RWCoeffNoPar(
      const ThinSliceEvent & event, int signal_index);

  static double GetSignalWeight_G4RWCoeffNoPar(
      const ThinSliceEvent & event, int signal_index);
  static double GetSignalWeight_TiedG4RWCoeffNoPar(
      const ThinSliceEvent & event, int signal_index);

  void SetupSyst_BeamShiftBins(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamShiftBin(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par);


  double GetSystWeight_TiedG4RWCoeff(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par);

  void SetupSyst_BeamShift(
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  double GetSystWeight_BeamShift(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par);
      //const std::map<std::string, ThinSliceSystematic> & pars);

  void SetupSyst_EDivWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EDiv(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par/*,
      //const std::map<std::string, ThinSliceSystematic> & pars);
      int ediv_selection_ID*/);
  static double GetSystWeight_EDivNoPar(
      const ThinSliceEvent & event, int signal_index);

  void SetupSyst_EndZNoTrackWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EndZNoTrack(
      const ThinSliceEvent & event,
      int signal_index,
      const ThinSliceSystematic & par/*,
      //const std::map<std::string, ThinSliceSystematic> & pars,
      int upstream_ID, int no_track_ID*/);

  static double GetSystWeight_EndZNoTrackNoPar(
      const ThinSliceEvent & event,
      int signal_index);

  void SetupSyst_BeamMatch(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamMatch(
      const ThinSliceEvent & event,
      int signal_index,
      const ThinSliceSystematic & par);

  static double GetSystWeight_BeamMatchNoPar(
      const ThinSliceEvent & event,
      int signal_index);

  void SetupSyst_BeamMatchLow(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamMatchLow(
      const ThinSliceEvent & event,
      int signal_index,
      const ThinSliceSystematic & par);
  static double GetSystWeight_BeamMatchLowNoPar(
      const ThinSliceEvent & event,
      int signal_index);

  void SetupSyst_BeamMatchHigh(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamMatchHigh(
      const ThinSliceEvent & event,
      int signal_index,
      const ThinSliceSystematic & par);
  static double GetSystWeight_BeamMatchHighNoPar(
      const ThinSliceEvent & event,
      int signal_index);


  void SetupSyst_BeamScraper(
      const std::map<std::string, ThinSliceSystematic> & pars);
  static double GetSystWeight_BeamScraper(
      const ThinSliceEvent & event,
      int signal_index);

  void SetupSyst_BeamScraperHighPMu(
      const std::map<std::string, ThinSliceSystematic> & pars);
  static double GetSystWeight_BeamScraperHighPMu(
      const ThinSliceEvent & event,
      int signal_index);
  void SetupSyst_BeamScraperLowPMu(
      const std::map<std::string, ThinSliceSystematic> & pars);
  static double GetSystWeight_BeamScraperLowPMu(
      const ThinSliceEvent & event,
      int signal_index);

  void SetupSyst_BeamScraperPi(
      const std::map<std::string, ThinSliceSystematic> & pars);
  static double GetSystWeight_BeamScraperPi(
      const ThinSliceEvent & event,
      int signal_index);

  double GetSystWeight_UpstreamInt(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par/*,
      //const std::map<std::string, ThinSliceSystematic> & pars,
      int upstream_ID*/);

  double GetSystWeight_BGPions(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par/*,
      //const std::map<std::string, ThinSliceSystematic> & pars,
      int past_FV_ID, int decay_ID*/);

  void SetupSyst_BoxBeam(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BoxBeam(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par/*,
      //const std::map<std::string, ThinSliceSystematic> & pars,
      int beam_cut_ID*/);

  void SetupSyst_TrueBeamShift(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_TrueBeamShift(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par);
      //const std::map<std::string, ThinSliceSystematic> & pars);

  void SetupSyst_QuadBeamShift(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_QuadBeamShift(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par);

  void SetupSyst_ELoss(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_ELoss(
      const ThinSliceEvent & event,
      int signal_index,
      const ThinSliceSystematic & par/*,
      const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID*/);

  void SetupSyst_ELossMuon(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_ELossMuon(
      const ThinSliceEvent & event,
      const ThinSliceSystematic & par/*,
      const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID*/);


  static double CheckAndReturn(double weight, std::string name,
                        const ThinSliceSystematic & par,
                        const ThinSliceEvent & event);
  //static double CheckAndReturn(double weight, std::string name,
  //                      ThinSliceSystematic * par,
  //                      const ThinSliceEvent & event);

  //EDiv
  static double fEDivF, fEDivCut;
  static const ThinSliceSystematic * fEDivPar;
  static int fPastFVSelectionID;

  //G4RW
  static std::vector<const ThinSliceSystematic *> fActiveG4RWPars;
  static std::vector<const ThinSliceSystematic *> fActiveTiedG4RWPars;
  static std::vector<const ThinSliceSystematic *> fActiveSignalG4RWPars;
  static std::vector<const ThinSliceSystematic *> fActiveTiedSignalG4RWPars;

  static int fUpstreamID, fNoTrackID;

  //EndZ
  static std::map<int, std::vector<double>> fEndZFractions;
  static double fEndZNoTrackCut;
  static const ThinSliceSystematic * fEndZPar;

  //BeamMatch
  static const ThinSliceSystematic * fBeamMatchPar;
  static std::map<int, std::vector<double>> fBeamMatchFractions;
  static double fBeamMatchFraction;
  static bool fBeamMatchUseSingleFrac;

  //BeamMatchLow
  static double fBeamMatchLowLimit, fBeamMatchLowFraction;
  static std::map<int, std::vector<double>> fBeamMatchLowFractions;
  static bool fBeamMatchLowUseSingleFrac;
  static const ThinSliceSystematic * fBeamMatchLowPar;

  //BeamMatchHigh
  static double fBeamMatchHighLimit, fBeamMatchHighFraction;
  static std::map<int, std::vector<double>> fBeamMatchHighFractions;
  static bool fBeamMatchHighUseSingleFrac;
  static const ThinSliceSystematic * fBeamMatchHighPar;

  static const ThinSliceSystematic * fBeamScraperPar;
  static const ThinSliceSystematic * fBeamScraperParLowPMu;
  static const ThinSliceSystematic * fBeamScraperParHighPMu;
  static const ThinSliceSystematic * fBeamScraperParPi;

  //CosThetaEff
  static std::vector<int> fCosThetaEffIDs;
  static std::vector<int> fCosThetaEffSels;
  static std::vector<int> fCosThetaEffPDGs;
  static std::vector<std::pair<double, double>> fCosThetaEffRanges;
  static std::vector<double> fCosThetaEffFracs;
  static std::vector<const ThinSliceSystematic *> fActiveCosThetaEffPars;
  //static double GetSystWeight_CosThetaEff(const ThinSliceEvent & event);
 private:

  static double GetFractionBySample(
    const std::map<int, std::vector<double>> & fractions, int sample_ID,
    int signal_index);

  //G4RW Coeff
  std::map<std::string, std::string> fG4RWCoeffBranches;

  //BeamShiftBins
  std::map<std::string, double> fBeamShiftBinMeans, fBeamShiftBinSigmas,
                                fBeamShiftBinRangeLows, fBeamShiftBinRangeHighs;

  //BeamShift
  std::pair<double, double> fSystBeamShiftLimits;
  TGraph * fSystBeamShiftMeans, * fSystBeamShiftWidths;
  bool fSetupSystBeamShift = false;
  double fSystBeamShiftWeightCap;

  //EDiv
  //static double fEDivF, fEDivCut;
  //static const ThinSliceSystematic * fEDivPar;

  //EndZ
  //std::map<int, std::vector<double>> fEndZFractions;
  //double fEndZNoTrackCut;

  //Beam Match
  //std::vector<double> fBeamMatchLimits;
  //std::map<int, std::vector<double>> fBeamMatchFractions;

  //Beam Match Low
  //double fBeamMatchLowLimit, fBeamMatchLowFraction;
  //std::map<int, std::vector<double>> fBeamMatchLowFractions;
  //bool fBeamMatchLowUseSingleFrac;

  //Beam Match High
  //double fBeamMatchHighLimit, fBeamMatchHighFraction;
  //std::map<int, std::vector<double>> fBeamMatchHighFractions;
  //bool fBeamMatchHighUseSingleFrac;

  //Box Beam
  std::vector<std::pair<double, double>> fBoxBeamRegions;
  double fBoxBeamFraction;

  //TrueBeamShift
  std::vector<double> fTrueBeamBins;
  std::vector<TSpline3 *> fTrueBeamSplines;

  //ELoss
  double fELossCut;
  //std::map<int, double> fELossFractions;
  std::map<int, std::vector<double>> fELossFractions;
  double fELossMuonCut;
  std::map<int, double> fELossMuonFractions;

  //std::vector<std::string> fActiveSysts;
  int /*fUpstreamID, fNoTrackID, */fDecayID, fPastFVID, fBeamCutID;
  //static int fPastFVSelectionID;

  //Quad Beam Shift
  std::vector<double> fQuadBeamBins, fQuadBeamMeans, fQuadBeamSigmas;


  std::vector<
      std::pair<std::string,
                std::function<double(const ThinSliceEvent & event,
                                     int signal_index)>>> fActiveSystematics;
};
}
#endif
