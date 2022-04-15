#ifndef ABSCEXDRIVER_hh
#define ABSCEXDRIVER_hh

#include "ThinSliceDriver.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSpline.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TRandom3.h"
#include <map>
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "PDSPSystematics.h"

namespace protoana {
class AbsCexDriver : public ThinSliceDriver {
 public:
  AbsCexDriver(const fhicl::ParameterSet & extra_options);
  virtual ~AbsCexDriver();

  void FillMCEvents(
    TTree * tree, std::vector<ThinSliceEvent> & events,
    std::vector<ThinSliceEvent> & fake_data_events,
    int & split_val, const bool & do_split, int max_entries,
    const bool & do_fake_data) override;

  void BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val = 0) override;
  void BuildFakeData(
    TTree * tree, const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    int split_val = 0) override;
  void FakeDataSampleScales(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataBinnedScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataG4RWGrid(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    int split_val = 0);
  void FakeDataEffVar(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);

  void FakeDataLowP(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales, int split_val = 0);

  void FakeDatadEdX(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);

  void FakeDataPionAngle(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);

  void FakeDataAngleVar(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataBeamWeight(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void BuildMCSamples(
      //TTree * tree,
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
      std::vector<double> & beam_energy_bins, bool use_beam_inst_P) override;

  void RefillMCSamples(
      //TTree * tree,
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<int, std::vector<double>> & signal_pars,
      const std::map<int, double> & flux_pars,
      const std::map<std::string, ThinSliceSystematic> & syst_pars,
      bool fit_under_over, bool tie_under_over,
      bool use_beam_inst_P,
      bool fill_incident = false) override;

  /*void BuildSystSamples(
      TTree * tree,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins) override;*/
  
  void SetupSyst_G4RWCoeff(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_G4RW(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);

  void SetupSyst_dEdX_Cal(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  void SetupSyst_BeamShiftSpline2(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);

  void SetupSyst_BeamShiftRatio(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);

  std::pair<double, size_t> CalculateChi2(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set) override;
  void CompareSelections(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned,
      bool post_fit, int nPars,
      TDirectory * plot_dir) override;

  void GetCurrentHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      std::map<int, std::vector<TH1*>> & throw_hists,
      bool plot_rebinned) override;

  virtual void GetCurrentTruthHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      std::map<int, std::vector<TH1*>> & hists,
      std::map<int, std::vector<TH1*>> & inc_hists,
      std::map<int, std::vector<TH1*>> & xsec_hists,
      const std::vector<int> & incident_samples,
      const std::map<int, std::vector<double>> & signal_bins) override;

  void PlotThrows(
      ThinSliceDataSet & data_set, std::map<int, std::vector<TH1*>> & throw_hists,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      size_t nThrows,
      std::map<int, std::vector<TH1*>> & truth_throw_hists,
      std::map<int, std::vector<TH1*>> & truth_inc_hists,
      std::map<int, std::vector<TH1*>> & truth_xsec_hists,
      std::map<int, TH1*> & best_fit_incs,
      std::map<int, TH1*> & best_fit_xsecs,
      std::map<int, TH1*> & nominal_incs,
      std::map<int, TH1*> & nominal_xsecs,
      TFile & output_file, bool plot_rebinned,
      std::map<int, std::vector<double>> * sample_scales = 0x0) override;

  void SetupSysts(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file) override;

  /*void SetupSyst_BeamRes(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);*/
  void SetupSyst_BeamShift(
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  void SetupSyst_BeamShiftSpline(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  void SetupSyst_EndZNoTrackWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_BeamMatch(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_LowP(
    const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_NPi0(
    const std::map<std::string, ThinSliceSystematic> & pars);
  /*void SetupSyst_BeamShift2D(
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);*/
  void SetupSyst_EffVar(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  void SetupSyst_EffVarWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_EDivWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_NoTrackWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_BeamEffsWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  void SetupSyst_BoxBeam(
      const std::map<std::string, ThinSliceSystematic> & pars);

  /*double GetSystWeight_BeamRes(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);*/
  double GetSystWeight_BeamShift(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  /*double GetSystWeight_BeamShift2D(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);*/
  double GetSystWeight_G4RW(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars,
      const ThinSliceSample & sample,
      int selection_ID, double val);
  double GetSystWeight_G4RWCoeff(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EffVar(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_LowP(
      const ThinSliceEvent & event,
      int signal_index,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_NPi0(
      const ThinSliceEvent & event,
      int signal_index,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EndZNoTrack(
      const ThinSliceEvent & event,
      int signal_index,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamMatch(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EDiv(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_NoTrack(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamEffs(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_UpstreamInt(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BoxBeam(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);
  void WrapUpSysts(TFile & output_file) override;
  int RecalculateSelectionID(
      const ThinSliceEvent & event,
      double C_cal,
      TProfile * prot_template);
  double TruncatedMean(const std::vector<double> & dEdX);
 private:
   TH2D * fEndSlices;
   TFile * fIn;
   std::map<int, double> fMeans;

   double fEnergyFix;
   bool fDoEnergyFix;

   double fPitch;
   double fZ0;
   bool fMultinomial;
   bool fSkipFirstLast;
   double fEndZCut;
   std::map<int, std::vector<double>> fEndZFractions;
   double fTrajZStart;
   std::string fSliceMethod;
   int fSliceCut;

   double fBetaP, fRho, fWion, fAlpha, fNominalCCal;

   //bool fStaticBeamResWidth = false;
   //bool fStaticBeamResMean = false;
   //double fBeamResMeanVal = 1.;
   //double fBeamResWidthVal = 1.;
   TTree /** fSystBeamResTree, */* fSystBeamShiftTree/*, * fSystBeamShift2DTree*/;
   //double fSystBeamResWeight, fSystBeamResMeanOutput, fSystBeamResWidthOutput;
   //double fSystBeamResWeightCap, fSystBeamResOutput;
   double fSystBeamShiftWeight, fSystBeamShiftVal, fSystBeamShiftR;
   bool /*fSetupSystBeamRes = false,*/ fSetupSystBeamShift = false,
        /*fSetupSystBeamShift2D = false, */fSetupSystEffVar = false,
        fSystBeamShiftTreeSave = false;
   //double fSystBeamShift2DWeight, fSystBeamShift2DBVal, fSystBeamShift2DVal,
   //       fSystBeamShift2DR;
  // double fEffVarSystVal;
   TGraph2D * fSystBeamShiftMap; // , * fSystBeam2DMeans, * fSystBeam2DStdDevs;
   TGraph * fSystBeamShiftMeans, * fSystBeamShiftWidths;
   //double fSystBeamShiftRatioLimitUp, fSystBeamShiftRatioLimitDown;
   std::pair<double, double> fSystBeamShiftLimits;
   double fSystBeamShiftWeightCap;

   std::map<std::string, std::map<int, std::vector<TH1D*>>> fFullSelectionVars;
   std::map<std::string, std::map<int, std::vector<TSpline3*>>> fFullSelectionSplines;

   std::map<std::string, std::map<int, std::vector<TH1D*>>> fG4RWSelectionVarsPlus;
   std::map<std::string, std::map<int, std::vector<TH1D*>>> fG4RWSelectionVarsMinus;
   std::vector<std::string> fActiveG4RWSysts;
   TRandom3 fRNG = TRandom3(0);

   double fEffVarF, fEffVarCut;
   std::vector<double> fLowPFractions, fNPi0Fractions;
   double fEndZNoTrackCut;
   double fEDivF, fEDivCut, fNoTrackF, fBeamCutF;
   ProtoDUNETrackUtils fTrackUtil;

   std::vector<double> MakeTrueIncidentEnergies(
     const std::vector<double> & true_beam_traj_Z,
     const std::vector<double> & true_beam_traj_KE,
     const std::vector<int> & true_beam_slices,
     const std::vector<double> & true_beam_incidentEnergies);

  
   int GetBeamBin(
     const std::vector<double> & beam_energy_bins,
     const double & true_beam_startP, bool restrict_P=false);

   TH1D fBeamShiftRatioNomHist;
   std::vector<TSpline3*> fBeamShiftRatioSplines;
   std::map<std::string, std::string> fG4RWCoeffBranches;
   std::vector<std::pair<double, double>> fBoxBeamRegions;
   double fBoxBeamFraction;

   std::vector<double> fBeamMatchLimits, fBeamMatchFractions;
   PDSPSystematics * fSystematics;

   bool fInclusive;
   std::vector<int> fERecoSelections, fEndZSelections, fOneBinSelections;
   double fBeamInstPScale;
   bool fRestrictBeamInstP, fDebugRestrictBeamP;
};
}
#endif
