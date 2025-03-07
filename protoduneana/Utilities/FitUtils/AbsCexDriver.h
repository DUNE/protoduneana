#ifndef ABSCEXDRIVER_hh
#define ABSCEXDRIVER_hh

#include "ThinSliceDriver.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSpline.h"
#include "TGraph2D.h"
#include "TGraph.h"
// #include "TRandom3.h"
#include <map>
#include "TF1.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "PDSPSystematics.h"
#include "TDecompChol.h"
#include <mutex>

namespace protoana {

struct ExtraHistDataVars {
  ExtraHistDataVars(int selid, double endz, double startx, double starty,
                    double startz, double michel, int nhits, double calo_endz)
    : selection_ID(selid), reco_endz(endz),
      reco_startx_sce(startx), reco_starty_sce(starty), reco_startz_sce(startz),
      vertex_michel_score(michel), vertex_nhits(nhits), reco_endz_sce(calo_endz)
  {};

  void AddRecoProductTrackScores(std::vector<double> & scores) {
    reco_product_track_scores = scores;
    //reco_product_track_scores.insert(reco_product_track_scores.begin(),
    //                                 scores.begin(),
    //                                 scores.end());
  };

  void AddRecoProductTrunc_dEdXs(std::vector<double> & vals) {
    reco_product_truncated_dEdXs = vals;
  };
  void AddRecoProductChi2sPerHit(std::vector<double> & vals) {
    reco_product_chi2s_per_hit = vals;
  };

  void AddRecoProductShowerEnergy(std::vector<double> & vals) {
    reco_daughter_shower_energy = vals;
  };

  int selection_ID;
  double reco_endz;
  double reco_startx_sce, reco_starty_sce, reco_startz_sce;
  double vertex_michel_score;
  int vertex_nhits;
  double reco_endz_sce;
  std::vector<double> reco_product_track_scores,
                      reco_product_truncated_dEdXs,
                      reco_product_chi2s_per_hit,
                      reco_daughter_shower_energy;
};

class AbsCexDriver : public ThinSliceDriver {
 public:
  AbsCexDriver(const fhicl::ParameterSet & extra_options);
  virtual ~AbsCexDriver();

  // void FillMCEvents(
  //   std::vector<ThinSliceEvent> & events,
  //   std::vector<ThinSliceEvent> & fake_data_events,
  //   int & split_val, const bool & do_split, const bool & shuffle,
  //   int max_entries, int max_fake_entries,
  //   const bool & do_fake_data) override;

  void BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val = 0) override;

  void FillDataHistsFromSamples(
    ThinSliceDataSet & data_set,
    const std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    double & flux, std::vector<double> & fluxes_by_beam, bool fluctuate
  ) override;

  void BuildFakeData(
    TTree * tree, const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val = 0, bool scale_to_data_beam_p = false) override;
  void FakeDataSampleScales(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val = 0, bool scale_to_data_beam_p = false);
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
    std::vector<double> & beam_fluxes,
    int split_val = 0, bool norm_to_data_beam_P = false);
  void FakeDataEffVar(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataMuonVar(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set, double & flux,
    std::vector<double> & beam_fluxes);

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
    //TTree * tree,
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0,
    bool norm_to_data_beam_P = false);
  void FakeDataBeamWeight(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    int split_val = 0);
  void FakeDataBeamScale(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins, std::vector<double> & beam_fluxes,
    std::map<int, std::vector<double>> & sample_scales, int split_val = 0,
    bool norm_to_data_beam_P = false);
  void BuildMCSamples(
      //TTree * tree,
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
      std::vector<double> & beam_energy_bins, bool use_beam_inst_P) override;

  void RefillSampleLoop(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<int, std::vector<double>> & signal_pars,
      const std::map<int, double> & flux_pars,
      const std::map<std::string, ThinSliceSystematic> & syst_pars,
      const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
      bool fit_under_over, bool tie_under_over, bool use_beam_inst_P,
      bool fill_incident, std::map<int, TH1 *> * fix_factors,
      size_t worker_id, std::vector<size_t> n_events);

  void RefillMCSamples(
      //TTree * tree,
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<int, std::vector<double>> & signal_pars,
      const std::map<int, double> & flux_pars,
      const std::map<std::string, ThinSliceSystematic> & syst_pars,
      const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
      bool fit_under_over, bool tie_under_over,
      bool use_beam_inst_P, bool fill_incident = false,
      std::map<int, TH1*> * fix_factors = 0x0) override;

  /*void BuildSystSamples(
      TTree * tree,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins) override;*/
  
/*
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
      TFile & output_file);*/

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
      const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
      TFile & output_file) override;

  /*void SetupSyst_BeamRes(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
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
  void SetupSyst_BeamShift2D(
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  
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
      const std::map<std::string, ThinSliceSystematic> & pars);*/

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
  //void WrapUpSysts(TFile & output_file) override;
  int RecalculateSelectionID(
      const ThinSliceEvent & event,
      double C_cal,
      TProfile * prot_template);
  double TruncatedMean(const std::vector<double> & dEdX);

  //  void ConstructCovariances(
  //   const std::vector<ThinSliceEvent> & events,
  //   std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
  //   std::map<int, std::vector<std::vector<ThinSliceSample>>> & covariance_samples,
  //   const std::map<int, bool> & signal_sample_checks,
  //   std::vector<double> & beam_energy_bins,
  //   std::map<int, double> & nominal_fluxes,
  //   std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
  //   const std::map<int, std::vector<double>> & signal_pars,
  //   const std::map<int, double> & flux_pars,
  //   const std::map<std::string, ThinSliceSystematic> & syst_pars,
  //   const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
  //   bool fit_under_over, bool tie_under_over, bool use_beam_inst_P

  //   /*   const std::vector<ThinSliceEvent> & events,
  //      std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
  //      std::map<int, std::vector<std::vector<ThinSliceSample>>> & covariance_samples,
  //      const std::map<int, bool> & signal_sample_checks,
  //      std::map<int, double> & nominal_fluxes,
  //      std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
  //      std::vector<double> & beam_energy_bins, bool use_beam_inst_P*/) override;
  void SetupExtraHists(
    ThinSliceDataSet & data_set, 
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & fake_samples) override;

  void TurnOnFakeData() override;
  void TurnOffFakeData() override;

 private:
   TH2D * fEndSlices;
   TFile * fIn;
   std::map<int, double> fMeans;

   double fEnergyFix;
   bool fDoEnergyFix, fDoEnergyByLen;

   double fPitch;
   double fZ0;
   bool fMultinomial;
   bool fSkipFirstLast;
   double fEndZCut;
   std::map<int, std::vector<double>> fEndZFractions;
  //  double fTrajZStart;
   std::string fSliceMethod;
   int fSliceCut;

   double fBetaP, fRho, fWion, fAlpha/*, fNominalCCal*/;

   //bool fStaticBeamResWidth = false;
   //bool fStaticBeamResMean = false;
   //double fBeamResMeanVal = 1.;
   //double fBeamResWidthVal = 1.;
   //double fSystBeamResWeight, fSystBeamResMeanOutput, fSystBeamResWidthOutput;
   //double fSystBeamResWeightCap, fSystBeamResOutput;
   //double fSystBeamShift2DWeight, fSystBeamShift2DBVal, fSystBeamShift2DVal,
   //       fSystBeamShift2DR;
  // double fEffVarSystVal;
   //double fSystBeamShiftRatioLimitUp, fSystBeamShiftRatioLimitDown;

   std::map<std::string, std::map<int, std::vector<TH1D*>>> fFullSelectionVars;
   std::map<std::string, std::map<int, std::vector<TSpline3*>>> fFullSelectionSplines;
   std::map<std::string, std::map<int, std::vector<TF1*>>> fFullSelectionFuncs;

   std::map<std::string, std::map<int, std::vector<TH1D*>>> fG4RWSelectionVarsPlus;
   std::map<std::string, std::map<int, std::vector<TH1D*>>> fG4RWSelectionVarsMinus;
   std::vector<std::string> fActiveG4RWSysts;
  //  TRandom3 fRNG = TRandom3(0);

   double fEffVarF, fEffVarCut;
   std::vector<double> fLowPFractions, fNPi0Fractions;
   double fEndZNoTrackCut;
   double fEDivF, fEDivCut, fNoTrackF, fBeamCutF;
   ProtoDUNETrackUtils fTrackUtil;

   std::vector<double> MakeTrueIncidentEnergies(
     const std::vector<double> & true_beam_traj_Z,
     const std::vector<double> & true_beam_traj_KE/*,
     const std::vector<int> & true_beam_slices,
     const std::vector<double> & true_beam_incidentEnergies*/);

  
   int GetBeamBin(
     const std::vector<double> & beam_energy_bins,
     const double & true_beam_startP, bool restrict_P=false);

   void CovarianceRoutineBeamShift(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & new_samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    const std::map<int, std::vector<double>> & signal_pars,
    const std::map<int, double> & flux_pars,
    const std::map<std::string, ThinSliceSystematic> & syst_pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    bool fit_under_over, bool tie_under_over, bool use_beam_inst_P);
    //   const std::vector<ThinSliceEvent> & events,
    //   std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
    //   std::map<int, std::vector<std::vector<ThinSliceSample>>> & new_samples,
    //   const std::map<int, bool> & signal_sample_checks,
    //   std::map<int, double> & nominal_fluxes,
    //   std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    //   std::vector<double> & beam_energy_bins, bool use_beam_inst_P);

   void OpenBeamShiftInput();
   void SetupBeamShiftCovRoutine(fhicl::ParameterSet & routine);
   //double GetBeamShiftDelta(const std::vector<double> & energies);
   double GetBeamShiftDelta(double energy);
   void GenerateBeamShiftUniverse();

   void ScaleSamples(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      double scale);

   void SetupFakeDataG4RW();

   static double GetFakeWeight_G4RWCoeff(
      const ThinSliceEvent & event/*,
      const std::vector<std::string> & branches,
      const std::vector<double> & vars*/);

   //void SetupExtraHists(ThinSliceDataSet & data_set);
   //void SetupExtraHistEndZ();
   static void FillExtraHistDataSCEEndZ(ThinSliceDataSet & data_set,
                                       const ExtraHistDataVars & vars);
   static void FillExtraHistDataSCEEndZAbs(ThinSliceDataSet & data_set,
                                       const ExtraHistDataVars & vars);
   static void FillExtraHistDataSCEEndZCex(ThinSliceDataSet & data_set,
                                       const ExtraHistDataVars & vars);
   static void FillExtraHistDataSCEEndZOther(ThinSliceDataSet & data_set,
                                       const ExtraHistDataVars & vars);
   static void FillExtraHistDataShowerEnergy(ThinSliceDataSet & data_set,
                                       const ExtraHistDataVars & vars);
   static void FillExtraHistDataEndZ(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataEndZGoodReco(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataStartX(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataStartY(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataStartZ(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataTrackScore(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataVertexMichel(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataTrunc_dEdX(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   static void FillExtraHistDataChi2PerHit(ThinSliceDataSet & data_set,
                                 const ExtraHistDataVars & vars);
   void FillExtraHistsData(ThinSliceDataSet & data_set,
                           const ExtraHistDataVars & vars);

   static void FillExtraHistMCSCEEndZ(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCSCEEndZAbs(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCSCEEndZCex(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCSCEEndZOther(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCShowerEnergy(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCEndZ(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCEndZGoodReco(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCStartX(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCStartY(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCStartZ(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCTrackScore(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCVertexMichel(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCTrunc_dEdX(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   static void FillExtraHistMCChi2PerHit(ThinSliceSample & sample,
                                   const ThinSliceEvent & event,
                                   double weight);
   void FillExtraHistsMC(ThinSliceSample & sample,
                         const ThinSliceEvent & event,
                         double weight);

   TH1D fBeamShiftRatioNomHist;
   std::vector<TSpline3*> fBeamShiftRatioSplines;
   std::map<std::string, std::string> fG4RWCoeffBranches;

   std::vector<double> fBeamMatchLimits, fBeamMatchFractions;
   PDSPSystematics * fSystematics = 0x0;
   PDSPSystematics * fG4RWPars = 0x0;

   bool fInclusive, fAltFluxVar/*, fNewFVSelection*/;
   std::vector<int> fERecoSelections, fEndZSelections, fOneBinSelections;
   double fBeamInstPScale/*, fMCBeamInstPShift*/;
   bool fRestrictBeamInstP, fDebugRestrictBeamP;
   bool fVaryDataCalibration, fVaryMCCalibration, fVaryMCCalSelection, fThrowCalibration,
        fVaryCalibrationFakeData;
   bool fMCFrontUpSelection, fMCFrontDownSelection,
        fMCBackUpSelection, fMCBackDownSelection;
   double fDataCalibrationFactor, fMCCalibrationFactor;
   bool fBarlowBeeston;
   std::vector<int> fToSkip;
   size_t fNWorkers;
   //bool fStatVar;

   std::vector<std::string> fCovarianceRoutines;
   bool fBeamShiftCovRoutineActive = false;
   bool fRandomFakeDataBeamShift, fRandomMCBeamShift, fUseBeamShift, fOpenedBeamShiftInput = false, fDoingCovCreate = false;
   std::string fBeamShiftCovOutput;
   size_t fNCovarianceGens;
   std::pair<double, double> fBeamShiftCovP0, fBeamShiftCovP1, fBeamShiftCovP2;
   double fBeamShiftNominalP0, fBeamShiftNominalP1, fBeamShiftNominalP2;
   TVectorD * fBeamShiftInputCentrals;
   TMatrixD * fBeamShiftInputCov;
   TDecompChol fBeamShiftInputChol;
   TMatrixD * fBeamShiftInputCovL;
   bool fBeamShiftCovUseInput;
   std::string fBeamShiftInputFileName;
   double fCurrentBeamShiftP0, fCurrentBeamShiftP1, fCurrentBeamShiftP2;
   std::mutex fRefillMutex, fFillMutex;
   int fNThreads;
   std::vector<fhicl::ParameterSet> fExtraHistSets;
   double fFakeResolution = -999.;
   std::vector<std::vector<double>> fStoredEnergies;
   bool fUseStoredEnergies = false;

   std::vector<std::pair<std::string, std::function<void(
       ThinSliceDataSet & data_set, const ExtraHistDataVars & vars)>>>
           fActiveExtraHistsData;
   std::vector<std::pair<std::string, std::function<void(
       ThinSliceSample & sample, const ThinSliceEvent & event, double weight)>>>
           fActiveExtraHistsMC;

   std::map<std::string, std::function<void(
       ThinSliceDataSet & data_set, const ExtraHistDataVars & vars)>>
           fStoredExtraHistsData = {
    {"EndZ", FillExtraHistDataEndZ},
    {"EndZGoodReco", FillExtraHistDataEndZGoodReco},
    {"StartX", FillExtraHistDataStartX},
    {"StartY", FillExtraHistDataStartY},
    {"StartZ", FillExtraHistDataStartZ},
    {"TrackScore", FillExtraHistDataTrackScore},
    {"Trunc_dEdX", FillExtraHistDataTrunc_dEdX},
    {"Chi2PerHit", FillExtraHistDataChi2PerHit}
  };
  std::map<std::string, std::function<void(
       ThinSliceSample & sample, const ThinSliceEvent & event, double weight)>>
           fStoredExtraHistsMC = {
    {"EndZ", FillExtraHistMCEndZ},
    {"EndZGoodReco", FillExtraHistMCEndZGoodReco},
    {"StartX", FillExtraHistMCStartX},
    {"StartY", FillExtraHistMCStartY},
    {"StartZ", FillExtraHistMCStartZ},
    {"TrackScore", FillExtraHistMCTrackScore},
    {"Trunc_dEdX", FillExtraHistMCTrunc_dEdX},
    {"Chi2PerHit", FillExtraHistMCChi2PerHit}
  };

  std::function<double(const ThinSliceEvent & event)> FakeDataWeight;
  static double DummyWeight(const ThinSliceEvent & event) {return 1.;};

  static double fStartXDataMean, fStartXMCMean, 
                fStartYDataMean, fStartYMCMean,
                fStartZDataMean, fStartZMCMean,
                fStartXDataSigma, fStartXMCSigma, 
                fStartYDataSigma, fStartYMCSigma,
                fStartZDataSigma, fStartZMCSigma;

  void SetupOutgoingVarSyst();
  double GetEndKE(const ThinSliceEvent & event);
  double GetOutgoingVarWeight(const ThinSliceEvent & event);
  double GetPartialOutgoingWeight(
    double leading_value, double ke, 
    const std::vector<TH1D*> & ratios, const std::vector<double> & limits);
  bool fOutgoingSystVaryMomentum, fOutgoingSystCombine, fOutgoingSystActive;
  std::map<std::pair<int, int>, std::vector<TH1D*>> fOutgoingAngleRatios,
                                                    fOutgoingMomentumRatios;
  std::map<int, std::vector<double>> fOutgoingSystLimits;
  std::vector<int> fOutgoingSystCheckPDGs;


  void SetupMultiplicitySyst();
  double GetMultiplicityWeight(const ThinSliceEvent & event);
  bool fMultiplicitySystActive;
  std::map<int, std::vector<TH1D*>> fMultiplicityRatios;
  int fMultiplicityPDG;
  std::vector<int> fMultiplicitySamples;
  std::map<int, std::vector<double>> fMultiplicityLimits;
  bool fPreScaleSelection;
  double fNoBeamPreScale, fNoBeamFrac;


  bool fSCERatioSystActive;
  std::map<int, std::vector<double>> fSCERatioLimits;
  void SetupSCERatioSyst();
  std::map<std::pair<int, int>, std::vector<TH1D*>> fSCESystRatios;
  double GetSCERatioWeight(const ThinSliceEvent & event, double val);
protected:
  virtual void SetupExtraBranches() override;
  virtual void FillExtraBranches() override;
};
}
#endif
