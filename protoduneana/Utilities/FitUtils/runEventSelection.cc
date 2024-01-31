#include <iostream>
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include <ROOT/RDataFrame.hxx>
#include "TInterpreter.h"

#include "SelectionDefinitions.h"
#include "TTree.h"

class testing {
 private:
  int num; 
 public:
  testing(int n) : num(n){}
  int operator()() const {
    //std::cout << num << std::endl;
    return num;
  }
};

auto DefineMC(ROOT::RDataFrame & frame, const fhicl::ParameterSet & pset) {
  testing testing1(1);
  testing testing2(2);
  std::string track_score_string = (
      pset.get<bool>("UseWeightedTrackScore") ?
      "reco_daughter_PFP_trackScore_collection_weight_by_charge" :
      "reco_daughter_PFP_trackScore_collection");
                                    
  bool use_FV_beamcut = pset.get<bool>("UseFVBeamCut", false);
  if (use_FV_beamcut) std::cout << "Using FV_beamcut" << std::endl;
  auto mc = frame.Define("testing1", testing(1)/*testing1*/)
           .Define("testing2", testing(2)/*testing2*/)
           .Define("beam_P_range",
                   beam_P_range(pset.get<double>("BeamPLow", 0.),
                                pset.get<double>("BeamPHigh", 1.e6)),
                   {"beam_inst_P"})
           .Define("beam_XY_cuts",
                   beam_XY_cuts(pset.get<double>("BeamXMean", 0.),
                                pset.get<double>("BeamYMean", 0.),
                                pset.get<double>("BeamXYRadius", 0.)),
                   {"beam_inst_X", "beam_inst_Y", "beam_inst_nTracks"})
           .Define("primary_isBeamType", isBeamType(pset.get<bool>("CheckCalo")),
                   {"reco_beam_type", "reco_beam_incidentEnergies"})
           .Define("primary_ends_inAPA3",
                   endAPA3(pset.get<double>("EndZHigh")),
                   {(use_FV_beamcut ? "reco_beam_calo_endZ" : "reco_beam_endZ")})
           .Define("shower_dists",
                   shower_dists(pset.get<double>("TrackScoreCut")),
                   {track_score_string,
                    "reco_daughter_allShower_startX",
                    "reco_daughter_allShower_startY",
                    "reco_daughter_allShower_startZ",
                    "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ"})
           .Define("has_shower_dist_energy",
                   has_shower_dist_energy(pset.get<double>("TrackScoreCut")),
                   {track_score_string,
                    "reco_daughter_allShower_startX",
                    "reco_daughter_allShower_startY",
                    "reco_daughter_allShower_startZ",
                    "reco_daughter_allShower_energy",
                    "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ"})
           .Define("reco_daughter_allTrack_truncLibo_dEdX_pos",
                   truncatedMean_pos(pset.get<double>("Limit")),
                   {"reco_daughter_allTrack_calibrated_dEdX_SCE"})
           .Define("has_noPion_daughter",
                   secondary_noPion(
                       pset.get<double>("TrackScoreCut"),
                       pset.get<double>("Chi2Cut"),
                       pset.get<double>("dEdXLow"),
                       pset.get<double>("dEdXMed"),
                       pset.get<double>("dEdXHigh")),
                   {track_score_string,
                    "reco_daughter_allTrack_ID", 
                    "reco_daughter_allTrack_truncLibo_dEdX_pos",
                    "reco_daughter_allTrack_Chi2_proton",
                    "reco_daughter_allTrack_Chi2_ndof"})
           .Define("new_interaction_topology",
                   new_interaction_topology(pset.get<double>("EndZLow"),
                                            pset.get<double>("EndZHigh"),
                                            pset.get<double>("Threshold"),
                                            pset.get<bool>("CexNPi0"),
                                            pset.get<bool>("SignalPastFV", true)),
                   {"true_beam_PDG",
                    "true_beam_endZ", "true_beam_endProcess", "true_daughter_nPi0",
                    "true_beam_daughter_PDG", "true_beam_daughter_startP"})
           .Define("inclusive_topology", inclusive_topology(),
                   {"new_interaction_topology"})
           .Define("beam_backtrack", backtrack_beam,
                   {"reco_beam_true_byHits_process", "reco_beam_true_byHits_matched",
                    "reco_beam_true_byHits_origin", "reco_beam_true_byHits_PDG"})
           .Define("daughter_categories", categorize_daughters,
                   {"true_beam_ID", "reco_daughter_PFP_true_byHits_origin",
                    "reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_PDG",
                    "reco_daughter_PFP_true_byHits_parID", "reco_daughter_PFP_true_byHits_parPDG",
                    "true_beam_daughter_ID",
                    "true_beam_grand_daughter_ID"})
           //.Define("leading_p_costheta", leading_p_costheta,
           //        {"true_beam_endPx", "true_beam_endPy", "true_beam_endPz",
           //         "true_beam_daughter_PDG", "true_beam_daughter_startPx", 
           //         "true_beam_daughter_startPy", "true_beam_daughter_startPz"})
           .Define("leading_p_costheta", leading_costheta(2212),
                   {"true_beam_endPx", "true_beam_endPy", "true_beam_endPz",
                    "true_beam_daughter_PDG", "true_beam_daughter_startPx", 
                    "true_beam_daughter_startPy", "true_beam_daughter_startPz"})
           .Define("leading_piplus_costheta", leading_costheta(211),
                   {"true_beam_endPx", "true_beam_endPy", "true_beam_endPz",
                    "true_beam_daughter_PDG", "true_beam_daughter_startPx", 
                    "true_beam_daughter_startPy", "true_beam_daughter_startPz"})
           .Define("leading_pi0_costheta", leading_costheta(111),
                   {"true_beam_endPx", "true_beam_endPy", "true_beam_endPz",
                    "true_beam_daughter_PDG", "true_beam_daughter_startPx", 
                    "true_beam_daughter_startPy", "true_beam_daughter_startPz"})
           .Define("leading_p_momentum", leading_momentum(2212),
                   {"true_beam_daughter_PDG", "true_beam_daughter_startP"})
           .Define("leading_piplus_momentum", leading_momentum(211),
                   {"true_beam_daughter_PDG", "true_beam_daughter_startP"})
           .Define("leading_pi0_momentum", leading_momentum(111),
                   {"true_beam_daughter_PDG", "true_beam_daughter_startP"})
           .Define("reco_beam_fixed_interactingEnergy",
                   fixed_interacting_energy(80.),
                   {"reco_beam_incidentEnergies", "reco_beam_interactingEnergy"})
           .Define("reco_beam_modified_interactingEnergy",
                   modified_interacting_energy(80.),
                   {"beam_inst_P", "reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
           .Define("reco_beam_modified2_interactingEnergy",
                   modified_interacting_energy(20.),
                   {"beam_inst_P", "reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
           .Define("daughter_PDGs_types", daughter_PDG_types,
                   {"reco_daughter_PFP_true_byHits_PDG"})
           .Define("FV_vals",
                   FV_values(pset.get<double>("EndZCutFV", 30.)),
                   {"reco_beam_calo_X", "reco_beam_calo_Y",
                    "reco_beam_calo_Z", "reco_beam_calo_endZ",
                    "beam_inst_X", "beam_inst_Y",
                    "beam_inst_dirX", "beam_inst_dirY", "beam_inst_dirZ"})
           .Define("passBeamCutFV",
                   beam_cut_FV(
                       pset.get<double>("EndZCutFV", 30.),
                       pset.get<double>("MCCosLow", .95),
                       pset.get<double>("MCMeanX"),
                       pset.get<double>("MCMeanY"),
                       pset.get<double>("MCMeanR"),
                       pset.get<double>("MCSigmaX"),
                       pset.get<double>("MCSigmaY"),
                       pset.get<double>("MCSigmaR")),
                   {/*"reco_beam_calo_X", "reco_beam_calo_Y",
                    "reco_beam_calo_Z",*/"FV_vals", "reco_beam_calo_endZ"/*,
                    "beam_inst_X", "beam_inst_Y", "beam_inst_Z",
                    "beam_inst_dirX", "beam_inst_dirY", "beam_inst_dirZ"*/});

  if(pset.get<bool>("UseBI")) {
    mc = mc.Define(
        "passBeamCut",
        beam_cut_BI(pset.get<double>("MCXLow"),
                    pset.get<double>("MCXHigh"),
                    pset.get<double>("MCYLow"),
                    pset.get<double>("MCYHigh"),
                    pset.get<double>("MCZLow"),
                    pset.get<double>("MCZHigh"),
                    pset.get<double>("MCCosLow")),
        {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
         "reco_beam_trackDirX", "reco_beam_trackDirY",
         "reco_beam_trackDirZ", "beam_inst_X", "beam_inst_Y",
         "beam_inst_dirX", "beam_inst_dirY", "beam_inst_dirZ"});
  }
  else {
    mc = mc .Define(
        "passBeamCut",
        beam_cut_TPC(pset.get<bool>("DoAngle"),
                     pset.get<double>("XYZCut"),
                     pset.get<double>("CosCut"),
                     pset.get<double>("MCMeanX"),
                     pset.get<double>("MCMeanY"),
                     pset.get<double>("MCMeanZ"),
                     pset.get<double>("MCSigmaX"),
                     pset.get<double>("MCSigmaY"),
                     pset.get<double>("MCSigmaZ"),
                     pset.get<double>("MCMeanThetaX"),
                     pset.get<double>("MCMeanThetaY"),
                     pset.get<double>("MCMeanThetaZ")),
        {"reco_beam_calo_startX", "reco_beam_calo_startY",
         "reco_beam_calo_startZ", "reco_beam_calo_endX",
         "reco_beam_calo_endY", "reco_beam_calo_endZ"});
  }
  mc = mc.Define("vertex_cut",
                 vertex_michel_cut(pset.get<double>("MichelCut")),
                 {"reco_beam_vertex_michel_score", "reco_beam_vertex_nHits"});
  mc = mc.Define("selection_ID", selection_ID(pset.get<bool>("DoMichel")),
                 {"primary_isBeamType", "primary_ends_inAPA3",
                  "has_noPion_daughter",
                  (use_FV_beamcut ? "passBeamCutFV" : "passBeamCut"),
                  "has_shower_dist_energy", "vertex_cut"}) 
         .Define("selection_ID_inclusive",
                 selection_ID_inclusive(pset.get<bool>("DoMichel")),
                 {"primary_isBeamType", "primary_ends_inAPA3",
                  "passBeamCut", "vertex_cut"});

  if (pset.get<bool>("DoFakeResAndSel", false)) {
    auto fake_res_sel = pset.get<fhicl::ParameterSet>("FakeResSelPars");
    double fake_res = fake_res_sel.get<double>("Resolution", 10.);
    double fake_sel = fake_res_sel.get<double>("FakeSel", .05);

    mc = mc.Redefine("reco_beam_interactingEnergy", fake_res_func(fake_res),
                     {"true_beam_PDG", "true_beam_endP"})
           .Redefine("selection_ID", fake_selection(fake_sel),
                     {"new_interaction_topology"});
  }

  std::cout << "Filtering MC" << std::endl;
  auto filtered = mc.Filter("true_beam_PDG == 211 || true_beam_PDG == -13");
  return filtered;
}

auto DefineData(ROOT::RDataFrame & frame, const fhicl::ParameterSet & pset) {
  testing testing1(3);
  testing testing2(4);
  std::string track_score_string = (
      pset.get<bool>("UseWeightedTrackScore") ?
      "reco_daughter_PFP_trackScore_collection_weight_by_charge" :
      "reco_daughter_PFP_trackScore_collection");

  bool use_FV_beamcut = pset.get<bool>("UseFVBeamCut", false);
  if (use_FV_beamcut) std::cout << "Using FV_beamcut" << std::endl;
  auto data = frame.Define("testing1", testing(1)/*testing1*/)
           .Define("testing2", testing(2)/*testing2*/)
           .Define("beamPID", data_beam_PID, {"beam_inst_PDG_candidates", "MC", "true_beam_PDG"})
           .Define("beam_P_range",
                   beam_P_range(pset.get<double>("BeamPLow", 0.),
                                pset.get<double>("BeamPHigh", 1.e6)),
                   {"beam_inst_P"})
           .Define("beam_XY_cuts",
                   beam_XY_cuts(pset.get<double>("DataBeamXMean", 0.),
                                pset.get<double>("DataBeamYMean", 0.),
                                pset.get<double>("BeamXYRadius", 0.)),
                   {"beam_inst_X", "beam_inst_Y", "beam_inst_nTracks"})
           .Define("passBeamQuality",
                   data_BI_quality(pset.get<bool>("DoNTracks")),
                   {"beam_inst_nMomenta", "beam_inst_nTracks"})
           .Define("primary_isBeamType", isBeamType(pset.get<bool>("CheckCalo")),
                   {"reco_beam_type", "reco_beam_incidentEnergies"})
           .Define("primary_ends_inAPA3",
                   endAPA3(pset.get<double>("EndZHigh")),
                   {(use_FV_beamcut ? "reco_beam_calo_endZ" : "reco_beam_endZ")})
           .Define("shower_dists",
                   shower_dists(pset.get<double>("TrackScoreCut")),
                   {track_score_string,
                    "reco_daughter_allShower_startX",
                    "reco_daughter_allShower_startY",
                    "reco_daughter_allShower_startZ",
                    "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ"})
           .Define("has_shower_dist_energy",
                   has_shower_dist_energy(pset.get<double>("TrackScoreCut")),
                   {track_score_string,
                    "reco_daughter_allShower_startX",
                    "reco_daughter_allShower_startY",
                    "reco_daughter_allShower_startZ",
                    "reco_daughter_allShower_energy",
                    "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ"})
           .Define("reco_daughter_allTrack_truncLibo_dEdX_pos",
                   truncatedMean_pos(pset.get<double>("Limit")),
                   {"reco_daughter_allTrack_calibrated_dEdX_SCE"})
           .Define("reco_beam_fixed_interactingEnergy",
                   fixed_interacting_energy(80.),
                   {"reco_beam_incidentEnergies", "reco_beam_interactingEnergy"})
           .Define("reco_beam_modified_interactingEnergy",
                   modified_interacting_energy(80.),
                   {"beam_inst_P", "reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
           .Define("reco_beam_modified2_interactingEnergy",
                   modified_interacting_energy(20.),
                   {"beam_inst_P", "reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
           .Define("has_noPion_daughter",
                   secondary_noPion(
                       pset.get<double>("TrackScoreCut"),
                       pset.get<double>("Chi2Cut"),
                       pset.get<double>("dEdXLow"),
                       pset.get<double>("dEdXMed"),
                       pset.get<double>("dEdXHigh")),
                   {track_score_string,
                    "reco_daughter_allTrack_ID", 
                    "reco_daughter_allTrack_truncLibo_dEdX_pos",
                    "reco_daughter_allTrack_Chi2_proton",
                    "reco_daughter_allTrack_Chi2_ndof"})
           .Define("good_run",
                   exclude_runs(
                       pset.get<std::vector<int>>("ExcludeRuns")),
                   {"run"})
           .Define("FV_vals",
                   FV_values(pset.get<double>("EndZCutFV", 30.)),
                   {"reco_beam_calo_X", "reco_beam_calo_Y",
                    "reco_beam_calo_Z", "reco_beam_calo_endZ",
                    "beam_inst_X", "beam_inst_Y",
                    "beam_inst_dirX", "beam_inst_dirY", "beam_inst_dirZ"})
           .Define("passBeamCutFV",
                   beam_cut_FV(
                       pset.get<double>("EndZCutFV", 30.),
                       pset.get<double>("DataCosLow", .95),
                       pset.get<double>("DataMeanX"),
                       pset.get<double>("DataMeanY"),
                       pset.get<double>("DataMeanR"),
                       pset.get<double>("DataSigmaX"),
                       pset.get<double>("DataSigmaY"),
                       pset.get<double>("DataSigmaR")),
                   {/*"reco_beam_calo_X", "reco_beam_calo_Y",
                    "reco_beam_calo_Z",*/ "FV_vals", "reco_beam_calo_endZ"/*,
                    "beam_inst_X", "beam_inst_Y", "beam_inst_Z",
                    "beam_inst_dirX", "beam_inst_dirY", "beam_inst_dirZ"*/});
           /*.Define("beam_inst_P_scaled", beam_inst_P_scaled())*/;

  if(pset.get<bool>("UseBI")) {
    data = data.Define(
        "passBeamCut",
        beam_cut_BI(pset.get<double>("DataXLow"),
                    pset.get<double>("DataXHigh"),
                    pset.get<double>("DataYLow"),
                    pset.get<double>("DataYHigh"),
                    pset.get<double>("DataZLow"),
                    pset.get<double>("DataZHigh"),
                    pset.get<double>("DataCosLow")),
        {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
         "reco_beam_trackDirX", "reco_beam_trackDirY",
         "reco_beam_trackDirZ", "beam_inst_X", "beam_inst_Y",
         "beam_inst_dirX", "beam_inst_dirY", "beam_inst_dirZ"});
  }
  else {
    data = data.Define(
        "passBeamCut",
        beam_cut_TPC(pset.get<bool>("DoAngle"),
                     pset.get<double>("XYZCut"),
                     pset.get<double>("CosCut"),
                     pset.get<double>("DataMeanX"),
                     pset.get<double>("DataMeanY"),
                     pset.get<double>("DataMeanZ"),
                     pset.get<double>("DataSigmaX"),
                     pset.get<double>("DataSigmaY"),
                     pset.get<double>("DataSigmaZ"),
                     pset.get<double>("DataMeanThetaX"),
                     pset.get<double>("DataMeanThetaY"),
                     pset.get<double>("DataMeanThetaZ")),
        {"reco_beam_calo_startX", "reco_beam_calo_startY",
         "reco_beam_calo_startZ", "reco_beam_calo_endX",
         "reco_beam_calo_endY", "reco_beam_calo_endZ"});
  }


  data = data.Define("vertex_cut",
                     vertex_michel_cut(pset.get<double>("MichelCut")),
                     {"reco_beam_vertex_michel_score",
                      "reco_beam_vertex_nHits"})/*;
  data = data*/.Define("selection_ID", selection_ID(pset.get<bool>("DoMichel")),
                 {"primary_isBeamType", "primary_ends_inAPA3",
                  "has_noPion_daughter",
                  (use_FV_beamcut ? "passBeamCutFV" : "passBeamCut"),
                  "has_shower_dist_energy", "vertex_cut"})
                .Define("selection_ID_inclusive",
                        selection_ID_inclusive(pset.get<bool>("DoMichel")),
                        {"primary_isBeamType", "primary_ends_inAPA3",
                         "passBeamCut", "vertex_cut"});
  auto filtered = data.Filter("beamPID == true");
  return filtered;
}

int main(int argc, char ** argv){
  //gInterpreter->GenerateDictionary("vector<vector<int>>",
  //                                 "vector");

  std::string fcl_file;
  std::string output_file;
  std::string mc_file, data_file;
  bool found_mc = false, found_data = false;
  // Options to run
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-m")) {
     mc_file = argv[++iArg];
     found_mc = true;
    }
    if (!strcasecmp(argv[iArg],"-d")) {
     data_file = argv[++iArg];
     found_data = true;
    }
    if (!strcasecmp(argv[iArg],"-o")) {
      output_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: runPDSPThinSliceFit -c fclfile.fcl " << 
                    "-o outputfile.root " << std::endl;
      return 1;
    }
  }

  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;
  if (fhicl_env == nullptr) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  auto const pset = fhicl::ParameterSet::make(fcl_file, lookupPolicy);
  std::string tree_name = pset.get<std::string>("TreeName"); 

  if (pset.get<bool>("UseMT"))
    ROOT::EnableImplicitMT(4);

  const auto & column_list
      = pset.get<std::vector<std::string>>("ColumnList");
  if (found_mc) {
    ROOT::RDataFrame frame(tree_name, mc_file);
    std::cout << "Calling DefineMC" << std::endl;
    auto mc = DefineMC(frame, pset);
    std::cout << "Snapshotting" << std::endl;
    if (pset.get<bool>("DoReconstructable"))
      mc = mc.Filter("reco_reconstructable_beam_event");
    if (pset.get<bool>("RestrictBeamP")) {
      mc = mc.Filter("beam_P_range");
    }

    if (pset.get<bool>("RestrictBeamXY", false)) {
      mc = mc.Filter("beam_XY_cuts");
    }

    auto time0 = std::chrono::high_resolution_clock::now();
    mc.Snapshot(tree_name, "eventSelection_mc_all.root", column_list);
    //mc.Snapshot(tree_name, "eventSelection_mc_reconstructable.root");
    auto time1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time: " <<
                 std::chrono::duration_cast<std::chrono::seconds>(time1 - time0).count() <<
                 std::endl;
  
  }

  if (found_data) {
    ROOT::RDataFrame data_frame(tree_name, data_file);
    std::cout << "Calling DefineData" << std::endl;
    auto data = DefineData(data_frame, pset);
    std::cout << "Snapshotting" << std::endl;
    //data.Filter("good_run");
    if (pset.get<bool>("SaveAllData", false))
      data.Snapshot(tree_name, "eventSelection_data_all.root");
    data = data.Filter("passBeamQuality && good_run");
    if (pset.get<bool>("DoReconstructable"))
      data = data.Filter("reco_reconstructable_beam_event");
    if (pset.get<bool>("RestrictBeamP")) {
      data = data.Filter("beam_P_range");
    }
    if (pset.get<bool>("RestrictBeamXY", false)) {
      data = data.Filter("beam_XY_cuts");
    }

    /*
    if (pset.get<bool>("RestrictBeamXY", false)) {
    }*/

    auto time0 = std::chrono::high_resolution_clock::now();
    data.Snapshot(tree_name, "eventSelection_data_BeamQuality.root", column_list);
    auto time1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time: " <<
                 std::chrono::duration_cast<std::chrono::seconds>(time1 - time0).count() <<
                 std::endl;
  }

  return 0;
}
