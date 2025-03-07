#include "ThinSliceDriver.h"

#include "THStack.h"
#include "TCanvas.h"

protoana::ThinSliceDriver::ThinSliceDriver(
    const fhicl::ParameterSet & extra_options)
    : fExtraOptions(extra_options),
      fNewFVSelection(extra_options.get<bool>("NewFVSelection", false)),
      fTrajZStart(extra_options.get<double>("TrajZStart", 0.)) {
  fFakeDataRoutine = fExtraOptions.get<std::string>("FakeDataRoutine", "");
}

protoana::ThinSliceDriver::~ThinSliceDriver() {}

void protoana::ThinSliceDriver::CompareDataMC(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    int nPars,
    TDirectory * plot_dir,
    bool plot_rebinned,
    bool post_fit) {

  for (auto it = samples.begin(); it != samples.end(); ++it) {
    std::vector<ThinSliceSample> & vec = it->second[0];
    for (size_t i = 0; i < vec.size(); ++i) {
      vec[i].RefillRebinnedHists();

    }

    for (size_t i = 1; i < it->second.size(); ++i) {
      std::vector<ThinSliceSample> & vec = it->second[i];
      for (size_t j = 0; j < vec.size(); ++j) {
        vec[j].RefillRebinnedHists();
      }
    }
  }

  CompareSelections(
      samples, data_set, output_file, plot_style, plot_rebinned,
      post_fit, nPars, plot_dir);
}

void protoana::ThinSliceDriver::FillMCEvents(
  std::vector<ThinSliceEvent> & events,
  std::vector<ThinSliceEvent> & fake_data_events,
  int & split_val, const bool & do_split, const bool & shuffle,
  int max_entries, int max_fake_entries, const bool & do_fake_data) {
std::cout << "Filling MC Events" << std::endl;

int sample_ID, selection_ID, event, run, subrun;
int cal_up_selection_ID, cal_down_selection_ID;
int selection_ID_front_shift_up, selection_ID_front_shift_down,
    selection_ID_back_shift_up, selection_ID_back_shift_down;
int true_beam_PDG;
double true_beam_interactingEnergy, reco_beam_interactingEnergy, reco_beam_alt_len;
double true_beam_endP, true_beam_mass, true_beam_endZ;
double reco_beam_endZ, true_beam_startP, reco_beam_startY;
double reco_beam_startX_SCE, reco_beam_startY_SCE, reco_beam_startZ_SCE;
double reco_beam_endZ_SCE;
double vertex_michel_score;
int vertex_nhits;
double beam_inst_P;
double leading_p_costheta, leading_piplus_costheta, leading_pi0_costheta;
double leading_p_momentum, leading_piplus_momentum, leading_pi0_momentum;
std::vector<double> * reco_beam_incidentEnergies = 0x0,
                    * true_beam_incidentEnergies = 0x0,
                    * true_beam_traj_Z = 0x0,
                    * true_beam_traj_KE = 0x0,
                    * reco_daughter_track_thetas = 0x0,
                    * reco_daughter_track_scores = 0x0;
fMCTree->SetBranchAddress("event", &event); //good
fMCTree->SetBranchAddress("subrun", &subrun); //good
fMCTree->SetBranchAddress("run", &run); //good

fMCTree->SetBranchAddress("true_beam_PDG", &true_beam_PDG); //good

// if (!fInclusive) { //good
  fMCTree->SetBranchAddress("new_interaction_topology", &sample_ID);

  if (fNewFVSelection) {
    fMCTree->SetBranchAddress("selection_ID_FV", &selection_ID);
  }
  else {
    fMCTree->SetBranchAddress("selection_ID", &selection_ID);
  }
// }
// else {
//   fMCTree->SetBranchAddress("inclusive_topology", &sample_ID);
//   fMCTree->SetBranchAddress("selection_ID_inclusive", &selection_ID);
// }
fMCTree->SetBranchAddress("cal_up_selection_ID", &cal_up_selection_ID);
fMCTree->SetBranchAddress("cal_down_selection_ID", &cal_down_selection_ID);

fMCTree->SetBranchAddress("selection_ID_front_shift_up",
                       &selection_ID_front_shift_up);
fMCTree->SetBranchAddress("selection_ID_front_shift_down",
                       &selection_ID_front_shift_down);

fMCTree->SetBranchAddress("selection_ID_back_shift_up",
                       &selection_ID_back_shift_up);
fMCTree->SetBranchAddress("selection_ID_back_shift_down",
                       &selection_ID_back_shift_down);

fMCTree->SetBranchAddress("true_beam_interactingEnergy", //good
                       &true_beam_interactingEnergy);
fMCTree->SetBranchAddress("true_beam_endP", &true_beam_endP); //good
fMCTree->SetBranchAddress("true_beam_endZ", &true_beam_endZ); //good
fMCTree->SetBranchAddress("true_beam_mass", &true_beam_mass); //good
fMCTree->SetBranchAddress("reco_beam_interactingEnergy", //good
                       &reco_beam_interactingEnergy);
fMCTree->SetBranchAddress("reco_beam_alt_len", &reco_beam_alt_len);
fMCTree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ); //good
fMCTree->SetBranchAddress("reco_beam_calo_startX", &reco_beam_startX_SCE);
fMCTree->SetBranchAddress("reco_beam_calo_startY", &reco_beam_startY_SCE);
fMCTree->SetBranchAddress("reco_beam_calo_startZ", &reco_beam_startZ_SCE);
fMCTree->SetBranchAddress("reco_beam_calo_endZ", &reco_beam_endZ_SCE);
fMCTree->SetBranchAddress("reco_beam_vertex_michel_score", &vertex_michel_score);
fMCTree->SetBranchAddress("reco_beam_vertex_nHits", &vertex_nhits);
fMCTree->SetBranchAddress("reco_beam_startY", &reco_beam_startY); //good
fMCTree->SetBranchAddress("reco_beam_incidentEnergies", //good
                       &reco_beam_incidentEnergies);
fMCTree->SetBranchAddress("true_beam_incidentEnergies", //good
                       &true_beam_incidentEnergies);

fMCTree->SetBranchAddress("true_beam_startP", &true_beam_startP); //good
fMCTree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z); //good
fMCTree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE); //good

fMCTree->SetBranchAddress("beam_inst_P", &beam_inst_P); //good
fMCTree->SetBranchAddress("reco_daughter_allTrack_Theta", &reco_daughter_track_thetas); //good
fMCTree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", //good
                          &reco_daughter_track_scores);
std::vector<int> * true_beam_daughter_PDG = 0x0;
fMCTree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);//good

int true_beam_ID, reco_beam_true_byHits_ID;
fMCTree->SetBranchAddress("true_beam_ID", &true_beam_ID); //good
fMCTree->SetBranchAddress("reco_beam_true_byHits_ID", &reco_beam_true_byHits_ID); //good

std::vector<std::vector<double>> * g4rw_full_grid_proton_coeffs = 0x0,
                                 * g4rw_full_grid_piplus_coeffs = 0x0,
                                 * g4rw_full_fine_piplus_coeffs = 0x0,
                                 * g4rw_full_grid_abscex_coeffs = 0x0,
                                 * g4rw_primary_grid_abscex_coeffs = 0x0,
                                 * g4rw_downstream_grid_piplus_coeffs = 0x0;
fMCTree->SetBranchAddress("g4rw_full_grid_proton_coeffs", //good
                       &g4rw_full_grid_proton_coeffs);
fMCTree->SetBranchAddress("g4rw_full_grid_piplus_coeffs",
                       &g4rw_full_grid_piplus_coeffs);
fMCTree->SetBranchAddress("g4rw_full_fine_piplus_coeffs",
                       &g4rw_full_fine_piplus_coeffs);
fMCTree->SetBranchAddress("g4rw_full_grid_abscex_coeffs",
                       &g4rw_full_grid_abscex_coeffs);
fMCTree->SetBranchAddress("g4rw_primary_grid_abscex_coeffs",
                       &g4rw_primary_grid_abscex_coeffs);
fMCTree->SetBranchAddress("g4rw_downstream_grid_piplus_coeffs", //good
                       &g4rw_downstream_grid_piplus_coeffs);

std::vector<std::vector<double>> * daughter_dQdXs = 0x0,
                                 * daughter_resRanges = 0x0,
                                 * daughter_EFields = 0x0;
fMCTree->SetBranchAddress(
    "reco_daughter_allTrack_calibrated_dQdX_SCE", &daughter_dQdXs);//good
fMCTree->SetBranchAddress(
    "reco_daughter_allTrack_resRange_SCE", &daughter_resRanges);//good
fMCTree->SetBranchAddress(
    "reco_daughter_allTrack_EField_SCE", &daughter_EFields);//good
bool has_pi0_shower;
fMCTree->SetBranchAddress("has_shower_dist_energy", &has_pi0_shower);//good 
std::vector<double> * true_beam_daughter_startP = 0x0;
fMCTree->SetBranchAddress("true_beam_daughter_startP",  &true_beam_daughter_startP);//good
fMCTree->SetBranchAddress("leading_p_costheta", &leading_p_costheta);//good
fMCTree->SetBranchAddress("leading_piplus_costheta", &leading_piplus_costheta);//good
fMCTree->SetBranchAddress("leading_pi0_costheta", &leading_pi0_costheta);//good
fMCTree->SetBranchAddress("leading_p_momentum", &leading_p_momentum);//good
fMCTree->SetBranchAddress("leading_piplus_momentum", &leading_piplus_momentum);//good
fMCTree->SetBranchAddress("leading_pi0_momentum", &leading_pi0_momentum);//good
bool is_beam_scraper;
fMCTree->SetBranchAddress("true_beam_is_scraper", &is_beam_scraper);//good

std::vector<double> * reco_daughter_chi2 = 0x0,
                    * reco_daughter_truncated_dEdX = 0x0;
std::vector<int> * reco_daughter_chi2_nhits = 0x0;
fMCTree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton",
                       &reco_daughter_chi2);
fMCTree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof",
                       &reco_daughter_chi2_nhits);
fMCTree->SetBranchAddress("reco_daughter_allTrack_truncLibo_dEdX_pos",
                       &reco_daughter_truncated_dEdX);
int true_n_neutrons, true_n_protons, true_n_piplus, true_n_piminus,
    true_n_pi0;
fMCTree->SetBranchAddress("true_daughter_nNeutron", &true_n_neutrons);
fMCTree->SetBranchAddress("true_daughter_nProton", &true_n_protons);
fMCTree->SetBranchAddress("true_daughter_nPiPlus", &true_n_piplus);
fMCTree->SetBranchAddress("true_daughter_nPiMinus", &true_n_piminus);
fMCTree->SetBranchAddress("true_daughter_nPi0", &true_n_pi0);

std::vector<double> * reco_daughter_shower_energy = 0x0;
fMCTree->SetBranchAddress("reco_daughter_allShower_energy",
                       &reco_daughter_shower_energy);

int nentries = (max_entries < 0 ? fMCTree->GetEntries() : max_entries);
if (max_entries > fMCTree->GetEntries()) {
    std::string message = "Requested more entries than in MC fMCTree";
    throw std::runtime_error(message);
}

split_val = nentries;

int events_end = nentries;
int fake_start = 0;

std::vector<int> event_list, fake_event_list;

if (do_split) {
  split_val = nentries/2;
  events_end = nentries/2;
  fake_start = events_end;
  std::cout << "Note: Splitting MC in half. " <<
               split_val << "/" << nentries <<std::endl;
  if (shuffle) {
    for (int i = 0; i < nentries; ++i) {
      event_list.push_back(i);
    }
    for (int i = 0; i < events_end; ++i) {
      int r = fRNG.Integer(event_list.size());
      fake_event_list.push_back(event_list[r]);
      event_list.erase(event_list.begin()+r);
    }

    std::sort(event_list.begin(), event_list.end());
    std::sort(fake_event_list.begin(), fake_event_list.end());

    std::cout << "Shuffled " << event_list.size() << " " <<
                 fake_event_list.size() << std::endl;
    for (int & i : event_list) {
      if (std::find(fake_event_list.begin(), fake_event_list.end(), i) !=
          fake_event_list.end()) {
        std::cout << "Incorrectly shuffled" << std::endl;
      }
    }
  }
  else {
    for (int i = 0; i < events_end; ++i) {
      event_list.push_back(i);
    }
    for (int i = fake_start; i < nentries; ++i) {
      fake_event_list.push_back(i);
    }
    std::cout << "Split " << event_list.size() << " " <<
                 fake_event_list.size() << std::endl;
  }
}
else {
  int nfake_entries = (max_fake_entries < 0 ? fMCTree->GetEntries() : max_fake_entries);
  if (max_fake_entries > fMCTree->GetEntries()) {
      std::string message = "Requested more fake_entries than in MC fMCTree";
      throw std::runtime_error(message);
  }
  for (int i = 0; i < events_end; ++i) {
    event_list.push_back(i);
  }
  for (int i = 0; i < nfake_entries; ++i) {
    fake_event_list.push_back(i);
  }
}

/*
TEntryList event_list, fake_event_list;
for (int i = events_start; i < events_end; ++i) {
  event_list.Enter(i);
}
for (int i = fake_start; i < fake_end; ++i) {
  fake_event_list.Enter(i);
}

ROOT::EnableImplicitMT(fNThreads);
ROOT::TTreeProcessorMT tp(*fMCTree, event_list);

auto fill_func = [&](TTreeReader & myReader) {
  TTreeReaderValue<int> runRV(myReader, "run");
  TTreeReaderValue<int> eventRV(myReader, "event");
  TTreeReaderValue<int> sample_IDRV(myReader, (!fInclusive ? 
                                               "new_interaction_topology" :
                                               "inclusive_topology"));

  TTreeReaderValue<int> selection_IDRV(myReader, (!fInclusive ? 
                                                  "selection_ID" :
                                                  "selection_ID_inclusive"));
  TTreeReaderValue<int> subrunRV(myReader, "subrun");
  TTreeReaderValue<int> true_beam_PDGRV(myReader, "true_beam_PDG");

  TTreeReaderValue<double> true_beam_interactingEnergyRV(myReader, "true_beam_interactingEnergy");
  TTreeReaderValue<double> reco_beam_interactingEnergyRV(myReader, "reco_beam_interactingEnergy");
  TTreeReaderValue<double> true_beam_endPRV(myReader, "true_beam_endP");
  TTreeReaderValue<double> true_beam_massRV(myReader, "true_beam_mass");
  TTreeReaderValue<double> true_beam_endZRV(myReader, "true_beam_endZ");
  TTreeReaderValue<double> reco_beam_endZRV(myReader, "reco_beam_endZ");
  TTreeReaderValue<double> true_beam_startPRV(myReader, "true_beam_startP");
  TTreeReaderValue<double> reco_beam_startYRV(myReader, "reco_beam_startY");
  TTreeReaderValue<double> beam_inst_PRV(myReader, "beam_inst_P");
  TTreeReaderValue<double> leading_p_costhetaRV(myReader, "leading_p_costheta");
  TTreeReaderValue<double> leading_piplus_costhetaRV(myReader, "leading_piplus_costheta");
  TTreeReaderValue<double> leading_pi0_costhetaRV(myReader, "leading_pi0_costheta");

  TTreeReaderValue<std::vector<double>> reco_beam_incidentEnergiesRV(myReader, "reco_beam_incidentEnergies");
  TTreeReaderValue<std::vector<double>> true_beam_incidentEnergiesRV(myReader, "true_beam_incidentEnergies");
  TTreeReaderValue<std::vector<double>> true_beam_traj_ZRV(myReader, "true_beam_traj_Z");
  TTreeReaderValue<std::vector<double>> true_beam_traj_KERV(myReader, "true_beam_traj_KE");
  TTreeReaderValue<std::vector<double>> reco_daughter_track_thetasRV(myReader, "reco_daughter_allTrack_Theta");
  TTreeReaderValue<std::vector<double>> reco_daughter_track_scoresRV(myReader, "reco_daughter_PFP_trackScore_collection");
  TTreeReaderValue<std::vector<int>> true_beam_slicesRV(myReader, "true_beam_slices");
  TTreeReaderValue<std::vector<double>> calibrated_dQdXRV(myReader, "reco_beam_calibrated_dQdX_SCE");
  TTreeReaderValue<std::vector<double>> beam_EFieldRV(myReader, "reco_beam_EField_SCE");
  TTreeReaderValue<std::vector<double>> track_pitchRV(myReader, "reco_beam_TrkPitch_SCE");
  TTreeReaderValue<std::vector<int>> true_beam_daughter_PDGRV(myReader, "true_beam_daughter_PDG");


  TTreeReaderValue<std::vector<std::vector<double>>> g4rw_full_grid_proton_coeffsRV(myReader, "g4rw_full_grid_proton_coeffs");
  TTreeReaderValue<std::vector<std::vector<double>>> g4rw_downstream_grid_piplus_coeffsRV(myReader, "g4rw_downstream_grid_piplus_coeffs");
  TTreeReaderValue<std::vector<std::vector<double>>> daughter_dQdXsRV(myReader, "reco_daughter_allTrack_calibrated_dQdX_SCE");
  TTreeReaderValue<std::vector<std::vector<double>>> daughter_resRangesRV(myReader, "reco_daughter_allTrack_resRange_SCE");
  TTreeReaderValue<std::vector<std::vector<double>>> daughter_EFieldsRV(myReader, "reco_daughter_allTrack_EField_SCE");

  TTreeReaderValue<std::vector<double>> true_beam_daughter_startPRV(myReader, "true_beam_daughter_startP");
  TTreeReaderValue<int> true_beam_IDRV(myReader, "true_beam_ID");
  TTreeReaderValue<int> reco_beam_true_byHits_IDRV(myReader, "reco_beam_true_byHits_ID");
  TTreeReaderValue<bool> has_pi0_showerRV(myReader, "has_shower_dist_energy");
  TTreeReaderValue<bool> is_beam_scraperRV(myReader, "true_beam_is_scraper");


  std::cout << "Starting thread" << std::endl;
  while (myReader.Next()) {
    auto run = *runRV;
    auto subrun = *subrunRV;
    auto event = *eventRV;

    ThinSliceEvent thin_slice_event(event, subrun, run);
    thin_slice_event.SetSampleID(*sample_IDRV);
    thin_slice_event.SetSelectionID(*selection_IDRV);
    thin_slice_event.SetTrueInteractingEnergy(*true_beam_interactingEnergyRV);
    thin_slice_event.SetRecoInteractingEnergy(*reco_beam_interactingEnergyRV);
    thin_slice_event.SetTrueEndP(*true_beam_endPRV);
    thin_slice_event.SetTrueEndZ(*true_beam_endZRV);
    thin_slice_event.SetTrueStartP(*true_beam_startPRV);
    thin_slice_event.SetTrueMass(*true_beam_massRV);
    thin_slice_event.SetRecoEndZ(*reco_beam_endZRV);
    thin_slice_event.SetRecoStartY(*reco_beam_startYRV);

    thin_slice_event.SetTrueID(*true_beam_IDRV);
    thin_slice_event.SetRecoToTrueID(*reco_beam_true_byHits_IDRV);

    thin_slice_event.SetRecoIncidentEnergies(*reco_beam_incidentEnergiesRV);
    thin_slice_event.SetTrueIncidentEnergies(*true_beam_incidentEnergiesRV);
    thin_slice_event.SetTrueTrajZ(*true_beam_traj_ZRV);
    thin_slice_event.SetTrueTrajKE(*true_beam_traj_KERV);
    thin_slice_event.SetTrueSlices(*true_beam_slicesRV);
    thin_slice_event.SetdQdXCalibrated(*calibrated_dQdXRV);
    thin_slice_event.SetEField(*beam_EFieldRV);
    thin_slice_event.SetTrackPitch(*track_pitchRV);
    thin_slice_event.SetBeamInstP(*beam_inst_PRV);
    thin_slice_event.SetPDG(*true_beam_PDGRV);
    thin_slice_event.SetRecoDaughterTrackThetas(*reco_daughter_track_thetasRV);
    thin_slice_event.SetRecoDaughterTrackScores(*reco_daughter_track_scoresRV);
    thin_slice_event.SetHasPi0Shower(*has_pi0_showerRV);
    thin_slice_event.SetLeadingPCostheta(*leading_p_costhetaRV);
    thin_slice_event.SetLeadingPiPlusCostheta(*leading_piplus_costhetaRV);
    thin_slice_event.SetLeadingPi0Costheta(*leading_pi0_costhetaRV);
    thin_slice_event.SetIsBeamScraper(*is_beam_scraperRV);
    for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {
      thin_slice_event.AddRecoDaughterTrackdQdX((*daughter_dQdXsRV)[j]);
      thin_slice_event.AddRecoDaughterTrackResRange((*daughter_resRangesRV)[j]);
      thin_slice_event.AddRecoDaughterEField((*daughter_EFieldsRV)[j]);
    }

    for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
      std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                   std::to_string(j);
      thin_slice_event.MakeG4RWCoeff(name_downstream,
                                  (*g4rw_downstream_grid_piplus_coeffsRV)[j]);
    }
    thin_slice_event.MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                                (*g4rw_full_grid_proton_coeffsRV)[0]);

    bool found_start = false;
    for (size_t j = 1; j < true_beam_traj_ZRV->size(); ++j) {
       
      if ((*true_beam_traj_ZRV)[j] < fTrajZStart) continue;

      double delta = (*true_beam_traj_KERV)[0] - (*true_beam_traj_KERV)[j];
      thin_slice_event.SetDeltaEToTPC(delta);
      found_start = true;
      break;
    }
    if (!found_start) {
      thin_slice_event.SetDeltaEToTPC(-999.);
    }

    std::lock_guard<std::mutex> guard(fFillMutex);
    events.push_back(thin_slice_event);
  }
};

tp.Process(fill_func);*/


for (int & i : event_list) {
  if (!(i % 20000)) std::cout << i << "/" << split_val << std::endl;
  fMCTree->GetEntry(i);

  //ThinSliceEvent thin_slice_event(event, subrun, run);
  events.push_back(ThinSliceEvent(event, subrun, run));
  events.back().SetMCStatVarWeight(fRNG.Poisson(1));
  events.back().SetSampleID(sample_ID);
  events.back().SetSelectionID(selection_ID);
  events.back().SetCalUpSelectionID(cal_up_selection_ID);
  events.back().SetCalDownSelectionID(cal_down_selection_ID);

  events.back().SetFrontUpSelectionID(selection_ID_front_shift_up);
  events.back().SetFrontDownSelectionID(selection_ID_front_shift_down);
  events.back().SetBackUpSelectionID(selection_ID_back_shift_up);
  events.back().SetBackDownSelectionID(selection_ID_back_shift_down);

  events.back().SetTrueInteractingEnergy(true_beam_interactingEnergy);
  // if (fDoEnergyByLen) {
  //   events.back().SetRecoInteractingEnergy((*reco_beam_incidentEnergies)[0] - 2.1*reco_beam_alt_len);
  // }
  // else {
    events.back().SetRecoInteractingEnergy(reco_beam_interactingEnergy);
  // }
  events.back().SetTrueEndP(true_beam_endP);
  events.back().SetTrueEndZ(true_beam_endZ);
  events.back().SetTrueStartP(true_beam_startP);
  events.back().SetTrueMass(true_beam_mass);
  events.back().SetRecoEndZ(reco_beam_endZ);
  events.back().SetRecoStartY(reco_beam_startY);
  events.back().SetRecoStartX_SCE(reco_beam_startX_SCE);
  events.back().SetRecoStartY_SCE(reco_beam_startY_SCE);
  events.back().SetRecoStartZ_SCE(reco_beam_startZ_SCE);
  events.back().SetRecoEndZ_SCE(reco_beam_endZ_SCE);

  events.back().SetVertexMichelScore(vertex_michel_score);
  events.back().SetVertexNHits(vertex_nhits);

  events.back().SetTrueID(true_beam_ID);
  events.back().SetRecoToTrueID(reco_beam_true_byHits_ID);

  events.back().SetRecoIncidentEnergies(*reco_beam_incidentEnergies);
  events.back().SetTrueInitEnergy(
    (true_beam_incidentEnergies->size() > 0 ?
     true_beam_incidentEnergies->at(0) : -999.)); // TODO -- CHANGE THIS WITH NEW VAL 
  events.back().SetTrueTrajZ(*true_beam_traj_Z);
  events.back().SetTrueTrajKE(*true_beam_traj_KE);
  events.back().SetTrueNNeutrons(true_n_neutrons);
  events.back().SetTrueNProtons(true_n_protons);
  events.back().SetTrueNPiPlus(true_n_piplus);
  events.back().SetTrueNPiMinus(true_n_piminus);
  events.back().SetTrueNPi0(true_n_pi0);

  events.back().SetBeamInstP(beam_inst_P);
  events.back().SetPDG(true_beam_PDG);
  events.back().SetRecoDaughterTrackScores(*reco_daughter_track_scores);
  events.back().SetHasPi0Shower(has_pi0_shower);
  events.back().SetLeadingPCostheta(leading_p_costheta);
  events.back().SetLeadingPiPlusCostheta(leading_piplus_costheta);
  events.back().SetLeadingPi0Costheta(leading_pi0_costheta);
  events.back().SetLeadingPMomentum(leading_p_momentum);
  events.back().SetLeadingPiPlusMomentum(leading_piplus_momentum);
  events.back().SetLeadingPi0Momentum(leading_pi0_momentum);
  events.back().SetIsBeamScraper(is_beam_scraper);

  events.back().SetRecoShowerEnergy(*reco_daughter_shower_energy);

  for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {
    events.back().AddRecoDaughterEField((*daughter_EFields)[j]);
    auto nhits = (*reco_daughter_chi2_nhits)[j];
    events.back().AddOneChi2PerHit(
        (nhits > 0) ?
        (*reco_daughter_chi2)[j]/nhits :
        -999.
    );
    events.back().AddOneTrunc_dEdX(
      (*reco_daughter_truncated_dEdX)[j]
    );
  }

  for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
    std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                 std::to_string(j);
    events.back().MakeG4RWCoeff(name_downstream,
                                (*g4rw_downstream_grid_piplus_coeffs)[j]);

    std::string name_full = "g4rw_full_grid_piplus_coeffs_" + std::to_string(j);
    events.back().MakeG4RWCoeff(name_full, (*g4rw_full_grid_piplus_coeffs)[j]);
  }

  for (size_t j = 0; j < g4rw_full_fine_piplus_coeffs->size(); ++j) {
    std::string fine_full = "g4rw_full_fine_piplus_coeffs_" + std::to_string(j);
    events.back().MakeG4RWCoeff(fine_full, (*g4rw_full_fine_piplus_coeffs)[j]);
  }

  for (size_t j = 0; j < g4rw_full_grid_abscex_coeffs->size(); ++j) {
    std::string abscex = "g4rw_full_grid_abscex_coeffs_" + std::to_string(j);
    events.back().MakeG4RWCoeff(abscex, (*g4rw_full_grid_abscex_coeffs)[j]);
  }


  events.back().MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                              (*g4rw_full_grid_proton_coeffs)[0]);

  bool found_start = false;
  for (size_t j = 1; j < true_beam_traj_Z->size(); ++j) {
     
    if ((*true_beam_traj_Z)[j] < fTrajZStart) continue;

    double delta = (*true_beam_traj_KE)[0] - (*true_beam_traj_KE)[j];
    events.back().SetDeltaEToTPC(delta);
    found_start = true;
    break;
  }
  if (!found_start) {
    events.back().SetDeltaEToTPC(-999.);
  }

}

if (do_fake_data) {
  std::cout << "Filling fake data " << fake_event_list.size() << std::endl;
  //for (int i = fake_start; i < fake_end; ++i) {
  auto start_time = std::chrono::high_resolution_clock::now();
  for (int & i : fake_event_list) {
    //std::cout << "Fake event " << i << std::endl;
    fMCTree->GetEntry(i);

    fake_data_events.push_back(ThinSliceEvent(event, subrun, run));
    if (!(fake_data_events.size() % 20000))
      std::cout << fake_data_events.size() << std::endl;
    fake_data_events.back().SetMCStatVarWeight(fRNG.Poisson(1));
    fake_data_events.back().SetSampleID(sample_ID);
    fake_data_events.back().SetSelectionID(selection_ID);
    fake_data_events.back().SetCalUpSelectionID(cal_up_selection_ID);
    fake_data_events.back().SetCalDownSelectionID(cal_down_selection_ID);

    fake_data_events.back().SetFrontUpSelectionID(selection_ID_front_shift_up);
    fake_data_events.back().SetFrontDownSelectionID(selection_ID_front_shift_down);
    fake_data_events.back().SetBackUpSelectionID(selection_ID_back_shift_up);
    fake_data_events.back().SetBackDownSelectionID(selection_ID_back_shift_down);

    fake_data_events.back().SetTrueInteractingEnergy(true_beam_interactingEnergy);
    // if (fDoEnergyByLen) {
    //   fake_data_events.back().SetRecoInteractingEnergy((*reco_beam_incidentEnergies)[0] - 2.1*reco_beam_alt_len);
    // }
    // else {
      fake_data_events.back().SetRecoInteractingEnergy(reco_beam_interactingEnergy);
    // }
    fake_data_events.back().SetTrueEndP(true_beam_endP);
    fake_data_events.back().SetTrueEndZ(true_beam_endZ);
    fake_data_events.back().SetTrueStartP(true_beam_startP);
    fake_data_events.back().SetTrueMass(true_beam_mass);
    fake_data_events.back().SetRecoEndZ(reco_beam_endZ);
    fake_data_events.back().SetRecoStartX_SCE(reco_beam_startX_SCE);
    fake_data_events.back().SetRecoStartY_SCE(reco_beam_startY_SCE);
    fake_data_events.back().SetRecoStartZ_SCE(reco_beam_startZ_SCE);
    fake_data_events.back().SetRecoEndZ_SCE(reco_beam_endZ_SCE);
    fake_data_events.back().SetVertexMichelScore(vertex_michel_score);
    fake_data_events.back().SetVertexNHits(vertex_nhits);

    fake_data_events.back().SetRecoIncidentEnergies(*reco_beam_incidentEnergies);
    fake_data_events.back().SetTrueInitEnergy(
        (true_beam_incidentEnergies->size() > 0 ?
         true_beam_incidentEnergies->at(0) : -999.));
    fake_data_events.back().SetTrueTrajZ(*true_beam_traj_Z);
    fake_data_events.back().SetTrueTrajKE(*true_beam_traj_KE);

    fake_data_events.back().SetBeamInstP(beam_inst_P);
    fake_data_events.back().SetPDG(true_beam_PDG);
    if (fFakeDataRoutine == "EffVar") {
      fake_data_events.back().SetRecoDaughterTrackThetas(*reco_daughter_track_thetas); }
    fake_data_events.back().SetTrueDaughterPDGs(*true_beam_daughter_PDG);
    fake_data_events.back().SetTrueDaughterStartPs(*true_beam_daughter_startP);
    fake_data_events.back().SetRecoDaughterTrackScores(*reco_daughter_track_scores);
    fake_data_events.back().SetHasPi0Shower(has_pi0_shower);
    fake_data_events.back().SetLeadingPCostheta(leading_p_costheta);
    fake_data_events.back().SetLeadingPiPlusCostheta(leading_piplus_costheta);
    fake_data_events.back().SetLeadingPi0Costheta(leading_pi0_costheta);
    fake_data_events.back().SetLeadingPMomentum(leading_p_momentum);
    fake_data_events.back().SetLeadingPiPlusMomentum(leading_piplus_momentum);
    fake_data_events.back().SetLeadingPi0Momentum(leading_pi0_momentum);
    fake_data_events.back().SetIsBeamScraper(is_beam_scraper);
    fake_data_events.back().SetRecoShowerEnergy(*reco_daughter_shower_energy);

    fake_data_events.back().SetTrueID(true_beam_ID);
    fake_data_events.back().SetTrueNNeutrons(true_n_neutrons);
    fake_data_events.back().SetTrueNProtons(true_n_protons);
    fake_data_events.back().SetTrueNPiPlus(true_n_piplus);
    fake_data_events.back().SetTrueNPiMinus(true_n_piminus);
    fake_data_events.back().SetTrueNPi0(true_n_pi0);
    fake_data_events.back().SetRecoToTrueID(reco_beam_true_byHits_ID);

    for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {

      fake_data_events.back().AddRecoDaughterEField((*daughter_EFields)[j]);
      auto nhits = (*reco_daughter_chi2_nhits)[j];
      fake_data_events.back().AddOneChi2PerHit(
          (nhits > 0) ?
          (*reco_daughter_chi2)[j]/nhits :
          -999.
      );
      fake_data_events.back().AddOneTrunc_dEdX(
        (*reco_daughter_truncated_dEdX)[j]
      );
    }

    for (size_t j = 0; j < g4rw_full_fine_piplus_coeffs->size(); ++j) {
      std::string fine_full = "g4rw_full_fine_piplus_coeffs_" + std::to_string(j);
      fake_data_events.back().MakeG4RWCoeff(fine_full, (*g4rw_full_fine_piplus_coeffs)[j]);
    }
    for (size_t j = 0; j < g4rw_full_grid_abscex_coeffs->size(); ++j) {
      std::string abscex = "g4rw_full_grid_abscex_coeffs_" + std::to_string(j);
      fake_data_events.back().MakeG4RWCoeff(abscex, (*g4rw_full_grid_abscex_coeffs)[j]);

    }

    for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
      std::string name_full = "g4rw_full_grid_piplus_coeffs_" + std::to_string(j);
      fake_data_events.back().MakeG4RWCoeff(name_full, (*g4rw_full_grid_piplus_coeffs)[j]);


      std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                   std::to_string(j);
      fake_data_events.back().MakeG4RWCoeff(name_downstream,
                                  (*g4rw_downstream_grid_piplus_coeffs)[j]);
    }

    fake_data_events.back().MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                                          (*g4rw_full_grid_proton_coeffs)[0]);
    bool found_start = false;
    for (size_t j = 1; j < true_beam_traj_Z->size(); ++j) {
      if ((*true_beam_traj_Z)[j-1] < fTrajZStart &&
          fTrajZStart < (*true_beam_traj_Z)[j+1]) {
        double delta = (*true_beam_traj_KE)[0] - (*true_beam_traj_KE)[j];
        fake_data_events.back().SetDeltaEToTPC(delta);
        found_start = true;
        break;
      }
    }
    if (!found_start) fake_data_events.back().SetDeltaEToTPC(-999.);
  }
  auto new_time = std::chrono::high_resolution_clock::now();
  auto delta =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          new_time - start_time).count();
  std::cout << "Filling " << fake_data_events.size() <<
               " fake data took " << delta << std::endl;
}

std::cout << "Filled MC Events" << std::endl;


}


std::pair<int, int> protoana::ThinSliceDriver::GetColorAndStyle(
    size_t i,
    const std::vector<std::pair<int, int>> & plot_style) {
  return {plot_style.at(i % plot_style.size()).first,
          (i < plot_style.size() ? 1001: 3244)};
}

void protoana::ThinSliceDriver::ResetSamples(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples) {
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].Reset();
      }
    }
  }
}
