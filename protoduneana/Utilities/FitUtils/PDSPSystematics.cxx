#include "PDSPSystematics.h"
#include <chrono>

double protoana::PDSPSystematics::fEDivF = 0.;
double protoana::PDSPSystematics::fEDivCut = 0.;
const protoana::ThinSliceSystematic * protoana::PDSPSystematics::fEDivPar = 0x0;
int protoana::PDSPSystematics::fPastFVSelectionID = 0;
std::vector<const protoana::ThinSliceSystematic *> protoana::PDSPSystematics::fActiveG4RWPars = {};
std::vector<const protoana::ThinSliceSystematic *> protoana::PDSPSystematics::fActiveTiedG4RWPars = {};

std::vector<const protoana::ThinSliceSystematic *> protoana::PDSPSystematics::fActiveSignalG4RWPars = {};
std::vector<const protoana::ThinSliceSystematic *> protoana::PDSPSystematics::fActiveTiedSignalG4RWPars = {};

std::vector<const protoana::ThinSliceSystematic *> protoana::PDSPSystematics::fActiveCosThetaEffPars = {};
std::vector<int> protoana::PDSPSystematics::fCosThetaEffIDs = {};
std::vector<int> protoana::PDSPSystematics::fCosThetaEffSels = {};
std::vector<std::pair<double, double>> protoana::PDSPSystematics::fCosThetaEffRanges = {};
std::vector<double> protoana::PDSPSystematics::fCosThetaEffFracs = {};
std::vector<int> protoana::PDSPSystematics::fCosThetaEffPDGs = {};


int protoana::PDSPSystematics::fUpstreamID = 0;
int protoana::PDSPSystematics::fNoTrackID = 0;
double protoana::PDSPSystematics::fEndZNoTrackCut = 0.;
const protoana::ThinSliceSystematic * protoana::PDSPSystematics::fEndZPar = 0x0;
std::map<int, std::vector<double>> protoana::PDSPSystematics::fEndZFractions = {};

std::map<int, std::vector<double>> protoana::PDSPSystematics::fBeamMatchFractions = {};
const protoana::ThinSliceSystematic * protoana::PDSPSystematics::fBeamMatchPar = 0x0;
double protoana::PDSPSystematics::fBeamMatchFraction = 0.;
bool protoana::PDSPSystematics::fBeamMatchUseSingleFrac = false;

double protoana::PDSPSystematics::fBeamMatchLowLimit = 0.;
double protoana::PDSPSystematics::fBeamMatchLowFraction = 0.;
std::map<int, std::vector<double>> protoana::PDSPSystematics::fBeamMatchLowFractions = {};
bool protoana::PDSPSystematics::fBeamMatchLowUseSingleFrac = false;
const protoana::ThinSliceSystematic * protoana::PDSPSystematics::fBeamMatchLowPar = 0x0;

double protoana::PDSPSystematics::fBeamMatchHighLimit = 0.;
double protoana::PDSPSystematics::fBeamMatchHighFraction = 0.;
std::map<int, std::vector<double>> protoana::PDSPSystematics::fBeamMatchHighFractions = {};
bool protoana::PDSPSystematics::fBeamMatchHighUseSingleFrac = false;
const protoana::ThinSliceSystematic * protoana::PDSPSystematics::fBeamMatchHighPar = 0x0;

const protoana::ThinSliceSystematic * protoana::PDSPSystematics::fBeamScraperPar = 0x0;


protoana::PDSPSystematics::PDSPSystematics(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<std::string, ThinSliceSystematic> & pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    TFile & output_file, int upstream_ID, int no_track_ID, int decay_ID,
    int past_FV_ID, int beam_cut_ID, int past_FV_sel_ID)
    : /*fUpstreamID(upstream_ID), fNoTrackID(no_track_ID),*/ fDecayID(decay_ID),
      fPastFVID(past_FV_ID), fBeamCutID(beam_cut_ID)/*,
      fPastFVSelectionID(past_FV_sel_ID)*/ {

  fPastFVSelectionID = past_FV_sel_ID;
  fUpstreamID = upstream_ID;
  fNoTrackID = no_track_ID;


  SetupSyst_G4RWCoeff(pars);
  SetupSyst_TiedG4RWCoeff(pars);
  SetupSyst_BeamShift(pars, output_file);//Trim
  SetupSyst_EDivWeight(pars);//Done
  SetupSyst_EndZNoTrackWeight(pars);//Done
  SetupSyst_CosThetaEff(pars);
  SetupSyst_BeamMatch(pars);//Done
  SetupSyst_BeamMatchHigh(pars);//Done
  SetupSyst_BeamMatchLow(pars);//Done
  SetupSyst_BoxBeam(pars);//Trim
  SetupSyst_ELoss(pars);//Trim
  SetupSyst_ELossMuon(pars);//Trim
  SetupSyst_QuadBeamShift(pars);//Trim

  SetupSyst_BeamShiftBins(pars);//Trim
  SetupSyst_BeamScraper(pars);

  SetupSignal_G4RWCoeff(g4rw_pars);
  SetupSignal_TiedG4RWCoeff(g4rw_pars);
}

double protoana::PDSPSystematics::GetEventWeight(
    const ThinSliceEvent & event,
    int signal_index,
    //int selection_bin,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  double weight = 1.;
  //auto begin_time = std::chrono::high_resolution_clock::now();

  /*
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->second.GetIsG4RWCoeff()) {
      weight *= GetSystWeight_G4RWCoeff(event, it->second);
    }
    else if (it->second.GetIsTiedG4RWCoeff()) {
      weight *= GetSystWeight_TiedG4RWCoeff(event, it->second);
    }
    else if (it->first == "ediv_weight") {
      weight *= GetSystWeight_EDiv(event, it->second);
    }
    else if (it->first == "end_z_no_track_weight") {
      weight *= GetSystWeight_EndZNoTrack(event, signal_index, it->second);
    }
    else if (it->first == "beam_match_weight") {
      weight *= GetSystWeight_BeamMatch(event, signal_index, it->second);
    }
    else if (it->first == "beam_match_low_weight") {
      weight *= GetSystWeight_BeamMatchLow(event, signal_index, it->second);
    }
    else if (it->first == "beam_match_high_weight") {
      weight *= GetSystWeight_BeamMatchHigh(event, signal_index, it->second);
    }
    else if (it->first == "eloss_weight") {
      weight *= GetSystWeight_ELoss(event, signal_index, it->second);
    }
    else if (it->first == "eloss_weight_muon") {
      weight *= GetSystWeight_ELossMuon(event, it->second);
    }
    else if (it->first == "quad_beam_shift") {
      weight *= GetSystWeight_QuadBeamShift(event, it->second);
    }
  }*/

  //auto end_time = std::chrono::high_resolution_clock::now();
  //auto delta =
  //    std::chrono::duration_cast<std::chrono::nanoseconds>(
  //        end_time - begin_time).count();
  //std::cout << "\tWeighting took " << delta << " ns" << std::endl;



  for (auto syst : fActiveSystematics) {
    double contribution = syst.second(event, signal_index);
    //std::cout << syst.first << " " << contribution <<
    //             std::endl;
    weight *= contribution;
  }

  return weight;
}

void protoana::PDSPSystematics::SetupSyst_CosThetaEff(
  const std::map<std::string, ThinSliceSystematic> & pars) {
  for (const auto & par : pars) {
    if (par.second.GetIsCosThetaEff()) {
      std::cout << "Setting up costheta eff syst: " << par.first << std::endl;
    }
    else continue;

    fCosThetaEffIDs.push_back(par.second.GetOption<int>("ID"));
    fCosThetaEffSels.push_back(par.second.GetOption<int>("Sel"));
    fCosThetaEffRanges.push_back(
        par.second.GetOption<std::pair<double, double>>("Range"));
    fCosThetaEffFracs.push_back(
        par.second.GetOption<double>("Frac"));
    fCosThetaEffPDGs.push_back(
        par.second.GetOption<int>("PDG"));
    fActiveCosThetaEffPars.push_back(&(par.second));
  }

  fActiveSystematics.push_back({"CosThetaEff", GetSystWeight_CosThetaEff});
}

double protoana::PDSPSystematics::GetSystWeight_CosThetaEff(
    const ThinSliceEvent & event, int signal_index) {
  double weight = 1.;
  size_t ipar = 0;
  for (const auto * par : fActiveCosThetaEffPars) {

    if (event.GetSampleID() != fCosThetaEffIDs[ipar]) continue;
    //std::cout << par->GetName() << " " << par->GetValue() << " " <<
    //           event.GetSampleID() << " " << fCosThetaEffIDs[ipar] << " ";

    double cos = -1.;
    if (fCosThetaEffPDGs[ipar] == 211) {
      cos = event.GetLeadingPiPlusCostheta();
    }
    else if (fCosThetaEffPDGs[ipar] == 2212) {
      cos = event.GetLeadingPCostheta();
    }
    else if (fCosThetaEffPDGs[ipar] == 111) {
      cos = event.GetLeadingPi0Costheta();
    }

    auto range = fCosThetaEffRanges[ipar];
    //std::cout << cos << " " << range.first << " " << range.second;
    if (cos < range.first || cos > range.second) continue;

    double c = par->GetValue();
    double f = fCosThetaEffFracs[ipar];
    double subweight = (
      event.GetSelectionID() == fCosThetaEffSels[ipar] ?
      c : (1. - c*f)/(1. - f)
    );
    //std::cout << " " << c << " " << f << " " << subweight << std::endl;

    weight *= CheckAndReturn(subweight, "CosThetaEff", (*par), event);
  }

  return weight;
}

void protoana::PDSPSystematics::SetupSyst_G4RWCoeff(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->second.GetIsG4RWCoeff()/*(it->first.find("g4rw") != std::string::npos ||
         it->first.find("G4RW") != std::string::npos) &&
        (it->first.find("oeff") != std::string::npos) &&
        (it->first.find("tied") == std::string::npos)*/) {
      std::cout << "Setting up g4rw coeffsyst: " << it->first << std::endl;
    }
    else {
      continue;
    }
    fG4RWCoeffBranches[it->first] = it->second.GetOption<std::string>("Branch");
    fActiveG4RWPars.push_back(&(it->second));
    std::cout << "Added " << it->first << " " << fG4RWCoeffBranches[it->first] << std::endl;
  }
  fActiveSystematics.push_back({"G4RWCoeff", GetSystWeight_G4RWCoeffNoPar});
}

void protoana::PDSPSystematics::SetupSignal_G4RWCoeff(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->second.GetIsG4RWCoeff()) {
      std::cout << "Setting up g4rw coeffsyst: " << it->first << std::endl;
    }
    else {
      continue;
    }
    fActiveSignalG4RWPars.push_back(&(it->second));
  }
}

void protoana::PDSPSystematics::SetupSyst_TiedG4RWCoeff(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->second.GetIsTiedG4RWCoeff()/*(it->first.find("g4rw") != std::string::npos ||
         it->first.find("G4RW") != std::string::npos) &&
        (it->first.find("oeff") != std::string::npos) &&
        (it->first.find("tied") == std::string::npos)*/) {
      std::cout << "Setting up g4rw tied coeffsyst: " << it->first << std::endl;
    }
    else {
      continue;
    }
    //fG4RWCoeffBranches[it->first] = it->second.GetOption<std::string>("Branches");
    fActiveTiedG4RWPars.push_back(&(it->second));
    std::cout << "Added " << it->first << " " << fG4RWCoeffBranches[it->first] << std::endl;
  }
  fActiveSystematics.push_back({"TiedG4RWCoeff", GetSystWeight_TiedG4RWCoeffNoPar});
}

void protoana::PDSPSystematics::SetupSignal_TiedG4RWCoeff(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->second.GetIsTiedG4RWCoeff()) {
      std::cout << "Setting up g4rw tied coeffsyst: " << it->first << std::endl;
    }
    else {
      continue;
    }
    fActiveTiedSignalG4RWPars.push_back(&(it->second));
  }
}

void protoana::PDSPSystematics::SetupSyst_BeamShiftBins(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->first.find("beam_shift_bin") != std::string::npos) {
      std::cout << "Setting up beam_shift_bin syst: " << it->first << std::endl;
    }
    else {
      continue;
    }

    fBeamShiftBinMeans[it->second.GetName()] = it->second.GetOption<double>("Mean");
    fBeamShiftBinSigmas[it->second.GetName()] = it->second.GetOption<double>("Sigma");
    fBeamShiftBinRangeLows[it->second.GetName()] = it->second.GetOption<double>("RangeLow");
    fBeamShiftBinRangeHighs[it->second.GetName()] = it->second.GetOption<double>("RangeHigh");
    std::cout << "Added " << it->second.GetName() << " " << std::endl;
  }
}

//optimize this -- TODO 
double protoana::PDSPSystematics::GetSystWeight_BeamShiftBin(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
  if (event.GetPDG() != 211 || event.GetTrueIncidentEnergies().size() == 0)
    return 1.;

  double beam_inst_KE = sqrt(std::pow(event.GetBeamInstP()*1.e3, 2) +
                             139.57*139.57) - 139.57;
  //std::cout << beam_inst_KE << " " << fBeamShiftBinRangeLows[par.GetName()] <<
  //             " " << fBeamShiftBinRangeHighs[par.GetName()] << " " <<
  //             par.GetName() << std::endl;

  if (beam_inst_KE < fBeamShiftBinRangeLows[par.GetName()] ||
      beam_inst_KE >= fBeamShiftBinRangeHighs[par.GetName()]) {
    //std::cout << "returning 1" << std::endl;
    return 1.;
  }

  double x = (beam_inst_KE - event.GetTrueIncidentEnergies()[0]);

  double mean = fBeamShiftBinMeans[par.GetName()];
  double stddev = fBeamShiftBinSigmas[par.GetName()];
  double weight = (exp(.5*((std::pow(x - mean, 2)) -
                   (std::pow(x - mean + par.GetValue(), 2)))/(stddev*stddev)));
  //std::cout << "Returning " << par.GetValue() << " " << weight << std::endl;
  return weight;

}

double protoana::PDSPSystematics::GetSystWeight_G4RWCoeff(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars) {
  double weight = 1.;

  //for (auto it = fG4RWCoeffBranches.begin(); it != fG4RWCoeffBranches.end();
  //     ++it) {
    weight *= event.GetG4RWCoeffWeight(
        par.GetG4RWCoeffBranch(), par.GetValue());
    //std::cout << event.GetEventID() << " " << event.GetSubrunID() << " " <<
    //             event.GetRunID() << " " << it->first << " " <<
    //             pars.at(it->first).GetValue() << " " << weight << std::endl;
  //}
  return weight;
  //return CheckAndReturn(weight, "G4RWCoeff", par);
}

double protoana::PDSPSystematics::GetSystWeight_G4RWCoeffNoPar(
    const ThinSliceEvent & event, int signal_index) {
  double weight = 1.;

  for (const auto * par : fActiveG4RWPars) {
    double subweight = event.GetG4RWCoeffWeight(
        par->GetG4RWCoeffBranch(), par->GetValue(), par->GetG4RWExtend());
    weight *= CheckAndReturn(subweight, "G4RWCoeff", (*par), event);
    //weight *= subweight;
    //std::cout << "\t" << par->GetG4RWCoeffBranch() << std::endl;
  }
  return weight;
  //return CheckAndReturn(weight, "G4RWCoeff", par, event);
}

double protoana::PDSPSystematics::GetSignalWeight_G4RWCoeffNoPar(
    const ThinSliceEvent & event, int signal_index) {
  double weight = 1.;

  for (const auto * par : fActiveSignalG4RWPars) {
    double subweight = event.GetG4RWCoeffWeight(
        par->GetG4RWCoeffBranch(), par->GetValue(), par->GetG4RWExtend());
    weight *= CheckAndReturn(subweight, "SignalG4RWCoeff", (*par), event);
    //weight *= subweight;
    //std::cout << "\t" << par->GetG4RWCoeffBranch() << std::endl;
  }
  return weight;
  //return CheckAndReturn(weight, "G4RWCoeff", par, event);
}

double protoana::PDSPSystematics::GetSystWeight_TiedG4RWCoeff(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
  double weight = 1.;

  const auto & branches = par.GetTiedG4RWBranches();
  for (const auto & br : branches) {
    //std::cout << br << " " << par.GetValue() << " " <<
    //             event.GetG4RWCoeffWeight(br, par.GetValue()) << std::endl;
    weight *= event.GetG4RWCoeffWeight(
        br, par.GetValue());
  }
  return weight;
  //return CheckAndReturn(weight, "G4RWCoeffTied", par);
}

double protoana::PDSPSystematics::GetSystWeight_TiedG4RWCoeffNoPar(
    const ThinSliceEvent & event, int signal_index) {
  double weight = 1.;

  for (const auto * par : fActiveTiedG4RWPars) {
    const auto & branches = par->GetTiedG4RWBranches();
    //std::cout << "Tied g4rw: " << par->GetName() << std::endl;
    for (const auto & br : branches) {
      //std::cout << "\t" << br << " ";
      double subweight = event.GetG4RWCoeffWeight(br, par->GetValue(), false);
      //std::cout << subweight << std::endl;
      weight *= CheckAndReturn(subweight, "G4RWTiedCoeff", (*par), event);
    }
  }

  //std::cout << "\tTotal: " << weight << std::endl;
  return weight;
}

double protoana::PDSPSystematics::GetSignalWeight_TiedG4RWCoeffNoPar(
    const ThinSliceEvent & event, int signal_index) {
  double weight = 1.;

  for (const auto * par : fActiveTiedSignalG4RWPars) {
    const auto & branches = par->GetTiedG4RWBranches();
    //std::cout << "Tied g4rw: " << par->GetName() << std::endl;
    for (const auto & br : branches) {
      //std::cout << "\t" << br << " ";
      double subweight = event.GetG4RWCoeffWeight(br, par->GetValue(), false);
      //std::cout << subweight << std::endl;
      weight *= CheckAndReturn(subweight, "SignalG4RWTiedCoeff", (*par), event);
    }
  }

  //std::cout << "\tTotal: " << weight << std::endl;
  return weight;
}


void protoana::PDSPSystematics::SetupSyst_BeamShift(
    const std::map<std::string, ThinSliceSystematic> & pars,
    TFile & output_file) {
  if (pars.find("beam_shift") == pars.end()) {
    return;
  }

  TFile shift_file(
      pars.at("beam_shift").GetOption<std::string>("ShiftFile").c_str());
  //fSystBeamShiftMap = (TGraph2D*)shift_file.Get("g2d");
  //fSystBeamShiftMap->SetDirectory(0);
  fSystBeamShiftLimits
      = pars.at("beam_shift").GetOption<std::pair<double, double>>("Limits");

  //fSystBeamShiftRatioLimitUp = pars.at("beam_shift").GetOption<double>("RatioLimitUp");
  //fSystBeamShiftRatioLimitDown = pars.at("beam_shift").GetOption<double>("RatioLimitDown");
  fSystBeamShiftMeans = (TGraph*)shift_file.Get("gMeans");
  //fSystBeamShiftMeans->SetDirectory(0);
  fSystBeamShiftWidths = (TGraph*)shift_file.Get("gWidths");
  //fSystBeamShiftWidths->SetDirectory(0);
  shift_file.Close();

  fSystBeamShiftWeightCap = pars.at("beam_shift").GetOption<double>("WeightCap");
  /*
  fSystBeamShiftTreeSave = pars.at("beam_shift").GetOption<bool>("SaveInfo");
  if (fSystBeamShiftTreeSave) {
    output_file.mkdir("SystBeamShift");
    output_file.cd("SystBeamShift");
    fSystBeamShiftTree = new TTree("tree", "");
    fSystBeamShiftTree->Branch("Weight", &fSystBeamShiftWeight);
    fSystBeamShiftTree->Branch("Val", &fSystBeamShiftVal);
    fSystBeamShiftTree->Branch("R", &fSystBeamShiftR);
  }
  */
  fSetupSystBeamShift = true;

}

double protoana::PDSPSystematics::GetSystWeight_BeamShift(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars) {
  //if (pars.find("beam_shift") == pars.end()) return 1.;
  if (event.GetPDG() != 211) return 1.;
  double x_val = par/*s.at("beam_shift")*/.GetValue();
  double y_val = (event.GetBeamInstP() - event.GetTrueStartP())/
                  event.GetTrueStartP();
  if (y_val < fSystBeamShiftLimits.first ||
      y_val > fSystBeamShiftLimits.second) {
    return 1.;
  }

  double nominal_mean = fSystBeamShiftMeans->Eval(0.);
  double nominal_width = fSystBeamShiftWidths->Eval(0.);
  double varied_mean = fSystBeamShiftMeans->Eval(x_val);
  double varied_width = fSystBeamShiftWidths->Eval(x_val);

  double weight = (nominal_width/varied_width)*
                  exp(.5*std::pow(((y_val - nominal_mean)/nominal_width), 2)
                      - .5*std::pow(((y_val - varied_mean)/varied_width), 2));

  if (weight > fSystBeamShiftWeightCap && fSystBeamShiftWeightCap > 0.) {
    weight = fSystBeamShiftWeightCap;
  }

  /*
  if (fSystBeamShiftTreeSave) {
    fSystBeamShiftWeight = weight;
    fSystBeamShiftVal = x_val;
    fSystBeamShiftR = y_val;
    fSystBeamShiftTree->Fill();
  }*/
  //return weight;
  return CheckAndReturn(weight, "BeamShift", par, event);
}

void protoana::PDSPSystematics::SetupSyst_EDivWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("ediv_weight") == pars.end()) {
    return;
  }
  fEDivF = pars.at("ediv_weight").GetOption<double>("F");
  fEDivCut = pars.at("ediv_weight").GetOption<double>("Cut");

  fActiveSystematics.push_back({"ediv_weight", GetSystWeight_EDivNoPar});
  fEDivPar = &(pars.at("ediv_weight"));
}

double protoana::PDSPSystematics::GetSystWeight_EDiv(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int ediv_selection_ID) {
  
  //if (pars.find("ediv_weight") == pars.end()) return 1.;

  const int selection_ID = event.GetSelectionID();
  if (selection_ID != fPastFVSelectionID/*ediv_selection_ID4*/) return 1.;

  const double endZ = event.GetRecoEndZ();
 
  double weight = 1.;
  double var = par/*s.at("ediv_weight")*/.GetValue();
  if (endZ < fEDivCut) {
    weight = var;
  }
  else {
    weight = (1. - var*fEDivF)/(1. - fEDivF);
  }

  if (weight < 0.) {
    std::cout << endZ << " " << fEDivCut << " " << var << " " << fEDivF <<
                 std::endl;
  }

  //return weight;
  return CheckAndReturn(weight, "EDivWeight", par, event);
}

double protoana::PDSPSystematics::GetSystWeight_EDivNoPar(
    const ThinSliceEvent & event, int signal_index) {
  
  const int selection_ID = event.GetSelectionID();
  if (selection_ID != fPastFVSelectionID) return 1.;

  const double endZ = event.GetRecoEndZ();
 
  double weight = 1.;
  double var = fEDivPar->GetValue();
  if (endZ < fEDivCut) {
    weight = var;
  }
  else {
    weight = (1. - var*fEDivF)/(1. - fEDivF);
  }

  if (weight < 0.) {
    std::cout << endZ << " " << fEDivCut << " " << var << " " << fEDivF <<
                 std::endl;
  }

  //return weight;
  //std::cout << "EDivNoPar " << weight << std::endl;
  //return weight;
  return CheckAndReturn(weight, "EDivWeight", (*fEDivPar), event);
}

void protoana::PDSPSystematics::SetupSyst_EndZNoTrackWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("end_z_no_track_weight") != pars.end()) {
    fEndZNoTrackCut = pars.at("end_z_no_track_weight").GetOption<double>("Cut");

    std::vector<std::pair<int, std::vector<double>>> temp
        = pars.at("end_z_no_track_weight").GetOption
            <std::vector<std::pair<int, std::vector<double>>>>("Fractions");
    fEndZFractions = std::map<int, std::vector<double>>(temp.begin(), temp.end());
    std::cout << "EndZ Fracs" << std::endl;
    for (auto it = fEndZFractions.begin(); it != fEndZFractions.end(); ++it) {
      std::cout << it->first << " ";
      for (const auto & f : it->second) {
        std::cout << f << " "; 
      }
      std::cout << std::endl;
    }

    fEndZPar = &(pars.at("end_z_no_track_weight"));
    fActiveSystematics.push_back(
        {"end_z_no_track_weight", GetSystWeight_EndZNoTrackNoPar});
  }
}

double protoana::PDSPSystematics::GetSystWeight_EndZNoTrack(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  if (fEndZFractions[event.GetSampleID()].size() == 0) return 1.;

  if (event.GetTrueEndZ() > fEndZNoTrackCut) return 1.;

  double variation = par.GetValue();

  double fraction = GetFractionBySample(fEndZFractions, event.GetSampleID(),
                                        signal_index);
  //double fraction = (signal_index > -1 ?
  //                   fEndZFractions[event.GetSampleID()][signal_index] :
  //                   fEndZFractions[event.GetSampleID()].back());

  if (fraction < 0.) {
    return 1.;
  }

  double weight = 
      (event.GetSelectionID() == fNoTrackID ? variation :
       (1. - variation*fraction)/(1. - fraction));
  if (weight < 0.) {
    std::cout << event.GetSelectionID() << " " << event.GetSampleID() << " " <<
                 signal_index << " " << variation << " " <<
                 fraction << std::endl;
  }
  return CheckAndReturn(weight, "EndZNoTrack", par, event);

}


double protoana::PDSPSystematics::GetSystWeight_EndZNoTrackNoPar(
    const ThinSliceEvent & event,
    int signal_index) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  if (fEndZFractions[event.GetSampleID()].size() == 0) return 1.;

  if (event.GetTrueEndZ() > fEndZNoTrackCut) return 1.;

  double variation = fEndZPar->GetValue();

  double fraction = GetFractionBySample(fEndZFractions, event.GetSampleID(),
                                        signal_index);
  //double fraction = (signal_index > -1 ?
  //                   fEndZFractions[event.GetSampleID()][signal_index] :
  //                   fEndZFractions[event.GetSampleID()].back());

  if (fraction < 0.) {
    return 1.;
  }

  double weight = 
      (event.GetSelectionID() == fNoTrackID ? variation :
       (1. - variation*fraction)/(1. - fraction));
  if (weight < 0.) {
    std::cout << event.GetSelectionID() << " " << event.GetSampleID() << " " <<
                 signal_index << " " << variation << " " <<
                 fraction << std::endl;
  }
  //return weight;
  return CheckAndReturn(weight, "EndZNoTrack", (*fEndZPar), event);

}


void protoana::PDSPSystematics::SetupSyst_BeamMatch(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_weight") == pars.end()) return; 

  //fBeamMatchLimits = pars.at("beam_match_weight")
  //    .GetOption<std::vector<double>>("Limits");
  //fBeamMatchFractions = pars.at("beam_match_weight")
  //    .GetOption<std::vector<double>>("Fractions");
  auto temp = pars.at("beam_match_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fBeamMatchFraction = pars.at("beam_match_weight")
      .GetOption<double>("Fraction");
  fBeamMatchFractions
      = std::map<int, std::vector<double>>(temp.begin(), temp.end()); 
  fBeamMatchUseSingleFrac = pars.at("beam_match_weight")
      .GetOption<bool>("UseSingleFrac");
  if (fBeamMatchUseSingleFrac) std::cout << "UseSingleFrac" << std::endl;
  fBeamMatchPar = &(pars.at("beam_match_weight"));
  fActiveSystematics.push_back(
      {"beam_match_weight", GetSystWeight_BeamMatchNoPar});
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatch(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  /*
  bool matched_to_cosmic = (event.GetRecoOrigin() == 2);
  double variation = par.GetValue();

  double fraction = GetFractionBySample(fBeamMatchFractions,
                                        event.GetSampleID(),
                                        signal_index);*/

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = par.GetValue();

  double fraction = GetFractionBySample(fBeamMatchFractions,
                                        event.GetSampleID(),
                                        signal_index);
  //std::cout << event.GetSampleID() << " " << signal_index << " " << fraction <<
  //             std::endl;
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));
  return CheckAndReturn(weight, "BeamMatch", par, event);
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatchNoPar(
    const ThinSliceEvent & event,
    int signal_index) {

  //Upstream Interactions
  //if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  bool matched = (event.GetRecoOrigin() == 2);
  //bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = fBeamMatchPar->GetValue();

  //double fraction = GetFractionBySample(fBeamMatchFractions,
  //                                      event.GetSampleID(),
  //                                      signal_index);
  double fraction = (fBeamMatchUseSingleFrac ?
                     fBeamMatchFraction :
                     GetFractionBySample(
                         fBeamMatchFractions, event.GetSampleID(),
                         signal_index));

  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));
  //return weight;
  return CheckAndReturn(weight, "BeamMatch", (*fBeamMatchPar), event);
}

void protoana::PDSPSystematics::SetupSyst_BeamMatchLow(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_low_weight") == pars.end()) return; 

  fBeamMatchLowLimit = pars.at("beam_match_low_weight")
      .GetOption<double>("Limit");
  fBeamMatchLowFraction = pars.at("beam_match_low_weight")
      .GetOption<double>("Fraction");
  auto temp = pars.at("beam_match_low_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fBeamMatchLowFractions
      = std::map<int, std::vector<double>>(temp.begin(), temp.end()); 
  fBeamMatchLowUseSingleFrac = pars.at("beam_match_low_weight")
      .GetOption<bool>("UseSingleFrac");
  fBeamMatchLowPar = &(pars.at("beam_match_low_weight"));
  fActiveSystematics.push_back(
      {"beam_match_low_weight", GetSystWeight_BeamMatchLowNoPar});
}

void protoana::PDSPSystematics::SetupSyst_BeamMatchHigh(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_high_weight") == pars.end()) return; 

  fBeamMatchHighLimit = pars.at("beam_match_high_weight")
      .GetOption<double>("Limit");
  fBeamMatchHighFraction = pars.at("beam_match_high_weight")
      .GetOption<double>("Fraction");
  auto temp = pars.at("beam_match_high_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fBeamMatchHighFractions
      = std::map<int, std::vector<double>>(temp.begin(), temp.end()); 
  fBeamMatchHighUseSingleFrac = pars.at("beam_match_high_weight")
      .GetOption<bool>("UseSingleFrac");
  fBeamMatchHighPar = &(pars.at("beam_match_high_weight"));
  fActiveSystematics.push_back(
      {"beam_match_high_weight", GetSystWeight_BeamMatchHighNoPar});
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatchLow(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  //upper limit
  if (event.GetTrueEndZ() > fBeamMatchLowLimit) return 1.;

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = par.GetValue();

  double fraction = (fBeamMatchLowUseSingleFrac ?
                     fBeamMatchLowFraction :
                     GetFractionBySample(
                         fBeamMatchLowFractions, event.GetSampleID(),
                         signal_index));
  //std::cout << event.GetSampleID() << " " << signal_index << " " << fraction <<
  //             std::endl;
  if (fraction < 0.) {
    return 1.; 
  }
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));

  /*if (fraction < 0.) {
    std::cout << "Low Fraction " << fraction << " weight " << weight <<
                 " matched? " << matched << " var " << variation << std::endl;
  }*/
  return CheckAndReturn(weight, "BeamMatchLow", par, event);
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatchLowNoPar(
    const ThinSliceEvent & event,
    int signal_index) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  //upper limit
  if (event.GetTrueEndZ() > fBeamMatchLowLimit) return 1.;

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = fBeamMatchLowPar->GetValue();

  double fraction = (fBeamMatchLowUseSingleFrac ?
                     fBeamMatchLowFraction :
                     GetFractionBySample(
                         fBeamMatchLowFractions, event.GetSampleID(),
                         signal_index));
  //std::cout << event.GetSampleID() << " " << signal_index << " " << fraction <<
  //             std::endl;
  if (fraction < 0.) {
    return 1.; 
  }
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));

  /*if (fraction < 0.) {
    std::cout << "Low Fraction " << fraction << " weight " << weight <<
                 " matched? " << matched << " var " << variation << std::endl;
  }*/
  //return weight;
  return CheckAndReturn(weight, "BeamMatchLow", (*fBeamMatchLowPar), event);
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatchHigh(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  //upper limit
  if (event.GetTrueEndZ() < fBeamMatchHighLimit) return 1.;

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = par.GetValue();

  double fraction = (fBeamMatchHighUseSingleFrac ?
                     fBeamMatchHighFraction :
                     GetFractionBySample(
                         fBeamMatchHighFractions, event.GetSampleID(),
                         signal_index));
  if (fraction < 0.) {
    return 1.; 
  }
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));
  /*
  if (fraction < 0.) {
    std::cout << "High Fraction " << fraction << " weight " << weight <<
                 " matched? " << matched << " var " << variation << std::endl;
  }*/
  if (weight < 0.) {
    std::cout << variation << " " << fraction << " " << event.GetSampleID() <<
                 " " << signal_index << std::endl;
  }
  return CheckAndReturn(weight, "BeamMatchHigh", par, event);
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatchHighNoPar(
    const ThinSliceEvent & event,
    int signal_index) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  //upper limit
  if (event.GetTrueEndZ() < fBeamMatchHighLimit) return 1.;

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = fBeamMatchHighPar->GetValue();

  double fraction = (fBeamMatchHighUseSingleFrac ?
                     fBeamMatchHighFraction :
                     GetFractionBySample(
                         fBeamMatchHighFractions, event.GetSampleID(),
                         signal_index));
  if (fraction < 0.) {
    return 1.; 
  }
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));
  if (weight < 0.) {
    std::cout << variation << " " << fraction << " " << event.GetSampleID() <<
                 " " << signal_index << std::endl;
  }
  //return weight;
  //return CheckAndReturn(weight, "BeamMatchHigh", par);
  return CheckAndReturn(weight, "BeamMatchHigh", (*fBeamMatchHighPar), event);
}



double protoana::PDSPSystematics::GetSystWeight_UpstreamInt(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int upstream_ID) {
  //if (pars.find("upstream_int_weight") == pars.end())
  //  return 1.;

  return ((event.GetSampleID() == fUpstreamID) ?
          par/*s.at("upstream_int_weight")*/.GetValue() :
          1.);
}

double protoana::PDSPSystematics::GetSystWeight_BGPions(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int past_FV_ID, int decay_ID) {
  //if (pars.find("BG_pions_weight") == pars.end())
  //  return 1.;

  bool match = (event.GetSampleID() == fPastFVID ||
                event.GetSampleID() == fDecayID);
  return (match ?
          par/*s.at("BG_pions_weight")*/.GetValue() :
          1.);
}

void protoana::PDSPSystematics::SetupSyst_BeamScraper(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_scraper_weight") == pars.end()) {
    return;
  }
  fBeamScraperPar = &(pars.at("beam_scraper_weight"));
  fActiveSystematics.push_back(
      {"beam_scraper_weight", GetSystWeight_BeamScraper});
}

double protoana::PDSPSystematics::GetSystWeight_BeamScraper(
    const ThinSliceEvent & event,
    int signal_index) {
  double weight = (event.GetIsBeamScraper() ?
          fBeamScraperPar->GetValue() :
          1.);
  return CheckAndReturn(weight, "BeamScraper", (*fBeamScraperPar), event);
}

//Trim
void protoana::PDSPSystematics::SetupSyst_BoxBeam(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("box_beam_weight") == pars.end()) {
    return;
  }
  fBoxBeamRegions
      = pars.at("box_beam_weight").GetOption<
          std::vector<std::pair<double, double>>>("Regions");
  fBoxBeamFraction = pars.at("box_beam_weight").GetOption<double>("Fraction");
}

//Trim
double protoana::PDSPSystematics::GetSystWeight_BoxBeam(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int beam_cut_ID) {
  //if (pars.find("box_beam_weight") == pars.end()) return 1.;

  if (event.GetSelectionID() != fBeamCutID) return 1.;

  double startY = event.GetRecoStartY();

  bool near_box_beam = false;
  for (const auto & region : fBoxBeamRegions) {
    if (region.first < startY && startY < region.second) {
      near_box_beam = true;
      break;
    }
  }

  double variation = par/*s.at("box_beam_weight")*/.GetValue();
 // std::cout << "near/var: " << near_box_beam << "/" << variation << std::endl;
  //return (near_box_beam ? variation :
  //        (1. - variation*fBoxBeamFraction)/(1. - fBoxBeamFraction));
  return (near_box_beam ? variation : 1.);
}

void protoana::PDSPSystematics::SetupSyst_TrueBeamShift(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("true_beam_shift") == pars.end()) {return;}

  fTrueBeamBins = pars.at("true_beam_shift").GetOption<std::vector<double>>(
      "Bins");

  TFile file(
      pars.at("true_beam_shift").GetOption<std::string>("File").c_str(), "open");
  for (size_t i = 1; i < fTrueBeamBins.size(); ++i) {
    std::string spline_name = "spline_" + std::to_string(i);
    fTrueBeamSplines.push_back((TSpline3*)file.Get(spline_name.c_str()));
  }
}

double protoana::PDSPSystematics::GetSystWeight_TrueBeamShift(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars) {
  //if (pars.find("true_beam_shift") == pars.end()) return 1.;
  if (event.GetPDG() != 211) return 1.;

  double true_startP = event.GetTrueStartP();
  int bin = 0;
  for (size_t i = 1; i < fTrueBeamBins.size(); ++i) {
    if (fTrueBeamBins[i-1] <= true_startP && true_startP <= fTrueBeamBins[i]) {
      bin = i-1;
      break;
    }
  }

  TSpline3 * spline = fTrueBeamSplines[bin];
  return spline->Eval(par/*s.at("true_beam_shift")*/.GetValue());
}

void protoana::PDSPSystematics::SetupSyst_QuadBeamShift(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("quad_beam_shift") == pars.end()) {return;}

  fQuadBeamBins = pars.at("quad_beam_shift").GetOption<std::vector<double>>(
      "Bins");

  fQuadBeamMeans = pars.at("quad_beam_shift").GetOption<std::vector<double>>(
      "Means");
  fQuadBeamSigmas = pars.at("quad_beam_shift").GetOption<std::vector<double>>(
      "Sigmas");

}

double protoana::PDSPSystematics::GetSystWeight_QuadBeamShift(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
  if (event.GetPDG() != 211 || event.GetTrueIncidentEnergies().size() == 0)
    return 1.;

  double beam_inst_KE = sqrt(std::pow(event.GetBeamInstP()*1.e3, 2) +
                             139.57*139.57) - 139.57;

  int bin = -1;
  for (size_t i = 1; i < fQuadBeamBins.size(); ++i) {
    if (fQuadBeamBins[i-1] <= beam_inst_KE &&
        beam_inst_KE <= fQuadBeamBins[i]) {
      bin = i-1;
      break;
    }
  }

  //std::cout << beam_inst_KE << " " << bin << " returning 1." << std::endl;
  if (bin == -1) {
    return 1.;
  }

  double x = (beam_inst_KE - event.GetTrueIncidentEnergies()[0]);

  double mean = fQuadBeamMeans[bin];
  double stddev = fQuadBeamSigmas[bin];
  double weight = (exp(.5*((std::pow(x - mean, 2)) -
                   (std::pow(x - par.GetValue()*mean, 2)))/(stddev*stddev)));
  //std::cout << "Returning " << par.GetValue() << " " << weight << std::endl;
  return weight;
}

double protoana::PDSPSystematics::GetSystWeight_ELoss(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID) {
  //if (pars.find("eloss_weight") == pars.end()) return 1.;
  if (event.GetSampleID() == fUpstreamID) return 1.;
  if (abs(event.GetPDG()) == 13) return 1.;

  double var = par/*s.at("eloss_weight")*/.GetValue();
  //double frac = fELossFractions[event.GetSampleID()];
  double frac = (signal_index > -1 ?
                 fELossFractions[event.GetSampleID()][signal_index] :
                 fELossFractions[event.GetSampleID()].back());

  double weight = ((event.GetDeltaEToTPC() > fELossCut) ?
                   var : (1. - var*frac)/(1. - frac));
  return CheckAndReturn(weight, "ELoss", par, event);
}

void protoana::PDSPSystematics::SetupSyst_ELoss(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("eloss_weight") == pars.end()) {return;}
  fELossCut = pars.at("eloss_weight").GetOption<double>("Cut");
  /*
  auto temp = pars.at("eloss_weight")
      .GetOption<std::vector<std::pair<int, double>>>("Fractions");
  fELossFractions = std::map<int, double>(temp.begin(), temp.end());
  */

  auto temp = pars.at("eloss_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fELossFractions = std::map<int, std::vector<double>>(temp.begin(), temp.end());
  for (auto it = fELossFractions.begin(); it != fELossFractions.end(); ++it) {
    std::cout << it->first << " ";
    for (auto f : it->second) std::cout << f << " ";
    std::cout << std::endl;
  }
}

double protoana::PDSPSystematics::GetSystWeight_ELossMuon(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID) {
  //if (pars.find("eloss_muon_weight") == pars.end()) return 1.;
  if (event.GetSampleID() == fUpstreamID) return 1.;
  if (abs(event.GetPDG()) != 13) return 1.;

  double var = par/*s.at("eloss_muon_weight")*/.GetValue();
  double frac = 1.;//fELossFractions[event.GetSampleID()];

  return ((event.GetDeltaEToTPC() > fELossCut) ?
          var : (1. - var*frac)/(1. - frac));
}

void protoana::PDSPSystematics::SetupSyst_ELossMuon(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("eloss_muon_weight") == pars.end()) {return;}
  fELossMuonCut = pars.at("eloss_muon_weight").GetOption<double>("Cut");
  auto temp = pars.at("eloss_muon_weight")
      .GetOption<std::vector<std::pair<int, double>>>("Fractions");
  fELossMuonFractions = std::map<int, double>(temp.begin(), temp.end());
}

double protoana::PDSPSystematics::CheckAndReturn(double weight,
                                                 std::string name,
                                                 const ThinSliceSystematic & par,
                                                 const ThinSliceEvent & event) {
  if (weight < 0.) {
    std::cout << "Returning negative weight " << weight << " for syst " <<
                 name << std::endl;
    std::cout << "\t" << par.GetValue() << " " << par.GetName() << std::endl;
    std::cout << "Event: " << event.GetRunID() << " " << event.GetSubrunID() <<
                 " " << event.GetEventID() << std::endl;
  }
  return weight;
}

/*
double protoana::PDSPSystematics::CheckAndReturn(double weight,
                                                 std::string name,
                                                 ThinSliceSystematic * par,
                                                 const ThinSliceEvent & event) {
  if (weight < 0.) {
    std::cout << "Returning negative weight " << weight << " for syst " <<
                 name << std::endl;
    std::cout << "\t" << par->GetValue() << " " << par->GetName() << std::endl;
    std::cout << "Event: " << event.GetRunID() << " " << event.GetSubrunID() <<
                 " " << event.GetEventID() << std::endl;
  }
  return weight;
}*/

double protoana::PDSPSystematics::GetFractionBySample(
    const std::map<int, std::vector<double>> & fractions, int sample_ID,
    int signal_index) {
  return (signal_index > -1 ?
          fractions.at(sample_ID)[signal_index] :
          fractions.at(sample_ID).back());
}
