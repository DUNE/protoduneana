#include "KaonExcReweighter.hh"
#include "geant4reweight/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/ReweightBase/G4ReweightStep.hh"
#include "Geant4/G4KaonPlus.hh"

KaonExcReweighter::KaonExcReweighter(
    TFile * FSInput,
    const std::map<std::string, TH1D*> &FSScales,
    const fhicl::ParameterSet & material_pars,
    G4ReweightManager * rw_manager,
    TH1D * inputElasticBiasHist, bool fix)
  : G4Reweighter(FSInput, FSScales, material_pars, rw_manager,
                 {"cex", "abs", "prod", "inel", "dcex"},
                 inputElasticBiasHist, fix) {
  part_def = kplus->Definition();
  fInelastic = "kaon+Inelastic";
  std::cout << "Part def: " << part_def << std::endl;
  SetupProcesses();
}

std::string KaonExcReweighter::GetInteractionSubtype(
    const G4ReweightTraj & theTraj) {

  int nK0     = (theTraj.HasChild(130).size()
                 + theTraj.HasChild(310).size());
  int nKPlus  = theTraj.HasChild(321).size();
  int nKMinus = theTraj.HasChild(-321).size();

  if( (nK0 + nKPlus + nKMinus) == 0){
    return "abs";
  }
  else if( (nKPlus + nKMinus) == 0 && nK0 == 1 ){
    return "cex";
  }
  else if( (nKPlus + nKMinus + nK0) > 1 ){
    return "prod";
  }
  else if( (nK0 + nKMinus) == 0 && nKPlus == 1 ){
    return "inel";
  }
  else if( (nK0 + nKPlus) == 0 && nKMinus == 1 ){
    return "dcex";
  }

  return "";
}

KaonExcReweighter::~KaonExcReweighter(){}
