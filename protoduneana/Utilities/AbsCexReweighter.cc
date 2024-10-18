#include "AbsCexReweighter.hh"
#include "geant4reweight/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/ReweightBase/G4ReweightStep.hh"
#include "Geant4/G4PionPlus.hh"

AbsCexReweighter::AbsCexReweighter(
    TFile * FSInput,
    const std::map<std::string, TH1D*> &FSScales,
    const fhicl::ParameterSet & material_pars,
    G4ReweightManager * rw_manager,
    TH1D * inputElasticBiasHist, bool fix)
  : G4Reweighter(FSInput, FSScales, material_pars, rw_manager,
                 {"cex", "abs", "other"},
                 inputElasticBiasHist, fix) {
  part_def = piplus->Definition();
  fInelastic = "pi+Inelastic";
  std::cout << "Part def: " << part_def << std::endl;
  SetupProcesses();
}

std::string AbsCexReweighter::GetInteractionSubtype(
    const G4ReweightTraj & theTraj) {

  int nPi0 = 0;
  int nPiCharged = 0;
  const auto & children = theTraj.GetChildren();
  //std::cout << "Traj children:" << std::endl;
  for (auto * child : children) {
    int pdg = child->GetPDG();
    //std::cout << pdg << std::endl;
    if (pdg == 111) {
      //std::cout << "Has 111" << std::endl;
      ++nPi0;
    }
    else if (abs(pdg) == 211) {

      double momentum = DBL_MAX;
      if (child->GetNSteps() > 0) {
        const auto & step = child->GetStep(0);
        momentum = sqrt(step.GetPreStepPx()*step.GetPreStepPx() +
                        step.GetPreStepPy()*step.GetPreStepPy() +
                        step.GetPreStepPz()*step.GetPreStepPz());
      }
      //std::cout << pdg << " " << momentum << std::endl;
      if (momentum > fThreshold) {
        //std::cout << "Above thresh" << std::endl;
        ++nPiCharged;
      }
    }
    else continue;
  }

  if( (nPi0 + nPiCharged) == 0){
    //std::cout << "abs" << std::endl;
    return "abs";
  }
  else if( (nPiCharged) == 0 && nPi0 > 0 ){
    //std::cout << "cex" << std::endl;
    return "cex";
  }

  //std::cout << "other" << std::endl;
  return "other";
}

AbsCexReweighter::~AbsCexReweighter(){}
