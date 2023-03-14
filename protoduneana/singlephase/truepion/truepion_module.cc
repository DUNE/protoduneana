////////////////////////////////////////////////////////////////////////
// Class:       truepion
// File:        truepion_module.cc
//
// Extract protoDUNE useful information, do a first pre-selection and save output to a flat tree
// 
// Some parts are copied from the beam module example
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
// Heng-Ye Liao modified for his protonanalysis - liao@phys.ksu.edu
//Ajib Paudel Modified for the pion truth cross-section studies
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
//#include "dunesim/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
//#include "dunesim/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TH2.h"
// C++ Includes
#include <stdio.h>
#include <stdlib.h> 
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "TComplex.h"
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include <iterator>
#include <map>

// Maximum number of beam particles to save
const int NMAXDAUGTHERS = 30;
double  m_pi=0.1395;//GeV/c^2
double prim_energy=0.0;
//double geant_energy=0;
namespace protoana {
  class truepion;
}


class protoana::truepion : public art::EDAnalyzer {
public:

  explicit truepion(fhicl::ParameterSet const & p);

  truepion(truepion const &) = delete;
  truepion(truepion &&) = delete;
  truepion & operator = (truepion const &) = delete;
  truepion & operator = (truepion &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & evt) override;

private:

  // Helper utility functions
  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNETruthUtils truthUtil;

  // Initialise tree variables
  void Initialise();

  // Fill cosmics tree
  void FillCosmicsTree(art::Event const & evt, std::string pfParticleTag);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  bool fVerbose;

  geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

  TTree *fPandoraBeam;
  //TTree *fPandoraCosmics;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
 

  double fbeamtrackMomentum;
  double fbeamtrackP[3]; //Px/Py/Pz
  double fbeamtrackEnergy;
  double fbeamtrackPos[3];
  double fbeamtrackDir[3];
  double fbeamtrackTime;
  int fbeamtrackPdg;
  int fbeamtrackID;

  //int NumberBeamTrajectoryPoints;
  std::vector<float> beamtrk_x;
  std::vector<float> beamtrk_y;
  std::vector<float> beamtrk_z;
  std::vector<float> beamtrk_Px;
  std::vector<float> beamtrk_Py;
  std::vector<float> beamtrk_Pz;
  std::vector<float> beamtrk_Eng;



  // Reconstructed tracks/showers information
  double fvertex[3];
  double fsecvertex[3];

  int fisprimarytrack;
  int fisprimaryshower;
  double fprimaryBDTScore;
  int fprimaryNHits;
  double fprimaryTheta;
  double fprimaryPhi;
  double fprimaryLength;
  double fprimaryMomentum;
  double fprimaryEndMomentum;
  double fprimaryEndPosition[3];
  double fprimaryStartPosition[3];
  double fprimaryEndDirection[3];
  double fprimaryStartDirection[3];
  double fprimaryOpeningAngle;
  int fprimaryShowerBestPlane;
  double fprimaryShowerEnergy;
  double fprimaryShowerMIPEnergy;
  double fprimaryShowerdEdx;
  double fprimaryMomentumByRangeProton;
  double fprimaryMomentumByRangeMuon;
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  int fprimaryID;




  //*********************truth information*******************************************//
  int fprimary_truth_TrackId;
  int fprimary_truth_Pdg;
  int ftruthpdg;
  double fprimary_truth_StartPosition[4];
  double fprimary_truth_EndPosition[4];
  std::string fprimary_truth_EndProcess;
  std::string truth_last_process;
  double fprimary_truth_P;
  double fprimary_truth_Momentum[4];
  double fprimary_truth_EndMomentum[4];
  double fprimary_truth_Pt;
  double fprimary_truth_Mass;
  double fprimary_truth_Theta;
  double fprimary_truth_Phi;
  int fprimary_truth_Process;
  int fprimary_truth_Isbeammatched;
  int fprimary_truth_NDaughters;
  double fprimary_truth_tracklength;
  

  //interaction point details
  std::vector<double> interactionX;
  std::vector<double> interactionY;
  std::vector<double> interactionZ;
  std::vector<std::string> interactionProcesslist;
  std::vector<double> interactionAngles;
  std::vector<double> Zintersection;
  std::vector<double> timeintersection;

  //calo info
  std::vector< std::vector<double> > primtrk_dqdx;
  std::vector< std::vector<double> > primtrk_resrange;
  std::vector< std::vector<double> > primtrk_dedx;
  std::vector<double> primtrk_range;
  std::vector< std::vector<double> > primtrk_hitx;
  std::vector< std::vector<double> > primtrk_hity;
  std::vector< std::vector<double> > primtrk_hitz;
  std::vector< std::vector<double> > primtrk_pitch;
  std::vector<std::vector<double> > primtrk_truth_Z;
  std::vector<std::vector<double> > primtrk_truth_Eng;
  std::vector<std::vector<double> > primtrk_truth_trkide;



  // Daughters from primary
  int fNDAUGHTERS;
  int fisdaughtertrack[NMAXDAUGTHERS];
  int fisdaughtershower[NMAXDAUGTHERS];
  int fdaughterNHits[NMAXDAUGTHERS];
  double fdaughterTheta[NMAXDAUGTHERS];
  double fdaughterPhi[NMAXDAUGTHERS];
  double fdaughterLength[NMAXDAUGTHERS];
  double fdaughterMomentum[NMAXDAUGTHERS];
  double fdaughterEndMomentum[NMAXDAUGTHERS];
  double fdaughterEndPosition[NMAXDAUGTHERS][3];
  double fdaughterStartPosition[NMAXDAUGTHERS][3];
  double fdaughterEndDirection[NMAXDAUGTHERS][3];
  double fdaughterStartDirection[NMAXDAUGTHERS][3];
  double fdaughterOpeningAngle[NMAXDAUGTHERS];
  double fdaughterShowerEnergy[NMAXDAUGTHERS];
  double fdaughterShowerMIPEnergy[NMAXDAUGTHERS];
  double fdaughterShowerdEdx[NMAXDAUGTHERS];
  int fdaughterShowerBestPlane[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeProton[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeMuon[NMAXDAUGTHERS];
  double fdaughterKineticEnergy[NMAXDAUGTHERS][3];
  double fdaughterRange[NMAXDAUGTHERS][3];
  int fdaughterID[NMAXDAUGTHERS];
  //double fdaughterT0[NMAXDAUGTHERS];

  int fdaughter_truth_TrackId[NMAXDAUGTHERS];
  int fdaughter_truth_Pdg[NMAXDAUGTHERS];
  double fdaughter_truth_StartPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_P[NMAXDAUGTHERS];
  double fdaughter_truth_Momentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndMomentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_Pt[NMAXDAUGTHERS];
  double fdaughter_truth_Mass[NMAXDAUGTHERS];
  double fdaughter_truth_Theta[NMAXDAUGTHERS];
  double fdaughter_truth_Phi[NMAXDAUGTHERS];
  int fdaughter_truth_Process[NMAXDAUGTHERS];

  double minX =  -360.0;
  double maxX = 360.0;
  double minY =0.0;
  double maxY = 600.0;
  double minZ =  0.0;
  double maxZ = 695.0;

};


protoana::truepion::truepion(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose"))
{

}

void protoana::truepion::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fPandoraBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
  fPandoraBeam->Branch("run",                           &fRun,                          "run/I");
  fPandoraBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
  fPandoraBeam->Branch("event",                         &fevent,                        "event/I");
  fPandoraBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
  fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
  fPandoraBeam->Branch("beamtrackP",                    &fbeamtrackP,                   "beamtrackP[3]/D");
  fPandoraBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
  fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
  fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
  fPandoraBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
  fPandoraBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
  fPandoraBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");

  //fPandoraBeam->Branch("NumberBeamTrajectoryPoints",    &NumberBeamTrajectoryPoints,    "NumberBeamTrajectoryPoints/I");
  fPandoraBeam->Branch("beamtrk_x",&beamtrk_x);
  fPandoraBeam->Branch("beamtrk_y",&beamtrk_y);
  fPandoraBeam->Branch("beamtrk_z",&beamtrk_z);
  fPandoraBeam->Branch("beamtrk_Px",&beamtrk_Px);
  fPandoraBeam->Branch("beamtrk_Py",&beamtrk_Py);
  fPandoraBeam->Branch("beamtrk_Pz",&beamtrk_Pz);
  fPandoraBeam->Branch("beamtrk_Eng",&beamtrk_Eng);


  fPandoraBeam->Branch("vertex",                        &fvertex,                       "vertex[3]/D");
  fPandoraBeam->Branch("secvertex",                     &fsecvertex,                    "secvertex[3]/D");
  fPandoraBeam->Branch("isprimarytrack",                &fisprimarytrack,               "isprimarytrack/I");
  fPandoraBeam->Branch("isprimaryshower",               &fisprimaryshower,              "isprimaryshower/I");
  fPandoraBeam->Branch("primaryBDTScore",               &fprimaryBDTScore,              "primaryBDTScore/D");
  fPandoraBeam->Branch("primaryNHits",                  &fprimaryNHits,                 "fprimaryNHits/I");
  fPandoraBeam->Branch("primaryTheta",                  &fprimaryTheta,                 "primaryTheta/D");
  fPandoraBeam->Branch("primaryPhi",                    &fprimaryPhi,                   "primaryPhi/D");
  fPandoraBeam->Branch("primaryLength",                 &fprimaryLength,                "primaryLength/D");
  fPandoraBeam->Branch("primaryMomentum",               &fprimaryMomentum,              "primaryMomentum/D");
  fPandoraBeam->Branch("primaryEndMomentum",            &fprimaryEndMomentum,           "primaryEndMomentum/D");
  fPandoraBeam->Branch("primaryEndPosition",            &fprimaryEndPosition,           "primaryEndPosition[3]/D");
  fPandoraBeam->Branch("primaryStartPosition",          &fprimaryStartPosition,         "primaryStartPosition[3]/D");
  fPandoraBeam->Branch("primaryEndDirection",           &fprimaryEndDirection,          "primaryEndDirection[3]/D");
  fPandoraBeam->Branch("primaryStartDirection",         &fprimaryStartDirection,        "primaryStartDirection[3]/D");
  fPandoraBeam->Branch("primaryOpeningAngle",           &fprimaryOpeningAngle,          "primaryOpeningAngle/D");
  fPandoraBeam->Branch("primaryID",                     &fprimaryID,                    "primaryID/I");
  fPandoraBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
  fPandoraBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/I");
  fPandoraBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/I");
  fPandoraBeam->Branch("primaryShowerdEdx",             &fprimaryShowerdEdx,            "primaryShowerdEdx/I");
  fPandoraBeam->Branch("primaryMomentumByRangeProton",  &fprimaryMomentumByRangeProton, "primaryMomentumByRangeProton/D");
  fPandoraBeam->Branch("primaryMomentumByRangeMuon",    &fprimaryMomentumByRangeMuon,   "primaryMomentumByRangeMuon/D");
  fPandoraBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
  fPandoraBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
 

  fPandoraBeam->Branch("primary_truth_TrackId",         &fprimary_truth_TrackId,         "primary_truth_TrackId/I");
  fPandoraBeam->Branch("primary_truth_Pdg",             &fprimary_truth_Pdg,             "primary_truth_Pdg/I");
  fPandoraBeam->Branch("truthpdg",                      &ftruthpdg,                      "truthpdg/I");


  fPandoraBeam->Branch("primary_truth_StartPosition",   &fprimary_truth_StartPosition,   "primary_truth_StartPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_EndPosition",     &fprimary_truth_EndPosition,     "primary_truth_EndPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_Momentum",        &fprimary_truth_Momentum,        "primary_truth_Momentum[4]/D");
  fPandoraBeam->Branch("primary_truth_EndMomentum",     &fprimary_truth_EndMomentum,     "primary_truth_EndMomentum[4]/D");
  fPandoraBeam->Branch("primary_truth_P",               &fprimary_truth_P,               "primary_truth_P/D");
  fPandoraBeam->Branch("primary_truth_Pt",              &fprimary_truth_Pt,              "primary_truth_Pt/D");
  fPandoraBeam->Branch("primary_truth_Mass",            &fprimary_truth_Mass,            "primary_truth_Mass/D");
  fPandoraBeam->Branch("primary_truth_Theta",           &fprimary_truth_Theta,           "primary_truth_Theta/D");
  fPandoraBeam->Branch("primary_truth_Phi",             &fprimary_truth_Phi,             "primary_truth_Phi/D");
  fPandoraBeam->Branch("primary_truth_Process",         &fprimary_truth_Process,         "primary_truth_Process/I");
  fPandoraBeam->Branch("primary_truth_Isbeammatched",   &fprimary_truth_Isbeammatched,   "primary_truth_Isbeammatched/I");
  fPandoraBeam->Branch("primary_truth_NDaughters",      &fprimary_truth_NDaughters,      "primary_truth_NDaughters/I");
  fPandoraBeam->Branch("primary_truth_EndProcess",      &fprimary_truth_EndProcess);
  fPandoraBeam->Branch("truth_last_process",            &truth_last_process);
  fPandoraBeam->Branch("primary_truth_tracklength",     &fprimary_truth_tracklength,      "primary_truth_tracklength/D");

  fPandoraBeam->Branch("interactionX",&interactionX);
  fPandoraBeam->Branch("interactionY",&interactionY);
  fPandoraBeam->Branch("interactionZ",&interactionZ);
  fPandoraBeam->Branch("interactionProcesslist",&interactionProcesslist);
  fPandoraBeam->Branch("interactionAngles",&interactionAngles);
  fPandoraBeam->Branch("Zintersection",&Zintersection);
  fPandoraBeam->Branch("timeintersection",&timeintersection);


  fPandoraBeam->Branch("NDAUGHTERS",                    &fNDAUGHTERS,                   "NDAUGHTERS/I");
  fPandoraBeam->Branch("isdaughtertrack",               &fisdaughtertrack,              "isdaughtertrack[NDAUGHTERS]/I");
  fPandoraBeam->Branch("isdaughtershower",              &fisdaughtershower,             "isdaughtershower[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterNHits",                 &fdaughterNHits,                "daughterNHits[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterTheta",                 &fdaughterTheta,                "daughterTheta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterPhi",                   &fdaughterPhi,                  "daughterPhi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterLength",                &fdaughterLength,               "daughterLength[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentum",              &fdaughterMomentum,             "daughterMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndMomentum",           &fdaughterEndMomentum,          "daughterEndMomentum[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndPosition",           &fdaughterEndPosition,          "daughterEndPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartPosition",         &fdaughterStartPosition,        "daughterStartPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartDirection",        &fdaughterStartDirection,       "daughterStartDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterEndDirection",          &fdaughterEndDirection,         "daughterEndDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterOpeningAngle",          &fdaughterOpeningAngle,         "daughterOpeningAngle[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerBestPlane",       &fdaughterShowerBestPlane,      "daughterShowerBestPlane[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerEnergy",          &fdaughterShowerEnergy,         "daughterShowerEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerMIPEnergy",       &fdaughterShowerMIPEnergy,      "daughterShowerMIPEnergy[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerdEdx",            &fdaughterShowerdEdx,           "daughterShowerdEdx[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeProton", &fdaughterMomentumByRangeProton,"daughterMomentumByRangeProton[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterMomentumByRangeMuon",   &fdaughterMomentumByRangeMuon,  "daughterMomentumByRangeMuon[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterKineticEnergy",         &fdaughterKineticEnergy,        "daughterKineticEnergy[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterRange",                 &fdaughterRange,                "daughterRange[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterID",                    &fdaughterID,                   "daughterID[NDAUGHTERS]/I");
  // fPandoraBeam->Branch("daughterT0",                    &fdaughterT0,                   "daughterT0[NDAUGHTERS]/D");

  fPandoraBeam->Branch("daughter_truth_TrackId",        &fdaughter_truth_TrackId,       "daughter_truth_TrackId[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_Pdg",            &fdaughter_truth_Pdg,           "daughter_truth_Pdg[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughter_truth_StartPosition",  &fdaughter_truth_StartPosition, "daughter_truth_StartPosition[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_EndPosition",    &fdaughter_truth_EndPosition,   "daughter_truth_EndPosition[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_Momentum",       &fdaughter_truth_Momentum,      "daughter_truth_Momentum[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_EndMomentum",    &fdaughter_truth_EndMomentum,   "daughter_truth_EndMomentum[NDAUGHTERS][4]/D");
  fPandoraBeam->Branch("daughter_truth_P",              &fdaughter_truth_P,             "daughter_truth_P[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Pt",             &fdaughter_truth_Pt,            "daughter_truth_Pt[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Mass",           &fdaughter_truth_Mass,          "daughter_truth_Mass[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Theta",          &fdaughter_truth_Theta,         "daughter_truth_Theta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Phi",            &fdaughter_truth_Phi,           "daughter_truth_Phi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughter_truth_Process",        &fdaughter_truth_Process,       "daughter_truth_Process[NDAUGHTERS]/I");
 

  // hKE_truth_reco = tfs->make<TH2D>("hKE_truth_reco","KE truth vs reco;true KE in MeV;reconstructed KE in MeV"  , 3000, -100, 1400  ,3000, -100, 1400); 
  // hEndZ_truth_reco = tfs->make<TH2D>("hEndZ_truth_reco","EndZ truth vs reco;true EndZ (cm);reconstructed EndZ(cm)"  , 695, 0, 695  ,695,0,695); 

  fPandoraBeam->Branch("primtrk_dqdx",&primtrk_dqdx);
  fPandoraBeam->Branch("primtrk_dedx",&primtrk_dedx);
  fPandoraBeam->Branch("primtrk_resrange",&primtrk_resrange);
  fPandoraBeam->Branch("primtrk_range",&primtrk_range);
  fPandoraBeam->Branch("primtrk_hitx",&primtrk_hitx);
  fPandoraBeam->Branch("primtrk_hity",&primtrk_hity);
  fPandoraBeam->Branch("primtrk_hitz",&primtrk_hitz);
  fPandoraBeam->Branch("primtrk_pitch",&primtrk_pitch);
  ///////////////////////////////////////////////////
  fPandoraBeam->Branch("primtrk_truth_Z",&primtrk_truth_Z);//primary track true Z positions 
  fPandoraBeam->Branch("primtrk_truth_Eng",&primtrk_truth_Eng);//primary track true Energy deposited for each Z position
  fPandoraBeam->Branch("primtrk_truth_trkide",&primtrk_truth_trkide);//primary track true Energy deposited for each Z position
}//begin job 

void protoana::truepion::analyze(art::Event const & evt){

  // Initialise tree parameters
  Initialise();

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  // const sim::ParticleList& plist=pi_serv->ParticleList();

  int beamid=-9999;
  int truthid=-999;

  fRun = evt.run();
  fSubRun = evt.subRun();
  fevent  = evt.id().event(); 
  art::Timestamp ts = evt.time();
  if (ts.timeHigh() == 0){
    TTimeStamp ts2(ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }
  else{
    TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
    fTimeStamp = ts2.AsDouble();
  }

  bool beamTriggerEvent = false;
  // If this event is MC then we can check what the true beam particle is
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  if(!evt.isRealData()){
    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    //std::cout<<"geantGoodParticle "<<geantGoodParticle.size()<<std::endl;
    if(geantGoodParticle != 0x0){
      std::cout<<"geant good particle loop "<<std::endl;
      std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() 
		<< " , track id = " << geantGoodParticle->TrackId()
		<< " , Vx/Vy/Vz = " << geantGoodParticle->Vx() << "/"<< geantGoodParticle->Vy() << "/" << geantGoodParticle->Vz() 

		<< std::endl;

      std::vector<double> tmp_primtrk_truth_Z;
      std::vector<double> tmp_primtrk_truth_Eng;
      std::vector<double> tmp_primtrk_truth_trkide;     
      
      beamTriggerEvent = true;
      fbeamtrackPos[0]   = geantGoodParticle->Vx();
      fbeamtrackPos[1]   = geantGoodParticle->Vy();
      fbeamtrackPos[2]   = geantGoodParticle->Vz();
      fbeamtrackMomentum = geantGoodParticle->P();
      fbeamtrackP[0]     = geantGoodParticle->Px();
      fbeamtrackP[1]     = geantGoodParticle->Py();
      fbeamtrackP[2]     = geantGoodParticle->Pz();
      // fbeamtrackEnergy   = geantGoodParticle->E();
      fbeamtrackPdg      = geantGoodParticle->PdgCode();
      fbeamtrackTime     = geantGoodParticle->T();
      fbeamtrackID       = geantGoodParticle->TrackId();
      fprimary_truth_TrackId          = geantGoodParticle->TrackId();
      fprimary_truth_Pdg              = geantGoodParticle->PdgCode();
      beamid                          = geantGoodParticle->TrackId();
      fprimary_truth_StartPosition[3] = geantGoodParticle->T();
      fprimary_truth_EndPosition[3]   = geantGoodParticle->EndT();
      fprimary_truth_P                = geantGoodParticle->P();
      fprimary_truth_Momentum[3]      = geantGoodParticle->E();
      fprimary_truth_Pt               = geantGoodParticle->Pt();
      fprimary_truth_Mass             = geantGoodParticle->Mass();
      fprimary_truth_EndMomentum[0]   = geantGoodParticle->EndPx();
      fprimary_truth_EndMomentum[1]   = geantGoodParticle->EndPy();
      fprimary_truth_EndMomentum[2]   = geantGoodParticle->EndPz();
      fprimary_truth_EndMomentum[3]   = geantGoodParticle->EndE();
      fprimary_truth_Theta            = geantGoodParticle->Momentum().Theta();
      fprimary_truth_Phi              = geantGoodParticle->Momentum().Phi();
      fprimary_truth_NDaughters       = geantGoodParticle->NumberDaughters();

      ////////////////////////////////////////
      truth_last_process=geantGoodParticle->EndProcess();
      ////here is another new parameter added

      fprimary_truth_Process = int(geantGoodParticle->Trajectory().ProcessToKey(geantGoodParticle->Process()));
      prim_energy=0;
      for(size_t i_s=0; i_s < geantGoodParticle->NumberTrajectoryPoints(); i_s++){ //loop over beam tracks
	//	if(geantGoodParticle->Position(i_s).Z()>0) break;
	beamtrk_x.push_back(geantGoodParticle->Position(i_s).X());
	beamtrk_y.push_back(geantGoodParticle->Position(i_s).Y());
	beamtrk_z.push_back(geantGoodParticle->Position(i_s).Z());

        beamtrk_Px.push_back(geantGoodParticle->Momentum(i_s).X());
        beamtrk_Py.push_back(geantGoodParticle->Momentum(i_s).Y());
        beamtrk_Pz.push_back(geantGoodParticle->Momentum(i_s).Z());
	
        beamtrk_Eng.push_back(geantGoodParticle->Momentum(i_s).E()-geantGoodParticle->Mass());
	if(geantGoodParticle->Position(i_s).Z()<0) prim_energy=1000*(geantGoodParticle->Momentum(i_s).E()-geantGoodParticle->Mass());
	if(geantGoodParticle->Position(i_s).Z()<0) fbeamtrackEnergy=prim_energy; //correct energy at the beginning of track
      } //loop over beam trks

      //new section 
      art::ServiceHandle<cheat::BackTrackerService> bt_serv;


      simb::MCTrajectory truetraj=geantGoodParticle->Trajectory();
      auto thisTrajectoryProcessMap1 =  truetraj.TrajectoryProcesses();
      if (thisTrajectoryProcessMap1.size()){
	for(auto const& couple: thisTrajectoryProcessMap1){
	  // int_label=truetraj.KeyToProcess(couple.second);
	  fprimary_truth_EndPosition[0]=((truetraj.at(couple.first)).first).X();
	  fprimary_truth_EndPosition[1]=((truetraj.at(couple.first)).first).Y();
	  fprimary_truth_EndPosition[2]=((truetraj.at(couple.first)).first).Z();
	  fprimary_truth_EndProcess=truetraj.KeyToProcess(couple.second);
	  fprimary_truth_Momentum[0]=((truetraj.at(couple.first)).second).X();
	  fprimary_truth_Momentum[1]= ((truetraj.at(couple.first)).second).Y();
	  fprimary_truth_Momentum[2]=((truetraj.at(couple.first)).second).Z();
	  break;
	}
      }

      ////saving the complete information of all the interactions
      std::cout<<"interaction map size "<<thisTrajectoryProcessMap1.size()<<std::endl;
      if (thisTrajectoryProcessMap1.size()){
	for(auto const& couple1: thisTrajectoryProcessMap1){
	 
	  if ((truetraj.KeyToProcess(couple1.second)).find("CoulombScat")!= std::string::npos) continue;
	    // Let's check if the interaction is in the the TPC
	  auto     interactionPos4D =  (truetraj.at(couple1.first)).first ;        
	  if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
	  else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
	  else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;
	  interactionX.push_back(((truetraj.at(couple1.first)).first).X());
	  interactionY.push_back(((truetraj.at(couple1.first)).first).Y());
	  interactionZ.push_back(((truetraj.at(couple1.first)).first).Z());
	  interactionProcesslist.push_back(truetraj.KeyToProcess(couple1.second));
	  std::cout<<"number of interactions "<<thisTrajectoryProcessMap1.size()<<std::endl;
	  std::cout<<"int X, Y, Z and process "<<((truetraj.at(couple1.first)).first).X()<<" "<<((truetraj.at(couple1.first)).first).Y()<<" "<<((truetraj.at(couple1.first)).first).Z()<<" "<<truetraj.KeyToProcess(couple1.second)<<std::endl;
	  ///get the interaction angle here
	  double interactionAngle = 999999.; // This needs to be changed
	  //--------------------- Int Angle ---------------------------
	  // Try to retreive the interaction angle
	  auto  prevInteractionPos4D = (truetraj.at(couple1.first-1)).first ;
	  auto  prevInteractionPos3D = prevInteractionPos4D.Vect() ;
	  auto  interactionPos3D     = interactionPos4D.Vect() ;
	  auto  distanceBtwPoint     = interactionPos3D - prevInteractionPos3D;
	  //Let's try to see if the next point exists
	  if (truetraj.size() > couple1.first + 1) {
	    // The particle doesn't die. No need to check for anything else.
	    auto nextInteractionPos4D =  (truetraj.at(couple1.first+1)).first ;
	    auto nextInteractionPos3D =  nextInteractionPos4D.Vect() ;
	    auto distanceBtwPointNext =  nextInteractionPos3D - interactionPos3D;
	    interactionAngles.push_back(TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  ));
	    continue;
	  }
	  interactionAngles.push_back(interactionAngle);
	}
      }

      geo::View_t view = fWireReadoutGeom.Plane({0, 0, 2}).View();
      auto simIDE_prim=bt_serv->TrackIdToSimIDEs_Ps(geantGoodParticle->TrackId(),view);
      std::map<double, sim::IDE> orderedSimIDE;
      for (auto& ide : simIDE_prim) orderedSimIDE[ide->z]= *ide;
      auto inTPCPoint  = truetraj.begin(); 
      auto Momentum0   = inTPCPoint->second;
      auto old_iter = orderedSimIDE.begin();
      double tlen=0.0;
      double xi=0.0;double yi=0.0;double zi=0.0;
      int count=0; 
     
    
      for ( auto iter= orderedSimIDE.begin(); iter!= orderedSimIDE.end(); iter++,old_iter++){
	auto currentIde = iter->second;
	if(currentIde.z<minZ) continue;
	else if (currentIde.x < minX || currentIde.x > maxX ) continue;
	else if (currentIde.y < minY || currentIde.y > maxY ) continue;
	tmp_primtrk_truth_Z.push_back(currentIde.z);
	tmp_primtrk_truth_Eng.push_back(currentIde.energy);
	tmp_primtrk_truth_trkide.push_back(currentIde.trackID);
	if(count==0){
	  fprimary_truth_StartPosition[0] = currentIde.x;
	  fprimary_truth_StartPosition[1] = currentIde.y;
	  fprimary_truth_StartPosition[2] = currentIde.z;
	}
	if(currentIde.trackID>=0){
	  if(count>0){
	    tlen=tlen+TMath::Sqrt(std::pow(currentIde.x-xi,2)+std::pow(currentIde.y-yi,2)+std::pow(currentIde.z-zi,2));
	  }//if count
	  xi=currentIde.x;yi=currentIde.y;zi=currentIde.z;
	  count++;
	}//trackid>0 loop
      }// iter loop
      fprimary_truth_tracklength=tlen;
      primtrk_truth_Z.push_back(tmp_primtrk_truth_Z);
      primtrk_truth_Eng.push_back(tmp_primtrk_truth_Eng);
      primtrk_truth_trkide.push_back(tmp_primtrk_truth_trkide);

      tmp_primtrk_truth_Z.clear();
      tmp_primtrk_truth_Eng.clear();
      tmp_primtrk_truth_trkide.clear();
      //new section 

    }// geantGoodParticle
  }//is real data loop

  /*
  // Now we want to access the output from Pandora. This comes in the form of particle flow objects (recob::PFParticle).
  // The primary PFParticles are those we want to consider and these PFParticles then have a hierarchy of daughters that
  // describe the whole interaction of a given primary particle
  //
  //                     / daughter track
  //                    /
  //  primary track    /   
  //  ---------------- ---- daughter track
  //                   \
  //                   /\-
  //                   /\\-- daughter shower
  //
  // The above primary PFParticle will have links to three daughter particles, two track-like and one shower-like
  */

  // Track momentum algorithm calculates momentum based on track range
  trkf::TrackMomentumCalculator trmom;
  //trmom.SetMinLength(100);

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag) << std::endl;

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  // std::vector<recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
 


  auto pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  //cluster information

  std::vector<art::Ptr<recob::Track> > tracklist;
  auto trackListHandle = evt.getHandle< std::vector<recob::Track> >("pandoraTrack");
  if(trackListHandle) art::fill_ptr_vector(tracklist, trackListHandle);

  std::vector<art::Ptr<recob::PFParticle> > pfplist;
  auto PFPListHandle = evt.getHandle< std::vector<recob::PFParticle> >("pandora");
  if(PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);
  
  // get information about the hits
  std::vector<art::Ptr<recob::Cluster>> clusterlist;
  auto clusterListHandle = evt.getHandle< std::vector<recob::Cluster> >("pandora");
  if(clusterListHandle)
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::FindManyP<recob::Cluster> fmcp(PFPListHandle,evt,"pandora");
  art::FindManyP<recob::Track> pftrack(PFPListHandle,evt,"pandoraTrack");

    std::cout<<"number of pfp_particles "<<pfplist.size()<<std::endl;
    std::cout<<" size of pfParticles "<<pfParticles.size()<<std::endl;

    /* for(size_t p1=0;p1<pfplist.size();p1++){
      std::vector<art::Ptr<recob::Track>> trk=pftrack.at(p1);
      if(trk.size()) std::cout<<" trk  key "<<trk[0].key()<<std::endl; 

      std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(p1);
      std::cout<<"cluster size for each particle "<<allClusters.size()<<std::endl;
      for(size_t c1=0;c1<allClusters.size();c1++){
      std::cout<<" cluster ID "<<allClusters[c1]->ID();
	std::cout<<" plane number "<<allClusters[c1]->Plane().Plane;
	std::cout<<" TPC number "<<allClusters[c1]->Plane().TPC;
	std::cout<<" start wire "<<allClusters[c1]->StartWire();
	std::cout<<" end wire "<<allClusters[c1]->EndWire();
	std::cout<<" start tick "<<allClusters[c1]->StartTick();
	std::cout<<" end tick "<<allClusters[c1]->EndTick();
      }
    }*/




  // We can now look at these particles
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  for(const recob::PFParticle* particle : pfParticles){
    // Pandora's BDT beam-cosmic score
    fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
    // NHits associated with this pfParticle
    fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();
    // of this particle might be more helpful. These return null pointers if not track-like / shower-like
    const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    /////new line added here
    // std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(particle);
    // std::cout<<allClusters.size();
    /* for(size_t c1=0;c1<allClusters.size();c1++){


       }*/


    if(thisTrack != 0x0){
      // Get the true mc particle
      const simb::MCParticle* mcparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, evt, fTrackerTag);
      if(mcparticle!=0x0){
	std::cout<<"ftruth pdg "<<mcparticle->PdgCode()<<std::endl;
	ftruthpdg=mcparticle->PdgCode();
     
	truthid=mcparticle->TrackId();
	fprimary_truth_Isbeammatched=0;
	if(beamid==truthid) fprimary_truth_Isbeammatched=1;
      }//mcparticle loop





      fisprimarytrack               = 1;
      fisprimaryshower              = 0;
      fprimaryID                    = thisTrack->ID();
      std::cout<<"this Track track ID "<<thisTrack->ID()<<std::endl;
      fprimaryTheta                 = thisTrack->Theta();
      fprimaryPhi                   = thisTrack->Phi();
      fprimaryLength                = thisTrack->Length();
      fprimaryMomentum              = thisTrack->StartMomentum();
      fprimaryEndMomentum           = thisTrack->EndMomentum();
      fprimaryEndPosition[0]        = thisTrack->End().X();
      fprimaryEndPosition[1]        = thisTrack->End().Y();
      fprimaryEndPosition[2]        = thisTrack->End().Z();
      fprimaryStartPosition[0]      = thisTrack->Start().X();
      fprimaryStartPosition[1]      = thisTrack->Start().Y();
      fprimaryStartPosition[2]      = thisTrack->Start().Z();
      fprimaryEndDirection[0]       = thisTrack->Trajectory().EndDirection().X();
      fprimaryEndDirection[1]       = thisTrack->Trajectory().EndDirection().Y();
      fprimaryEndDirection[2]       = thisTrack->Trajectory().EndDirection().Z();
      fprimaryStartDirection[0]     = thisTrack->Trajectory().StartDirection().X();
      fprimaryStartDirection[1]     = thisTrack->Trajectory().StartDirection().Y();
      fprimaryStartDirection[2]     = thisTrack->Trajectory().StartDirection().Z();
      fprimaryMomentumByRangeMuon   = trmom.GetTrackMomentum(thisTrack->Length(),13);
      fprimaryMomentumByRangeProton = trmom.GetTrackMomentum(thisTrack->Length(),2212); 
      std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
      if(calovector.size() != 3)
	std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;
      std::vector<double> tmp_primtrk_dqdx;	
      std::vector<double> tmp_primtrk_resrange;	
      std::vector<double> tmp_primtrk_dedx;	
      std::vector<double> tmp_primtrk_hitx;	
      std::vector<double> tmp_primtrk_hity;	
      std::vector<double> tmp_primtrk_hitz;
      std::vector<double> tmp_primtrk_pitch;
      for (auto & calo : calovector) {
	if (calo.PlaneID().Plane == 2){ //only collection plane
	  primtrk_range.push_back(calo.Range());
	  for (size_t ihit = 0; ihit < calo.dQdx().size(); ++ihit){ //loop over hits
	    tmp_primtrk_dqdx.push_back(calo.dQdx()[ihit]);
	    tmp_primtrk_resrange.push_back(calo.ResidualRange()[ihit]);
	    tmp_primtrk_dedx.push_back(calo.dEdx()[ihit]);
	    tmp_primtrk_pitch.push_back(calo.TrkPitchVec()[ihit]);

	    const auto &primtrk_pos=(calo.XYZ())[ihit];
	    tmp_primtrk_hitx.push_back(primtrk_pos.X());
	    tmp_primtrk_hity.push_back(primtrk_pos.Y());
	    tmp_primtrk_hitz.push_back(primtrk_pos.Z());
	  } //loop over hits
	} //only collection plane
      }//calovector

      if(tmp_primtrk_dqdx.size()!=0){
	if(tmp_primtrk_hitz[0]>tmp_primtrk_hitz[tmp_primtrk_hitz.size()-1]){
	  std::reverse(tmp_primtrk_hitz.begin(),tmp_primtrk_hitz.end());
	  std::reverse(tmp_primtrk_pitch.begin(),tmp_primtrk_pitch.end());
	  std::reverse(tmp_primtrk_dedx.begin(),tmp_primtrk_dedx.end());
	  std::reverse(tmp_primtrk_dqdx.begin(),tmp_primtrk_dqdx.end());
	  std::reverse(tmp_primtrk_resrange.begin(),tmp_primtrk_resrange.end());
	}
	primtrk_dqdx.push_back(tmp_primtrk_dqdx);
	primtrk_resrange.push_back(tmp_primtrk_resrange);
	primtrk_dedx.push_back(tmp_primtrk_dedx);
	primtrk_hitx.push_back(tmp_primtrk_hitx);
	primtrk_hity.push_back(tmp_primtrk_hity);
	primtrk_hitz.push_back(tmp_primtrk_hitz);
	primtrk_pitch.push_back(tmp_primtrk_pitch);
      }
      tmp_primtrk_dqdx.clear();
      tmp_primtrk_resrange.clear();
      tmp_primtrk_dedx.clear();
      tmp_primtrk_hitx.clear();
      tmp_primtrk_hity.clear();
      tmp_primtrk_hitz.clear();
      tmp_primtrk_pitch.clear();
     
      for(size_t k = 0; k < calovector.size() && k<3; k++){
	fprimaryKineticEnergy[k] = calovector[k].KineticEnergy();
	fprimaryRange[k] = calovector[k].Range();
	//const std::vector< double > & dedxvec = calovector[k].dEdx();
      }

      //**************Section for dE/dx correction**************************//
      std::vector<float> Stw, Endw, Stt, Endt, Stwires, Endwires, Stticks, Endticks, TPCb, TPCcl;
      Stw.clear(); Endw.clear(); Stt.clear();  Endt.clear();  Stwires.clear();  Endwires.clear(); Stticks.clear(); Endticks.clear(); TPCb.clear(); TPCcl.clear();
      float den;
      float numw, numt,wire_no,ticks_no;
      for(size_t p1=0;p1<pfplist.size();p1++){
	std::vector<art::Ptr<recob::Track>> trk=pftrack.at(p1);
	std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(p1);
	for(size_t c1=0;c1<allClusters.size();c1++){
	  if(allClusters[c1]->Plane().Plane!=2) continue;
	  if(trk.size() && int(trk[0].key())==fprimaryID){
	    Stw.push_back(allClusters[c1]->StartWire());
	    Endw.push_back(allClusters[c1]->EndWire());
	    Stt.push_back(allClusters[c1]->StartTick());
	    Endt.push_back(allClusters[c1]->EndTick());
	    TPCb.push_back(allClusters[c1]->Plane().TPC);
	  }
	  else{
	    Stwires.push_back(allClusters[c1]->StartWire());
	    Endwires.push_back(allClusters[c1]->EndWire());
	    Stticks.push_back(allClusters[c1]->StartTick());
	    Endticks.push_back(allClusters[c1]->EndTick());
	    TPCcl.push_back(allClusters[c1]->Plane().TPC);
	  }
	}
      }
      for(size_t clt=0;clt<Stw.size();clt++){
	for(size_t cl1=0;cl1<Stwires.size();cl1++){
	  if(TPCcl[cl1]!=TPCb[clt]) continue;
	  std::cout<<"tpc are equal "<<std::endl;
	  den=(Stw[clt]-Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]-Endwires[cl1]);
	  if(den==0) continue;
	  numw=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stwires[cl1]-Endwires[cl1])-(Stw[clt]-Endw[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
	  numt=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
	  wire_no=numw/den;
	  ticks_no=numt/den;
	  //  std::cout<<"wireno and ticks not solution "<<wire_no<<"  "<<ticks_no<<std::endl;
	  if(((Stw[clt]<wire_no && Endw[clt]>wire_no)||(Stw[clt]>wire_no && Endw[clt]<wire_no))&&((Stt[clt]<ticks_no && Endt[clt]>ticks_no)||(Stt[clt]>ticks_no && Endt[clt]<ticks_no)) && ((Stwires[cl1]<wire_no && Endwires[cl1]>wire_no)||(Stwires[cl1]>wire_no && Endwires[cl1]<wire_no)) && ((Stticks[cl1]<ticks_no && Endticks[cl1]>ticks_no)||(Stticks[cl1]>ticks_no && Endticks[cl1]<ticks_no)))
	    { 
 std::cout<<"intersection wire and ticks are "<<std::round(wire_no)<<"  "<<ticks_no<<" Stw Endw StT EndT "<<Stwires[cl1]<<" "<<Endwires[cl1]<<" "<<Stticks[cl1]<<" "<<Endticks[cl1]<<std::endl;
 double xyzStart[3];
 double xyzEnd[3];
 unsigned int wireno=std::round(wire_no);
 geo::WireID wireid(0,TPCb[clt],2,wireno);
 fWireReadoutGeom.WireEndPoints(wireid, xyzStart, xyzEnd);
 std::cout<<"Z position of intersection = "<<xyzStart[2]<<" "<<xyzEnd[2]<<"  "<<wireno<<std::endl;
 Zintersection.push_back(xyzStart[2]);
 timeintersection.push_back(ticks_no);
	    }

	}
      }



      ////*****************section for dE/dx correction*******************************////

	
    }//this track
     
    else if(thisShower != 0x0){
      fisprimarytrack               = 0;
      fisprimaryshower              = 1;

      fprimaryID                    = thisShower->ID();
      fprimaryLength                = thisShower->Length();
      fprimaryShowerBestPlane       = thisShower->best_plane();
      fprimaryOpeningAngle          = thisShower->OpenAngle();
      fprimaryStartPosition[0]      = thisShower->ShowerStart().X();
      fprimaryStartPosition[1]      = thisShower->ShowerStart().Y();
      fprimaryStartPosition[2]      = thisShower->ShowerStart().Z();
      fprimaryStartDirection[0]     = thisShower->Direction().X();
      fprimaryStartDirection[1]     = thisShower->Direction().Y();
      fprimaryStartDirection[2]     = thisShower->Direction().Z();
      if( (thisShower->Energy()).size() > 0 )
	fprimaryShowerEnergy = thisShower->Energy()[0]; // thisShower->best_plane()
      if( (thisShower->MIPEnergy()).size() > 0 )
	fprimaryShowerMIPEnergy = thisShower->MIPEnergy()[0];
      if( (thisShower->dEdx()).size() > 0 )
	fprimaryShowerdEdx = thisShower->dEdx()[0];
    }
    else{
      std::cout << "INFO::Primary pfParticle is not track or shower. Skip!" << std::endl;
      continue;
    }
    
    // Find the particle vertex. We need the tracker tag here because we need to do a bit of
    // additional work if the PFParticle is track-like to find the vertex. 
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fvertex[0] = vtx.X(); fvertex[1] = vtx.Y(); fvertex[2] = vtx.Z();

    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fsecvertex[0] = interactionVtx.X(); fsecvertex[1] = interactionVtx.Y(); fsecvertex[2] = interactionVtx.Z();

    // Maximum number of daugthers to be processed
    if(particle->NumDaughters() > NMAXDAUGTHERS)
      std::cout << "INFO::Number of daughters is " << particle->NumDaughters() << ". Only the first NMAXDAUGTHERS are processed." << std::endl;

    // Let's get the daughter PFParticles... we can do this simply without the utility
    for(const int daughterID : particle->Daughters()){
      // Daughter ID is the element of the original recoParticle vector
      const recob::PFParticle *daughterParticle      = &(recoParticles->at(daughterID));
      std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;
      
      const recob::Track* daughterTrack              = pfpUtil.GetPFParticleTrack(*daughterParticle,evt,fPFParticleTag,fTrackerTag);
      const recob::Shower* daughterShower            = pfpUtil.GetPFParticleShower(*daughterParticle,evt,fPFParticleTag,fShowerTag);
  
      if(daughterTrack != 0x0){
	fisdaughtertrack[fNDAUGHTERS]                = 1;
	fisdaughtershower[fNDAUGHTERS]               = 0;
	fdaughterTheta[fNDAUGHTERS]                  = daughterTrack->Theta();
	fdaughterPhi[fNDAUGHTERS]                    = daughterTrack->Phi();
	fdaughterLength[fNDAUGHTERS]                 = daughterTrack->Length();
	fdaughterMomentum[fNDAUGHTERS]               = daughterTrack->StartMomentum();
	fdaughterEndMomentum[fNDAUGHTERS]            = daughterTrack->EndMomentum();
	fdaughterStartPosition[fNDAUGHTERS][0]       = daughterTrack->Trajectory().Start().X();
	fdaughterStartPosition[fNDAUGHTERS][1]       = daughterTrack->Trajectory().Start().Y();
	fdaughterStartPosition[fNDAUGHTERS][2]       = daughterTrack->Trajectory().Start().Z();
	fdaughterEndPosition[fNDAUGHTERS][0]         = daughterTrack->Trajectory().End().X();
	fdaughterEndPosition[fNDAUGHTERS][1]         = daughterTrack->Trajectory().End().Y();
	fdaughterEndPosition[fNDAUGHTERS][2]         = daughterTrack->Trajectory().End().Z();
	fdaughterStartDirection[fNDAUGHTERS][0]      = daughterTrack->Trajectory().StartDirection().X();
	fdaughterStartDirection[fNDAUGHTERS][1]      = daughterTrack->Trajectory().StartDirection().Y();
	fdaughterStartDirection[fNDAUGHTERS][2]      = daughterTrack->Trajectory().StartDirection().Z();
	fdaughterEndDirection[fNDAUGHTERS][0]        = daughterTrack->Trajectory().EndDirection().X();
	fdaughterEndDirection[fNDAUGHTERS][1]        = daughterTrack->Trajectory().EndDirection().Y();
	fdaughterEndDirection[fNDAUGHTERS][2]        = daughterTrack->Trajectory().EndDirection().Z();

	fdaughterMomentumByRangeMuon[fNDAUGHTERS]    = trmom.GetTrackMomentum(daughterTrack->Length(),13);
	fdaughterMomentumByRangeProton[fNDAUGHTERS]  = trmom.GetTrackMomentum(daughterTrack->Length(),2212);

	std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
	if(daughtercalovector.size() != 3)
	  std::cerr << "WARNING::Calorimetry vector size for daughter is = " << daughtercalovector.size() << std::endl;

	for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
	  fdaughterKineticEnergy[fNDAUGHTERS][k] = daughtercalovector[k].KineticEnergy();
	  fdaughterRange[fNDAUGHTERS][k] = daughtercalovector[k].Range();
	}

	// Get the true mc particle
        const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *daughterTrack, evt, fTrackerTag);
	if(mcdaughterparticle != 0x0){
	  fdaughter_truth_TrackId[fNDAUGHTERS]          = mcdaughterparticle->TrackId();
	  fdaughter_truth_Pdg[fNDAUGHTERS]              = mcdaughterparticle->PdgCode();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][0] = mcdaughterparticle->Vx();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][1] = mcdaughterparticle->Vy();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][2] = mcdaughterparticle->Vz();
	  fdaughter_truth_StartPosition[fNDAUGHTERS][3] = mcdaughterparticle->T();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][0]   = mcdaughterparticle->EndX();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][1]   = mcdaughterparticle->EndY();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][2]   = mcdaughterparticle->EndZ();
	  fdaughter_truth_EndPosition[fNDAUGHTERS][3]   = mcdaughterparticle->EndT();
	  fdaughter_truth_P[fNDAUGHTERS]                = mcdaughterparticle->P();
	  fdaughter_truth_Momentum[fNDAUGHTERS][0]      = mcdaughterparticle->Px();
	  fdaughter_truth_Momentum[fNDAUGHTERS][1]      = mcdaughterparticle->Py();
	  fdaughter_truth_Momentum[fNDAUGHTERS][2]      = mcdaughterparticle->Pz();
	  fdaughter_truth_Momentum[fNDAUGHTERS][3]      = mcdaughterparticle->E();
	  fdaughter_truth_Pt[fNDAUGHTERS]               = mcdaughterparticle->Pt();
	  fdaughter_truth_Mass[fNDAUGHTERS]             = mcdaughterparticle->Mass();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][0]   = mcdaughterparticle->EndPx();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][1]   = mcdaughterparticle->EndPy();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][2]   = mcdaughterparticle->EndPz();
	  fdaughter_truth_EndMomentum[fNDAUGHTERS][3]   = mcdaughterparticle->EndE();
	  fdaughter_truth_Theta[fNDAUGHTERS]            = mcdaughterparticle->Momentum().Theta();
	  fdaughter_truth_Phi[fNDAUGHTERS]              = mcdaughterparticle->Momentum().Phi();
	  fdaughter_truth_Process[fNDAUGHTERS]          = int(mcdaughterparticle->Trajectory().ProcessToKey(mcdaughterparticle->Process()));
	  std::cout << "Daughter Process = " << (mcdaughterparticle->Process()).c_str() 
		    << " , mother = " << mcdaughterparticle->Mother() 
		    << std::endl;
	}
      }
      else if(daughterShower != 0x0){
	fisdaughtertrack[fNDAUGHTERS]                = 0;
	fisdaughtershower[fNDAUGHTERS]               = 1;
	fdaughterLength[fNDAUGHTERS]                 = daughterShower->Length();
	fdaughterShowerBestPlane[fNDAUGHTERS]        = daughterShower->best_plane();
	fdaughterOpeningAngle[fNDAUGHTERS]           = daughterShower->OpenAngle();
	fdaughterStartPosition[fNDAUGHTERS][0]       = daughterShower->ShowerStart().X();
	fdaughterStartPosition[fNDAUGHTERS][1]       = daughterShower->ShowerStart().Y();
	fdaughterStartPosition[fNDAUGHTERS][2]       = daughterShower->ShowerStart().Z();
	fdaughterStartDirection[fNDAUGHTERS][0]      = daughterShower->Direction().X();
	fdaughterStartDirection[fNDAUGHTERS][1]      = daughterShower->Direction().Y();
	fdaughterStartDirection[fNDAUGHTERS][2]      = daughterShower->Direction().Z();
	if( (daughterShower->Energy()).size() > 0 )
	  fdaughterShowerEnergy[fNDAUGHTERS] = daughterShower->Energy()[0]; // thisShower->best_plane()
	if( (daughterShower->MIPEnergy()).size() > 0 )
	  fdaughterShowerMIPEnergy[fNDAUGHTERS] = daughterShower->MIPEnergy()[0];
	if( (daughterShower->dEdx()).size() > 0 )
	  fdaughterShowerdEdx[fNDAUGHTERS] = daughterShower->dEdx()[0];
      }
      else{
	std::cout << "INFO::Daughter pfParticle is not track or shower. Skip!" << std::endl;
	continue;
      }

      fdaughterID[fNDAUGHTERS]                       = daughterID;
      // NHits associated with this pfParticle
      fdaughterNHits[fNDAUGHTERS]                    = (pfpUtil.GetPFParticleHits(*daughterParticle,evt,fPFParticleTag)).size();
      // T0
      // std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*daughterParticle,evt,fPFParticleTag);
      // if(!pfT0vec.empty())
      //	fdaughterT0[fNDAUGHTERS] = pfdaughterT0vec[0].Time();

      fNDAUGHTERS++;

      // Only process NMAXDAUGTHERS
      if(fNDAUGHTERS > NMAXDAUGTHERS) break;

    }

    break;
  }//particle loop pfparticle

  // Fill trees
  if(beamTriggerEvent)
    fPandoraBeam->Fill();

  //fPandoraCosmics->Fill();
}//analyzer

void protoana::truepion::endJob(){

}

void protoana::truepion::FillCosmicsTree(art::Event const & evt, std::string pfParticleTag){

  // To fill

}



void protoana::truepion::Initialise(){
  
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;

  for(int k=0; k < 3; k++){
    fvertex[k] = -999.0;
    fsecvertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
  }
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -9999;
  fbeamtrackTime = -999.0;
  fbeamtrackID = -999;
  for(int l=0; l < 3; l++){
    fbeamtrackP[l] = -999.0;
    fbeamtrackPos[l] = -999.0;
    fbeamtrackDir[l] = -999.0;
  }

  //NumberBeamTrajectoryPoints=0; 
  beamtrk_x.clear();
  beamtrk_y.clear();
  beamtrk_z.clear();
  beamtrk_Px.clear();
  beamtrk_Py.clear();
  beamtrk_Pz.clear();
  beamtrk_Eng.clear();


 
  fisprimarytrack = 0;
  fisprimaryshower = 0;
  fNDAUGHTERS = 0;

  fprimaryBDTScore = -999.0;
  fprimaryNHits = -999;
  fprimaryTheta = -999.0;
  fprimaryPhi = -999.0;
  fprimaryLength = -999.0;
  fprimaryMomentum = -999.0;
  fprimaryEndMomentum = -999.0;
  fprimaryOpeningAngle = -999.0;
  fprimaryShowerBestPlane = -999;
  fprimaryShowerEnergy = -999;
  fprimaryShowerMIPEnergy = -999;
  fprimaryShowerdEdx = -999;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryMomentumByRangeMuon = -999.0;

  fprimary_truth_TrackId = -999;
  fprimary_truth_Pdg = -999;
  ftruthpdg=-999;
  fprimary_truth_P = -999.0;
  fprimary_truth_Pt = -999.0;
  fprimary_truth_Mass = -999.0;
  fprimary_truth_Theta = -999.0;
  fprimary_truth_Phi = -999.0;
  fprimary_truth_Process = -999;
  fprimary_truth_Isbeammatched = -999;
  fprimary_truth_NDaughters = -999;
  fprimary_truth_tracklength=-9999;
  fprimary_truth_EndProcess="";
  truth_last_process="";
  for(int l=0; l < 4; l++){
    fprimary_truth_StartPosition[l] = -999.0;
    fprimary_truth_EndPosition[l] = -999.0;
    fprimary_truth_Momentum[l] = -999.0;
    fprimary_truth_EndMomentum[l] = -999.0;
  }

  for(int k=0; k < NMAXDAUGTHERS; k++){
    fisdaughtertrack[k] = -999;
    fisdaughtershower[k] = -999;
    fdaughterNHits[k] = -999;
    fdaughterTheta[k] = -999.0;
    fdaughterPhi[k] = -999.0;
    fdaughterLength[k] = -999.0;
    fdaughterMomentum[k] = -999.0;
    fdaughterEndMomentum[k] = -999.0;
    for(int l=0; l < 3; l++){
      fdaughterEndPosition[k][l] = -999.0;
      fdaughterStartPosition[k][l] = -999.0;
      fdaughterEndDirection[k][l] = -999.0;
      fdaughterStartDirection[k][l] = -999.0;
      fdaughterKineticEnergy[k][l] = -999.0;
      fdaughterRange[k][l] = -999.0;
    }
    fdaughterOpeningAngle[k] = -999.0;
    fdaughterShowerBestPlane[k] = -999;
    fdaughterShowerEnergy[k] = -999;
    fdaughterShowerMIPEnergy[k] = -999;
    fdaughterShowerdEdx[k] = -999;
    fdaughterMomentumByRangeProton[k] = -999.0;
    fdaughterMomentumByRangeMuon[k] = -999.0;
    fdaughterID[k] = -999;
    // fdaughterT0[k] = -999;

    fdaughter_truth_TrackId[k] = -999;
    fdaughter_truth_Pdg[k] = -999;
    fdaughter_truth_P[k] = -999.0;
    fdaughter_truth_Pt[k] = -999.0;
    fdaughter_truth_Mass[k] = -999.0;
    fdaughter_truth_Theta[k] = -999.0;
    fdaughter_truth_Phi[k] = -999.0;
    fdaughter_truth_Process[k] = -999;
    for(int l=0; l < 4; l++){
      fdaughter_truth_StartPosition[k][l] = -999.0;
      fdaughter_truth_EndPosition[k][l] = -999.0;
      fdaughter_truth_Momentum[k][l] = -999.0;
      fdaughter_truth_EndMomentum[k][l] = -999.0;
    }
  }

  primtrk_dqdx.clear();
  primtrk_resrange.clear();
  primtrk_dedx.clear();
  primtrk_range.clear();
  primtrk_hitx.clear();
  primtrk_hity.clear();
  primtrk_hitz.clear();
  primtrk_pitch.clear();
  primtrk_truth_Z.clear();
  primtrk_truth_Eng.clear();
  primtrk_truth_trkide.clear();
  interactionX.clear();
  interactionY.clear();
  interactionZ.clear();
  interactionProcesslist.clear();
  interactionAngles.clear();
  Zintersection.clear();
  timeintersection.clear();
}

DEFINE_ART_MODULE(protoana::truepion)
