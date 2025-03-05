////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEelectronAnaTree
// File:        ProtoDUNEelectronAnaTree_module.cc
//
// Extract protoDUNE useful information, do a firs tpre-selection and save output to a flat tree
// 
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
// modified by Aaron Higuera ahiguera@central.uh.edu
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Wire.h"
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
#include "larcore/Geometry/WireReadout.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dunecalib/CalibServices/XYZCalibServiceProtoDUNE.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"

// C++ Includes
#include <stdio.h>
#include <stdlib.h> 
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// Maximum number of beam particles to save
const int NMAXDAUGTHERS = 25;
const int MAXHits = 6000;
namespace protoana {
  class ProtoDUNEelectronAnaTree;
}


class protoana::ProtoDUNEelectronAnaTree : public art::EDAnalyzer {
public:

  explicit ProtoDUNEelectronAnaTree(fhicl::ParameterSet const & p);

  ProtoDUNEelectronAnaTree(ProtoDUNEelectronAnaTree const &) = delete;
  ProtoDUNEelectronAnaTree(ProtoDUNEelectronAnaTree &&) = delete;
  ProtoDUNEelectronAnaTree & operator = (ProtoDUNEelectronAnaTree const &) = delete;
  ProtoDUNEelectronAnaTree & operator = (ProtoDUNEelectronAnaTree &&) = delete;

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
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNEBeamlineUtils beamlineUtil;


  // Track momentum algorithm calculates momentum based on track range
  trkf::TrackMomentumCalculator trmom;
  geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

  // Initialize tree variables
  void Initialize();
  void FillPrimaryPFParticle(art::Event const & evt,
                             detinfo::DetectorClocksData const& clockData,
                             detinfo::DetectorPropertiesData const& detProp,
                             const recob::PFParticle* particle);
  void FillPrimaryDaughterPFParticle(art::Event const & evt, const recob::PFParticle* daughterParticle, int daughterID);

  // Fill cosmics tree

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fCalorimetryTag;
  std::string fParticleIDTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fHitTag;
  std::string fShowerCaloTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  int fVerbose;

  TTree *fPandoraBeam;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];

  // Beam track
  int fbeamtrigger;
  int fbeamCheckIsMatched;
  double ftof;
  int fcerenkovStatus[2];
  double fcerenkovTime[2];
  double fcerenkovPressure[2];
  double fbeamtrackMomentum;
  double fbeamtrackEnergy;
  double fbeamtrackPos[3];
  double fbeamtrackEndPos[3];
  double fbeamtrackDir[3];
  double fbeamtrackTime;
  int fbeam_ntrjPoints;
  int fbeamtrackPdg;
  int fbeamtrackID;
  double fbeamtrackPos_at[1000][3];
  double fbeamtrackMom_at[1000][4];
  std::vector<std::string> fbeamtrackEndProcess;
  int fbeamtrackNDaughters;

  // Reconstructed tracks/showers
  double fprimaryVertex[3];
  int fprimaryIstrack;
  int fprimaryIsshower;
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
  double fprimaryShowerCharge;
  double fprimaryShowerMIPEnergy;
  double fprimaryMomentumByRangeProton;
  double fprimaryMomentumByRangeMuon;
  int fprimaryIsBeamparticle;
  int    fprimaryTruth_trkID;
  int    fprimaryTruth_pdg;
  double fprimaryTruth_Edepo;
  double fprimaryTruth_purity;
  double fprimaryTruth_E;
  double fprimaryTruth_vtx[3];
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  int    fprimarynCal;
  double fprimarydEdx[MAXHits];
  double fprimarydQdx[MAXHits];
  double fprimary_calX[MAXHits];
  double fprimary_calY[MAXHits];
  double fprimary_calZ[MAXHits];
  double fprimary_cal_pitch[MAXHits];
  double fprimaryResidualRange[MAXHits];
  int    fprimaryTruthShower_nHits;
  double fprimaryShowerTruth_Charge;
  int    fprimaryShower_nHits; //collection only
  int    fprimaryShower_hit_w[MAXHits];
  double fprimaryShower_hit_q[MAXHits];
  double fprimaryShower_hit_t[MAXHits]; 
  double fprimaryShower_hit_X[MAXHits];
  int    fprimaryShower_hit_ch[MAXHits];
  double fprimaryShower_hit_Y[MAXHits];
  double fprimaryShower_hit_Z[MAXHits];
  double fprimaryShower_hit_pitch[MAXHits]; 
  double fprimaryShower_hit_cnn[MAXHits];
  int    fprimaryTruthShower_hit_w[MAXHits];
  double fprimaryTruthShower_hit_q[MAXHits];
  double fprimaryTruthShower_hit_t[MAXHits]; 
  double fprimaryTruthShower_hit_X[MAXHits];
  double fprimaryTruthShower_hit_Y[MAXHits];
  double fprimaryTruthShower_hit_Z[MAXHits];


  int fprimaryID;
  double fprimaryT0;

  int fNDAUGHTERS =0;
  double fdaughterVertex[3];
  int fdaughterIstrack[NMAXDAUGTHERS];
  int fdaughterIsshower[NMAXDAUGTHERS];
  int fdaughterNHits[NMAXDAUGTHERS];
  double fdaughterTheta[NMAXDAUGTHERS];
  double fdaughterPhi[NMAXDAUGTHERS];
  double fdaughterLength[NMAXDAUGTHERS];
  double fdaughterEndPosition[NMAXDAUGTHERS][3];
  double fdaughterStartPosition[NMAXDAUGTHERS][3];
  double fdaughterEndDirection[NMAXDAUGTHERS][3];
  double fdaughterStartDirection[NMAXDAUGTHERS][3];
  double fdaughterOpeningAngle[NMAXDAUGTHERS];
  int fdaughterShowerBestPlane[NMAXDAUGTHERS];
  int fdaughterID[NMAXDAUGTHERS];
  double fdaughterT0[NMAXDAUGTHERS];

};


protoana::ProtoDUNEelectronAnaTree::ProtoDUNEelectronAnaTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  beamlineUtil(p.get<fhicl::ParameterSet>("BeamLineUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fShowerCaloTag(p.get<std::string>("ShowerCalorimetryTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),

  fVerbose(p.get<int>("Verbose"))
{

}

void protoana::ProtoDUNEelectronAnaTree::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  fPandoraBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
  fPandoraBeam->Branch("run",                           &fRun,                          "run/I");
  fPandoraBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
  fPandoraBeam->Branch("event",                         &fevent,                        "event/I");
  fPandoraBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
  fPandoraBeam->Branch("Nactivefembs",                  &fNactivefembs,                 "Nactivefembs[5]/I");

  fPandoraBeam->Branch("beamtrigger",                   &fbeamtrigger,                  "beamtrigger/I");
  fPandoraBeam->Branch("beamCheckIsMatched",            &fbeamCheckIsMatched,           "beamCheckIsMatched/I");
  fPandoraBeam->Branch("tof",                           &ftof,                          "tof/D");
  fPandoraBeam->Branch("cerenkovStatus",                &fcerenkovStatus,               "cerenkovStatus[2]/I");
  fPandoraBeam->Branch("cerenkovTime",                  &fcerenkovTime,                 "cerenkovTime[2]/D");
  fPandoraBeam->Branch("cerenkovPressure",              &fcerenkovPressure,             "cerenkovPressure[2]/D");
  fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
  fPandoraBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
  fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
  fPandoraBeam->Branch("beamtrackEndPos",               &fbeamtrackEndPos,              "beamtrackEndPos[3]/D");
  fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
  fPandoraBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
  fPandoraBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
  fPandoraBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");
  fPandoraBeam->Branch("beam_ntrjPoints",               &fbeam_ntrjPoints,              "beam_ntrjPoints/I");
  fPandoraBeam->Branch("beamtrackNDaughters",           &fbeamtrackNDaughters,          "beamtrackNDaughters/I");
  fPandoraBeam->Branch("beamtrackPos_at",               &fbeamtrackPos_at,              "beamtrackPos_at[beam_ntrjPoints][3]/D");
  fPandoraBeam->Branch("beamtrackMom_at",               &fbeamtrackMom_at,              "beamtrackMom_at[beam_ntrjPoints][4]/D");
  fPandoraBeam->Branch("beamtrackEndProcess",           &fbeamtrackEndProcess);

  fPandoraBeam->Branch("primaryVertex",                 &fprimaryVertex,                "primaryVertex[3]/D");
  fPandoraBeam->Branch("primaryIsBeamparticle",         &fprimaryIsBeamparticle,        "primaryIsBeamparticle/I");
  fPandoraBeam->Branch("primaryIstrack",                &fprimaryIstrack,               "primaryIstrack/I");
  fPandoraBeam->Branch("primaryIsshower",               &fprimaryIsshower,              "primaryIsshower/I");
  fPandoraBeam->Branch("primaryBDTScore",               &fprimaryBDTScore,              "primaryBDTScore/D");
  fPandoraBeam->Branch("primaryNHits",                  &fprimaryNHits,                 "primaryNHits/I");
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
  fPandoraBeam->Branch("primaryTruth_Edepo",            &fprimaryTruth_Edepo,           "primaryTruth_Edepo/D");
  fPandoraBeam->Branch("primaryTruth_purity",           &fprimaryTruth_purity,          "primaryTruth_purity/D");
  fPandoraBeam->Branch("primaryTruth_E",                &fprimaryTruth_E,               "primaryTruth_E/D");
  fPandoraBeam->Branch("primaryTruth_vtx",              &fprimaryTruth_vtx,             "primaryTruth_vtx[3]/D");
  fPandoraBeam->Branch("primaryTruth_pdg",              &fprimaryTruth_pdg,             "primaryTruth_pdg/I");
  fPandoraBeam->Branch("primaryTruth_trkID",            &fprimaryTruth_trkID,           "primaryTruth_trkID/I");
  fPandoraBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
  fPandoraBeam->Branch("primaryShowerCharge",           &fprimaryShowerCharge,          "primaryShowerCharge/D");
  fPandoraBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/D");
  fPandoraBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/D");

  fPandoraBeam->Branch("primaryShower_nHits",        &fprimaryShower_nHits,       "primaryShower_nHits/I");
  fPandoraBeam->Branch("primaryTruthShower_nHits",        &fprimaryTruthShower_nHits,       "primaryTruthShower_nHits/I");
  fPandoraBeam->Branch("primaryShowerTruth_Charge",           &fprimaryShowerTruth_Charge,          "primaryShowerTruth_Charge/D");
  fPandoraBeam->Branch("primaryShower_hit_q",        &fprimaryShower_hit_q,       "primaryShower_hit_q[primaryShower_nHits]/D");
  fPandoraBeam->Branch("primaryShower_hit_w",        &fprimaryShower_hit_w,       "primaryShower_hit_w[primaryShower_nHits]/I");
  fPandoraBeam->Branch("primaryShower_hit_t",        &fprimaryShower_hit_t,       "primaryShower_hit_t[primaryShower_nHits]/D");
  fPandoraBeam->Branch("primaryShower_hit_X",        &fprimaryShower_hit_X,       "primaryShower_hit_X[primaryShower_nHits]/D");
  fPandoraBeam->Branch("primaryShower_hit_Y",        &fprimaryShower_hit_Y,       "primaryShower_hit_Y[primaryShower_nHits]/D");
  fPandoraBeam->Branch("primaryShower_hit_Z",        &fprimaryShower_hit_Z,       "primaryShower_hit_Z[primaryShower_nHits]/D");
  fPandoraBeam->Branch("primaryShower_hit_ch",        &fprimaryShower_hit_ch,       "primaryShower_hit_ch[primaryShower_nHits]/I");
 
  fPandoraBeam->Branch("primaryTruthShower_hit_q",        &fprimaryTruthShower_hit_q,       "primaryTruthShower_hit_q[primaryTruthShower_nHits]/D");
  fPandoraBeam->Branch("primaryTruthShower_hit_w",        &fprimaryTruthShower_hit_w,       "primaryTruthShower_hit_w[primaryTruthShower_nHits]/I");
  fPandoraBeam->Branch("primaryTruthShower_hit_t",        &fprimaryTruthShower_hit_t,       "primaryTruthShower_hit_t[primaryTruthShower_nHits]/D");
  fPandoraBeam->Branch("primaryTruthShower_hit_X",        &fprimaryTruthShower_hit_X,       "primaryTruthShower_hit_X[primaryTruthShower_nHits]/D");
  fPandoraBeam->Branch("primaryTruthShower_hit_Y",        &fprimaryTruthShower_hit_Y,       "primaryTruthShower_hit_Y[primaryTruthShower_nHits]/D");
  fPandoraBeam->Branch("primaryTruthShower_hit_Z",        &fprimaryTruthShower_hit_Z,       "primaryTruthShower_hit_Z[primaryTruthShower_nHits]/D");
 
  fPandoraBeam->Branch("primaryShower_hit_pitch",    &fprimaryShower_hit_pitch,   "primaryShower_hit_pitch[primaryShower_nHits]/D");
  fPandoraBeam->Branch("primaryShower_hit_cnn",    &fprimaryShower_hit_cnn,   "primaryShower_hit_cnn[primaryShower_nHits]/D");

  fPandoraBeam->Branch("primaryMomentumByRangeProton",  &fprimaryMomentumByRangeProton, "primaryMomentumByRangeProton/D");
  fPandoraBeam->Branch("primaryMomentumByRangeMuon",    &fprimaryMomentumByRangeMuon,   "primaryMomentumByRangeMuon/D");
  fPandoraBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
  fPandoraBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
  fPandoraBeam->Branch("primarynCal",                   &fprimarynCal,                  "primarynCal/I");
  fPandoraBeam->Branch("primarydQdx",              	&fprimarydQdx,	 	        "primarydQdx[primarynCal]/D");
  fPandoraBeam->Branch("primary_calX",              	&fprimary_calX,	 	        "primary_calX[primarynCal]/D");
  fPandoraBeam->Branch("primary_calY",              	&fprimary_calY,	 	        "primary_calY[primarynCal]/D");
  fPandoraBeam->Branch("primary_calZ",              	&fprimary_calZ,	 	        "primary_calZ[primarynCal]/D");
  fPandoraBeam->Branch("primary_cal_pitch",             &fprimary_cal_pitch,	 	"primary_cal_pitch[primarynCal]/D");
  fPandoraBeam->Branch("primarydEdx",              	&fprimarydEdx,	 	        "primarydEdx[primarynCal]/D");
  fPandoraBeam->Branch("primaryResidualRange",          &fprimaryResidualRange,	 	"primaryResidualRange[primarynCal]/D");
  fPandoraBeam->Branch("primaryT0",                     &fprimaryT0,                    "primaryT0/D");

  fPandoraBeam->Branch("NDAUGHTERS",                    &fNDAUGHTERS,                    "NDAUGHTERS/I");
  fPandoraBeam->Branch("daughterVertex",                &fdaughterVertex,                "daughterVertex[3]/D");
  fPandoraBeam->Branch("daughterIstrack",               &fdaughterIstrack,               "daughterIstrack[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterIsshower",              &fdaughterIsshower,              "daughterIsshower[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterNHits",                 &fdaughterNHits,                 "daughterNHits[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterTheta",                 &fdaughterTheta,                 "daughterTheta[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterPhi",                   &fdaughterPhi,                   "daughterPhi[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterLength",                &fdaughterLength,                "daughterLength[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterEndPosition",           &fdaughterEndPosition,           "daughterEndPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartPosition",         &fdaughterStartPosition,         "daughterStartPosition[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterStartDirection",        &fdaughterStartDirection,        "daughterStartDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterEndDirection",          &fdaughterEndDirection,          "daughterEndDirection[NDAUGHTERS][3]/D");
  fPandoraBeam->Branch("daughterOpeningAngle",          &fdaughterOpeningAngle,          "daughterOpeningAngle[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterShowerBestPlane",       &fdaughterShowerBestPlane,       "daughterShowerBestPlane[NDAUGHTERS]/D");
  fPandoraBeam->Branch("daughterID",                    &fdaughterID,                    "daughterID[NDAUGHTERS]/I");
  fPandoraBeam->Branch("daughterT0",                    &fdaughterT0,                    "daughterT0[NDAUGHTERS]/D");

}

void protoana::ProtoDUNEelectronAnaTree::analyze(art::Event const & evt){

  // Initialize tree parameters
  Initialize();
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

  // Get number of active fembs
  if(!evt.isRealData()){
    for(int k=0; k < 6; k++)
      fNactivefembs[k] = 20;
  }
  else{
    for(int k=0; k < 6; k++)
     fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
  }

  bool beamTriggerEvent = false;
  if(!evt.isRealData()){
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);

    // Also get the reconstructed beam information in the MC - TO DO
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    if(geantGoodParticle != 0x0){
      beamTriggerEvent = true;
      fbeamtrigger       = 12;
      fbeamtrackPos[0]   = geantGoodParticle->Vx();
      fbeamtrackPos[1]   = geantGoodParticle->Vy();
      fbeamtrackPos[2]   = geantGoodParticle->Vz();
      fbeamtrackEndPos[0]= geantGoodParticle->EndX();
      fbeamtrackEndPos[1]= geantGoodParticle->EndY();
      fbeamtrackEndPos[2]= geantGoodParticle->EndZ();
      fbeamtrackMomentum = geantGoodParticle->P();
      fbeamtrackEnergy   = geantGoodParticle->E();
      fbeamtrackPdg      = geantGoodParticle->PdgCode();
      fbeamtrackTime     = geantGoodParticle->T();
      fbeamtrackID       = geantGoodParticle->TrackId();
      fbeamtrackEndProcess.push_back(geantGoodParticle->EndProcess()); 
      fbeam_ntrjPoints = geantGoodParticle->NumberTrajectoryPoints();
      fbeamtrackNDaughters = geantGoodParticle->NumberDaughters(); 
      for( size_t i=0; i<geantGoodParticle->NumberTrajectoryPoints() && i<1000; ++i){
         fbeamtrackPos_at[i][0] = geantGoodParticle->Position(i).X();
         fbeamtrackPos_at[i][1] = geantGoodParticle->Position(i).Y();
         fbeamtrackPos_at[i][2] = geantGoodParticle->Position(i).Z();
         fbeamtrackMom_at[i][0] = geantGoodParticle->Momentum(i).Px();
         fbeamtrackMom_at[i][1] = geantGoodParticle->Momentum(i).Py();
         fbeamtrackMom_at[i][2] = geantGoodParticle->Momentum(i).Pz();
         fbeamtrackMom_at[i][3] = geantGoodParticle->Momentum(i).E();
      } 
    }

    if( abs(geantGoodParticle->PdgCode()) != 11 ) return;
  } // MC
  else{ //data
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);
    if( !beamTriggerEvent ) return;

    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    auto pdbeamHandle = evt.getHandle< std::vector<beam::ProtoDUNEBeamEvent> >(fBeamModuleLabel);
    if (pdbeamHandle)
      art::fill_ptr_vector(beaminfo, pdbeamHandle);
    else{
       std::cout<<"No beam information from "<<fBeamModuleLabel<<std::endl;
    } 
    if(beaminfo.size()){
      if( beamlineUtil.IsGoodBeamlineTrigger( evt ) ){  
        fbeamCheckIsMatched = beaminfo[0]->CheckIsMatched();
        fbeamtrigger = beaminfo[0]->GetTimingTrigger();
        if(beaminfo[0]->GetTOFChan() != -1){//if TOFChan == -1, then there was not a successful match, if it's 0, 1, 2, or 3, then there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
          ftof = beaminfo[0]->GetTOF();
        }
        // Get Cerenkov
        fcerenkovStatus[0]   = beaminfo[0]->GetCKov0Status();
        fcerenkovStatus[1]   = beaminfo[0]->GetCKov1Status();
        fcerenkovTime[0]     = beaminfo[0]->GetCKov0Time();
        fcerenkovTime[1]     = beaminfo[0]->GetCKov1Time();
        fcerenkovPressure[0] = beaminfo[0]->GetCKov0Pressure();
        fcerenkovPressure[1] = beaminfo[0]->GetCKov1Pressure();
        auto & tracks = beaminfo[0]->GetBeamTracks();
        if(!tracks.empty()){
          fbeamtrackPos[0] = tracks[0].Start().X();
          fbeamtrackPos[1] = tracks[0].Start().Y();
          fbeamtrackPos[2] = tracks[0].Start().Z();
          fbeamtrackEndPos[0] = tracks[0].End().X();
          fbeamtrackEndPos[1] = tracks[0].End().Y();
          fbeamtrackEndPos[2] = tracks[0].End().Z();
          fbeamtrackDir[0] = tracks[0].StartDirection().X();
          fbeamtrackDir[1] = tracks[0].StartDirection().Y();
          fbeamtrackDir[2] = tracks[0].StartDirection().Z();
        }
        auto & beammom = beaminfo[0]->GetRecoBeamMomenta();
        if(!beammom.empty()) fbeamtrackMomentum = beammom[0];
      } //good beam trigger
    }
  }//for data

  //check for reco pandora stuff
  auto recoParticleHandle = evt.getHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  if( !recoParticleHandle ) return;

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
  for(const recob::PFParticle* particle : pfParticles){

    FillPrimaryPFParticle(evt, clockData, detProp, particle);
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
     
    fprimaryVertex[0] = vtx.X(); fprimaryVertex[1] = vtx.Y(); fprimaryVertex[2] = vtx.Z();
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fdaughterVertex[0] = interactionVtx.X(); fdaughterVertex[1] = interactionVtx.Y(); fdaughterVertex[2] = interactionVtx.Z();
    
    fNDAUGHTERS =0;
    for( const int didx :  particle->Daughters() ){
      const recob::PFParticle *daughterParticle      =&(recoParticles->at(didx));
      FillPrimaryDaughterPFParticle(evt, daughterParticle, didx);
      // Only process NMAXDAUGTHERS
      if(fNDAUGHTERS > NMAXDAUGTHERS) break;
      fNDAUGHTERS ++;
    } 
    
    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  }
  if(beamTriggerEvent)   fPandoraBeam->Fill();
}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEelectronAnaTree::endJob(){

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEelectronAnaTree::FillPrimaryPFParticle(art::Event const & evt,
                                                               detinfo::DetectorClocksData const& clockData,
                                                               detinfo::DetectorPropertiesData const& detProp,
                                                               const recob::PFParticle* particle){

  // Pandora's BDT beam-cosmic score
  fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
  // NHits associated with this pfParticle
  fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();
  
  // Get the T0 for this pfParticle
  std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
  if(!pfT0vec.empty())
    fprimaryT0 = pfT0vec[0].Time();

  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
  const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle, evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);

  const simb::MCParticle* mcparticle = NULL; 
  if(thisTrack != 0x0){

    mcparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, evt, fTrackerTag);
    if(mcparticle){
      fprimaryTruth_pdg = mcparticle->PdgCode();        
      fprimaryTruth_trkID = mcparticle->TrackId();
      fprimaryTruth_vtx[0] =mcparticle->Vx();
      fprimaryTruth_vtx[1] =mcparticle->Vy();
      fprimaryTruth_vtx[2] =mcparticle->Vx();
      fprimaryTruth_E   = mcparticle->E();
      if( fbeamtrackID != -999 && fbeamtrackID == fprimaryTruth_trkID ) 
        fprimaryIsBeamparticle = 1;
    }
    fprimaryIstrack                    = 1;
    fprimaryIsshower                   = 0;
    fprimaryID                         = thisTrack->ParticleId();
    fprimaryTheta                      = thisTrack->Theta();
    fprimaryPhi                        = thisTrack->Phi();
    fprimaryLength                     = thisTrack->Length();
    fprimaryMomentum                   = thisTrack->StartMomentum();
    fprimaryEndMomentum                = thisTrack->EndMomentum();
    fprimaryEndPosition[0]             = thisTrack->Trajectory().End().X();
    fprimaryEndPosition[1]             = thisTrack->Trajectory().End().Y();
    fprimaryEndPosition[2]             = thisTrack->Trajectory().End().Z();
    fprimaryStartPosition[0]           = thisTrack->Trajectory().Start().X();
    fprimaryStartPosition[1]           = thisTrack->Trajectory().Start().Y();
    fprimaryStartPosition[2]           = thisTrack->Trajectory().Start().Z();
    fprimaryEndDirection[0]            = thisTrack->Trajectory().EndDirection().X();
    fprimaryEndDirection[1]            = thisTrack->Trajectory().EndDirection().Y();
    fprimaryEndDirection[2]            = thisTrack->Trajectory().EndDirection().Z();
    fprimaryStartDirection[0]          = thisTrack->Trajectory().StartDirection().X();
    fprimaryStartDirection[1]          = thisTrack->Trajectory().StartDirection().Y();
    fprimaryStartDirection[2]          = thisTrack->Trajectory().StartDirection().Z();
    
    fprimaryMomentumByRangeMuon        = trmom.GetTrackMomentum(thisTrack->Length(),13);
    fprimaryMomentumByRangeProton      = trmom.GetTrackMomentum(thisTrack->Length(),2212);      

    // Calorimetry only colleciton plane
    std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
    for(size_t k = 0; k < calovector.size() && k<3; k++){
      int plane = calovector[k].PlaneID().Plane;
      if(plane !=2 ) continue;
      fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
      fprimaryRange[plane]              = calovector[k].Range();
      fprimarynCal = calovector[k].dEdx().size();
      for(size_t l=0; l<calovector[k].dEdx().size() && l<MAXHits; ++l){
         fprimarydEdx[l]= calovector[k].dEdx()[l];
         fprimarydQdx[l]= calovector[k].dQdx()[l];
         fprimaryResidualRange[l]= calovector[k].ResidualRange()[l];
         const auto &pos=(calovector[k].XYZ())[l];
         fprimary_calX[l] = pos.X();
         fprimary_calY[l] = pos.Y();
         fprimary_calZ[l] = pos.Z();
         fprimary_cal_pitch[l] =calovector[k].TrkPitchVec()[l];
      }
    }
  } // end is track
  else if(thisShower != 0x0){

    mcparticle = truthUtil.GetMCParticleFromRecoShower(clockData, *thisShower, evt, fShowerTag);
    if(mcparticle){
      fprimaryTruth_pdg = mcparticle->PdgCode();        
      fprimaryTruth_trkID = mcparticle->TrackId();
      fprimaryTruth_vtx[0] =mcparticle->Vx();
      fprimaryTruth_vtx[1] =mcparticle->Vy();
      fprimaryTruth_vtx[2] =mcparticle->Vx();
      fprimaryTruth_E   = mcparticle->E();
      if( fbeamtrackID != -999 && fbeamtrackID == fprimaryTruth_trkID ) 
        fprimaryIsBeamparticle = 1;

      double Edepo = truthUtil.GetDepEnergyMC(evt,fWireReadoutGeom, mcparticle->TrackId(), 2 );
      double purity = truthUtil.GetPurity(clockData, *thisShower, evt, fShowerTag);
      fprimaryTruth_Edepo = Edepo;
      fprimaryTruth_purity = purity; //not really useful in this case
      
      double tot_ch =0;
      std::vector<const recob::Hit*> hitsFromMCPart = truthUtil.GetMCParticleHits(clockData, *mcparticle, evt, fHitTag);
      art::FindManyP<recob::SpacePoint> spFromMCPartHits(hitsFromMCPart,evt,fPFParticleTag);

      if( hitsFromMCPart.size()){
        int n_hits =0;
        for( size_t i=0; i<hitsFromMCPart.size(); ++i){
           if( hitsFromMCPart[i]->WireID().Plane != 2 ) continue;
           const geo::WireGeo* pwire = fWireReadoutGeom.WirePtr(hitsFromMCPart[i]->WireID());
           auto xyzWire = pwire->GetCenter();
           tot_ch += hitsFromMCPart[i]->Integral(); 
           fprimaryTruthShower_hit_w[n_hits]=hitsFromMCPart[i]->WireID().Wire;
           fprimaryTruthShower_hit_t[n_hits]=hitsFromMCPart[i]->PeakTime();
           fprimaryTruthShower_hit_q[n_hits]=hitsFromMCPart[i]->Integral(); 
           fprimaryTruthShower_hit_X[n_hits]=detProp.ConvertTicksToX(hitsFromMCPart[i]->PeakTime(),hitsFromMCPart[i]->WireID().Plane,hitsFromMCPart[i]->WireID().TPC,0);
           fprimaryTruthShower_hit_Z[n_hits] = xyzWire.Z();
           std::vector<art::Ptr<recob::SpacePoint>> sp = spFromMCPartHits.at(i); 
           if(!sp.empty() ){
             //fprimaryShower_hit_X[idx]= sp[0]->XYZ()[0];
             fprimaryTruthShower_hit_Y[n_hits]= sp[0]->XYZ()[1];
             //fprimaryTruthShower_hit_Z[n_hits]= sp[0]->XYZ()[2];
           }
           n_hits ++;
        }
        fprimaryShowerTruth_Charge = tot_ch; 
        fprimaryTruthShower_nHits = n_hits; //collection only
      }
    }
    fprimaryIstrack                     = 0;
    fprimaryIsshower                    = 1;
    fprimaryID                          = thisShower->ID();
    fprimaryLength                      = thisShower->Length();
    fprimaryShowerBestPlane             = thisShower->best_plane();
    fprimaryOpeningAngle                = thisShower->OpenAngle();
    fprimaryStartPosition[0]            = thisShower->ShowerStart().X();
    fprimaryStartPosition[1]            = thisShower->ShowerStart().Y();
    fprimaryStartPosition[2]            = thisShower->ShowerStart().Z();
    fprimaryStartDirection[0]           = thisShower->Direction().X();
    fprimaryStartDirection[1]           = thisShower->Direction().Y();
    fprimaryStartDirection[2]           = thisShower->Direction().Z();
   
    const std::vector<const recob::Hit*> sh_hits = showerUtil.GetRecoShowerHits(*thisShower, evt, fShowerTag);
    art::FindManyP<recob::SpacePoint> spFromShowerHits(sh_hits,evt,fPFParticleTag);

    //work around to save the CNN score 
    auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
    art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);
    anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );
    //int trk_idx = hitResults.getIndex("track");
    //int em_idx = hitResults.getIndex("em");
    std::vector< art::Ptr< recob::Hit > > tmp_sh_hits = findHitsFromShowers.at( thisShower->ID() );

    int idx =0;
    fprimaryShowerCharge =0.0;
    for( size_t j=0; j<sh_hits.size() && j<MAXHits; ++j){
       if( sh_hits[j]->WireID().Plane != 2 ) continue;
       art::FindManyP<recob::Wire> wFromHits(sh_hits,evt,"hitpdune");
       std::vector<art::Ptr<recob::Wire>> wires = wFromHits.at(j);
       const geo::WireGeo* pwire = fWireReadoutGeom.WirePtr(sh_hits[j]->WireID());
       auto xyzWire = pwire->GetCenter();
       std::array<float,4> cnn_out = hitResults.getOutput( tmp_sh_hits[j] );
       double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]; 
       double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh; 
       fprimaryShower_hit_cnn[idx] = cnn_score;
       fprimaryShowerCharge += sh_hits[j]->Integral();
       fprimaryShower_hit_w[idx]=sh_hits[j]->WireID().Wire;
       fprimaryShower_hit_t[idx]=sh_hits[j]->PeakTime();
       fprimaryShower_hit_q[idx]=sh_hits[j]->Integral(); 
       fprimaryShower_hit_X[idx]=detProp.ConvertTicksToX(sh_hits[j]->PeakTime(),sh_hits[j]->WireID().Plane,sh_hits[j]->WireID().TPC,0);
       fprimaryShower_hit_Z[idx]= xyzWire.Z();  
       fprimaryShower_hit_ch[idx]= wires[0]->Channel(); //only one channel per wire at collection plane
       std::vector<art::Ptr<recob::SpacePoint>> sp = spFromShowerHits.at(j); 

       if(!sp.empty() ){
         //fprimaryShower_hit_X[idx]= sp[0]->XYZ()[0];
         fprimaryShower_hit_Y[idx]= sp[0]->XYZ()[1];
         //fprimaryShower_hit_Z[idx]= sp[0]->XYZ()[2];
       }
    
       idx ++;
    } 
    fprimaryShower_nHits = idx; //only collection hits
    // Calorimetry only colleciton plane
    std::vector<anab::Calorimetry> calovector = showerUtil.GetRecoShowerCalorimetry(*thisShower, evt, fShowerTag, fShowerCaloTag);
    if(calovector.size() != 3 && fVerbose > 0)
      std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;
    
    for(size_t k = 0; k < calovector.size() && k<3; k++){
      int plane = calovector[k].PlaneID().Plane;
      if(plane !=2 ) continue;
      fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
      fprimaryRange[plane]              = calovector[k].Range();
      fprimarynCal = calovector[k].dEdx().size();
      for(size_t l=0; l<calovector[k].dEdx().size() && l<MAXHits; ++l){
         fprimarydEdx[l]= calovector[k].dEdx()[l];
         fprimarydQdx[l]= calovector[k].dQdx()[l];
         fprimaryResidualRange[l]= calovector[k].ResidualRange()[l];
         const auto &pos=(calovector[k].XYZ())[l];
         fprimary_calX[l] = pos.X();
         fprimary_calY[l] = pos.Y();
         fprimary_calZ[l] = pos.Z();
         fprimary_cal_pitch[l] =calovector[k].TrkPitchVec()[l];
      }
    }

  } // end is shower

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEelectronAnaTree::FillPrimaryDaughterPFParticle(art::Event const & evt, const recob::PFParticle* daughterParticle, int daughterID){
  art::Handle<std::vector<recob::Hit>> HitHandle;
 
  const recob::Track* daughterTrack              = pfpUtil.GetPFParticleTrack(*daughterParticle,evt, fPFParticleTag,fTrackerTag);
  const recob::Shower* daughterShower            = pfpUtil.GetPFParticleShower(*daughterParticle,evt,fPFParticleTag,fShowerTag);
  if(daughterTrack != 0x0){
    fdaughterIstrack[fNDAUGHTERS]                = 1;
    fdaughterIsshower[fNDAUGHTERS]               = 0;
    fdaughterTheta[fNDAUGHTERS]                  = daughterTrack->Theta();
    fdaughterPhi[fNDAUGHTERS]                    = daughterTrack->Phi();
    fdaughterLength[fNDAUGHTERS]                 = daughterTrack->Length();
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
    
  }
  else if(daughterShower != 0x0){
    fdaughterIstrack[fNDAUGHTERS]                   = 0;
    fdaughterIsshower[fNDAUGHTERS]                  = 1;
    fdaughterLength[fNDAUGHTERS]                    = daughterShower->Length();
    fdaughterShowerBestPlane[fNDAUGHTERS]           = daughterShower->best_plane();
    fdaughterOpeningAngle[fNDAUGHTERS]              = daughterShower->OpenAngle();
    fdaughterStartPosition[fNDAUGHTERS][0]          = daughterShower->ShowerStart().X();
    fdaughterStartPosition[fNDAUGHTERS][1]          = daughterShower->ShowerStart().Y();
    fdaughterStartPosition[fNDAUGHTERS][2]          = daughterShower->ShowerStart().Z();
    fdaughterStartDirection[fNDAUGHTERS][0]         = daughterShower->Direction().X();
    fdaughterStartDirection[fNDAUGHTERS][1]         = daughterShower->Direction().Y();
    fdaughterStartDirection[fNDAUGHTERS][2]         = daughterShower->Direction().Z();
  }
  fdaughterID[fNDAUGHTERS]                          = daughterID;
  // NHits associated with this pfParticle
  fdaughterNHits[fNDAUGHTERS]                       = (pfpUtil.GetPFParticleHits(*daughterParticle,evt,fPFParticleTag)).size();
  // T0
  std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*daughterParticle,evt,fPFParticleTag);
  if(!pfdaughterT0vec.empty())
    fdaughterT0[fNDAUGHTERS] = pfdaughterT0vec[0].Time();
  
  // Increment counter
  fNDAUGHTERS++;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void protoana::ProtoDUNEelectronAnaTree::Initialize(){
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  for(int k=0; k < 3; k++){
    fbeamtrackPos[k] = -999.0;
    fbeamtrackEndPos[k] = -999.0;
    fbeamtrackDir[k] = -999.0;
    fprimaryTruth_vtx[k]= -999.0; 
    fprimaryVertex[k] = -999.0;
    fdaughterVertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
  }
  for( int l=0; l<1000; l ++) {
     fbeamtrackPos_at[l][0] = -999.; 
     fbeamtrackPos_at[l][1] = -999.; 
     fbeamtrackPos_at[l][2] = -999.; 
     fbeamtrackMom_at[l][0] = -999.; 
     fbeamtrackMom_at[l][1] = -999.; 
     fbeamtrackMom_at[l][2] = -999.; 
     fbeamtrackMom_at[l][3] = -999.; 
  }
  fprimaryShower_nHits =0;
  fprimaryTruthShower_nHits =0;
  for( int m=0; m<MAXHits; m ++){
     fprimarydEdx[m]= -999.0;
     fprimarydQdx[m]= -999.0;
     fprimary_calX[m] =-999.0;
     fprimary_calY[m] =-999.0;
     fprimary_calZ[m] =-999.0;
     fprimary_cal_pitch[m]=-999.0;
     fprimaryResidualRange[m] = -999.0;
   
     fprimaryShower_hit_w[m] =-999.0;
     fprimaryShower_hit_q[m] =-999.0;
     fprimaryShower_hit_t[m] =-999.0;
     fprimaryShower_hit_X[m] =-999.0;
     fprimaryShower_hit_Y[m] =-999.0;
     fprimaryShower_hit_Z[m] =-999.0;
     fprimaryShower_hit_ch[m] =-999.0;
     fprimaryTruthShower_hit_w[m] =-999.0;
     fprimaryTruthShower_hit_q[m] =-999.0;
     fprimaryTruthShower_hit_t[m] =-999.0;
     fprimaryTruthShower_hit_X[m] =-999.0;
     fprimaryTruthShower_hit_Y[m] =-999.0;
     fprimaryTruthShower_hit_Z[m] =-999.0;
     fprimaryShower_hit_cnn[m] = -999.0;
  }
  fprimaryTruth_trkID =-999;
  fprimaryTruth_pdg = -999;
  fprimaryTruth_E = -999;
  fprimaryTruth_Edepo  = -999;
  fprimaryTruth_purity = -999;
  fprimarynCal = 0;
  fbeamCheckIsMatched = -999;
  fbeamtrigger = -999;
  ftof = -999.0;
  for(int k=0; k < 2; k++){
    fcerenkovPressure[k] = -999.0;
    fcerenkovTime[k] = -999.0; 
    fcerenkovStatus[k] = -999;
  }
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -999;
  fbeamtrackTime = -999.0;
  fbeamtrackID = -999;
  fbeam_ntrjPoints =0; 
  fbeamtrackNDaughters= 0;
  fprimaryIstrack = -999;
  fprimaryIsshower = -999;

  fprimaryBDTScore = -999.0;
  fprimaryNHits = 0;
  fprimaryTheta = -999.0;
  fprimaryPhi = -999.0;
  fprimaryLength = -999.0;
  fprimaryMomentum = -999.0;
  fprimaryEndMomentum = -999.0;
  fprimaryOpeningAngle = -999.0;
  fprimaryShowerBestPlane = -999;
  fprimaryShowerEnergy = -999.0;
  fprimaryShowerCharge = -999.0;
  fprimaryShowerTruth_Charge = -999.0;
  fprimaryShowerMIPEnergy = -999.0;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryMomentumByRangeMuon = -999.0;
  fprimaryT0 = -999.0;

  fNDAUGHTERS = 0;

  for(int k=0; k < NMAXDAUGTHERS; k++){
    fdaughterIstrack[k] = -999;
    fdaughterIsshower[k] = -999;
    fdaughterNHits[k] = -999;
    fdaughterTheta[k] = -999.0;
    fdaughterPhi[k] = -999.0;
    fdaughterLength[k] = -999.0;
    for(int l=0; l < 3; l++){
      fdaughterEndPosition[k][l] = -999.0;
      fdaughterStartPosition[k][l] = -999.0;
      fdaughterEndDirection[k][l] = -999.0;
      fdaughterStartDirection[k][l] = -999.0;

    }
    fdaughterOpeningAngle[k] = -999.0;
    fdaughterShowerBestPlane[k] = -999;
    fdaughterID[k] = -999;
    fdaughterT0[k] = -999;

  }

}

DEFINE_ART_MODULE(protoana::ProtoDUNEelectronAnaTree)
