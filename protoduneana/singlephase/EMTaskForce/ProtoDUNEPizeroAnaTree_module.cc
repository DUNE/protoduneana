////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEPizeroAnaTree
// File:        ProtoDUNEPizeroAnaTree_module.cc
//
// Extract protoDUNE useful information, do a first pre-selection and
// save output to a flat tree
//
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
// modified by Aaron Higuera ahiguera@central.uh.edu
// and Milo Vermeulen milov@nikhef.nl
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
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
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcore/Geometry/WireReadout.h"
// #include "larcore/Geometry/Geometry.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dunecalib/Calib/XYZCalib.h"
#include "dunecalib/CalibServices/XYZCalibService.h"
#include "dunecalib/CalibServices/XYZCalibServiceProtoDUNE.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "PiZeroProcess.h"

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

// Maximum number of beam particle daughters to save
const int NMAXDAUGHTERS = 25;
// // Maximum number of objects in APA3 to save
// const int NMAXOBJECTS = 100;
// Maximum number of primary MC pi0s to save
const int NMAXPIZEROS = 25;
// Maximum number of hits to save
const int NMAXHITS = 5000;

namespace protoana {
  class ProtoDUNEPizeroAnaTree;
}


class protoana::ProtoDUNEPizeroAnaTree : public art::EDAnalyzer {
public:

  explicit ProtoDUNEPizeroAnaTree(fhicl::ParameterSet const & p);

  ProtoDUNEPizeroAnaTree(ProtoDUNEPizeroAnaTree const &) = delete;
  ProtoDUNEPizeroAnaTree(ProtoDUNEPizeroAnaTree &&) = delete;
  ProtoDUNEPizeroAnaTree & operator = (ProtoDUNEPizeroAnaTree const &) = delete;
  ProtoDUNEPizeroAnaTree & operator = (ProtoDUNEPizeroAnaTree &&) = delete;

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
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;


  // Initialise tree variables
  void Initialise();
  void ResetPi0Vars();
  void FillPrimaryPFParticle(art::Event const & evt,
                             detinfo::DetectorClocksData const& clockData,
                             detinfo::DetectorPropertiesData const& detProp,
                             const recob::PFParticle* particle);
  void FillAPA3Object(art::Event const & evt, art::Handle<std::vector<recob::PFParticle>> pfpHandle);
  void setPiZeroInfo(art::Event const & evt,
                     detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     const std::vector<const simb::MCParticle*>& pi0s);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fSimulationTag;
  std::string fCalorimetryTag;
  calo::CalorimetryAlg fCalorimetryAlg;
  std::string fParticleIDTag;
  std::string fTrackTag;
  std::string fShowerTag;
  std::string fShowerCaloTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fHitTag;
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

  // Reconstructed primary particle
  double fprimaryVertex[3];
  int fprimaryIstrack;
  int fprimaryIsshower;
  double fprimaryBDTScore;
  double fprimaryCNNScore;
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
  int fprimaryIsBeamparticle;
  int    fprimaryTruth_trkID;
  int    fprimaryTruth_pdg;
  double fprimaryTruth_E;
  double fprimaryTruth_vtx[3];
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  int    fprimarynCal;
  double fprimarydEdx[NMAXHITS];
  double fprimarydQdx[NMAXHITS];
  double fprimary_cal_pos[NMAXHITS][3];
  double fprimary_cal_pitch[NMAXHITS];
  double fprimaryResidualRange[NMAXHITS];
  int    fprimaryShower_nHits; //collection only
  int    fprimaryShower_hit_w[NMAXHITS];
  double fprimaryShower_hit_q[NMAXHITS];
  double fprimaryShower_hit_t[NMAXHITS];
  double fprimaryShower_hit_pos[NMAXHITS][3];
  int fprimaryID;
  double fprimaryT0;
  int fprimaryNDaughters;

  // Primary particle daughters
  int fprimaryDaughterID[NMAXDAUGHTERS];
  int fprimaryDaughterIstrack[NMAXDAUGHTERS];
  int fprimaryDaughterIsshower[NMAXDAUGHTERS];
  double fprimaryDaughterBDTScore[NMAXDAUGHTERS];
  double fprimaryDaughterCNNScore[NMAXDAUGHTERS];
  int fprimaryDaughterNCollHits[NMAXDAUGHTERS];
  int fprimaryDaughterNHits[NMAXDAUGHTERS];
  double fprimaryDaughterEnergy[NMAXDAUGHTERS];
  double fprimaryDaughterEnergyFromHits[NMAXDAUGHTERS];
  double fprimaryDaughterMomentum[NMAXDAUGHTERS];
  double fprimaryDaughterEndMomentum[NMAXDAUGHTERS];
  double fprimaryDaughterLength[NMAXDAUGHTERS];
  double fprimaryDaughterStartPosition[NMAXDAUGHTERS][3];
  double fprimaryDaughterEndPosition[NMAXDAUGHTERS][3];
  double fprimaryDaughterStartDirection[NMAXDAUGHTERS][3];
  double fprimaryDaughterEndDirection[NMAXDAUGHTERS][3];
  int fprimaryDaughterParentPdg[NMAXDAUGHTERS];
  int fprimaryDaughterGrandparentPdg[NMAXDAUGHTERS];
  int fprimaryDaughterParentID[NMAXDAUGHTERS];
  int fprimaryDaughterGrandparentID[NMAXDAUGHTERS];
  double fprimaryDaughterParentE[NMAXDAUGHTERS];
  double fprimaryDaughterGrandparentE[NMAXDAUGHTERS];
  double fprimaryDaughterParentStart[NMAXDAUGHTERS][3];
  double fprimaryDaughterGrandparentStart[NMAXDAUGHTERS][3];
  double fprimaryDaughterParentEnd[NMAXDAUGHTERS][3];
  double fprimaryDaughterGrandparentEnd[NMAXDAUGHTERS][3];
  int fprimaryDaughterParentEndProcess[NMAXDAUGHTERS];
  // double fprimaryDaughterHitdQdx[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterCaloHitdQdx[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterCaloHitPitch[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterOwnHitPitch[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterHitCharge[NMAXDAUGHTERS][NMAXHITS];
  // double fprimaryDaughterHitdEdxArea[NMAXDAUGHTERS][NMAXHITS];
  int fprimaryDaughterHitWire[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterHitTime[NMAXDAUGHTERS][NMAXHITS];
  double fprimaryDaughterHitPosition[NMAXDAUGHTERS][NMAXHITS][3];
  // double fprimaryDaughterCaloHitPosition[NMAXDAUGHTERS][NMAXHITS][3];

  // // Instead of primary daughters, look at all objects in APA3
  // int fObjID[NMAXOBJECTS];
  // int fObjIsTrack[NMAXOBJECTS];
  // int fObjIsShower[NMAXOBJECTS];
  // int fObjIsPrimaryDaughter[NMAXOBJECTS];
  // double fObjBDTscore[NMAXOBJECTS];
  // double fObjCNNscore[NMAXOBJECTS];
  // double fObjMomentum[NMAXOBJECTS];
  // double fObjEndMomentum[NMAXOBJECTS];
  // double fObjLength[NMAXOBJECTS];
  // double fObjStartPosition[NMAXOBJECTS][3];
  // double fObjEndPosition[NMAXOBJECTS][3];
  // double fObjStartDirection[NMAXOBJECTS][3];
  // double fObjEndDirection[NMAXOBJECTS][3];
  // int fObjMCParentPdg[NMAXOBJECTS];
  // int fObjMCGrandparentPdg[NMAXOBJECTS];
  // int fObjMCParentID[NMAXOBJECTS];
  // int fObjMCGrandparentID[NMAXOBJECTS];
  // // Hit information
  // int fObjNumHits[NMAXOBJECTS];
  // double fObjHitPosition[NMAXOBJECTS][NMAXHITS][3];
  // double fObjHitTime[NMAXOBJECTS][NMAXHITS];
  // int fObjHitWire[NMAXOBJECTS][NMAXHITS];
  // double fObjHitCharge[NMAXOBJECTS][NMAXHITS];
  // double fObjHitPitch[NMAXOBJECTS][NMAXHITS];
  // double fObjHitdEdx[NMAXOBJECTS][NMAXHITS];
  // double fObjHitdQdx[NMAXOBJECTS][NMAXHITS];

  // MC pi0 variables
  int fMCPi0ID[NMAXPIZEROS];
  double fMCPi0Energy[NMAXPIZEROS];
  double fMCPi0StartPosition[NMAXPIZEROS][3];
  double fMCPi0EndPosition[NMAXPIZEROS][3];
  int fMCPhoton1ID[NMAXPIZEROS];
  int fMCPhoton2ID[NMAXPIZEROS];
  double fMCPhoton1Energy[NMAXPIZEROS];
  double fMCPhoton2Energy[NMAXPIZEROS];
  double fMCPhoton1StartPosition[NMAXPIZEROS][3];
  double fMCPhoton2StartPosition[NMAXPIZEROS][3];
  double fMCPhoton1EndPosition[NMAXPIZEROS][3];
  double fMCPhoton2EndPosition[NMAXPIZEROS][3];
  // Reco hits related to photons
  int fMCPhoton1NumHits[NMAXPIZEROS];
  int fMCPhoton2NumHits[NMAXPIZEROS];
  int fMCPhoton1_hit_w[NMAXPIZEROS][NMAXHITS];
  int fMCPhoton2_hit_w[NMAXPIZEROS][NMAXHITS];
  double fMCPhoton1_hit_t[NMAXPIZEROS][NMAXHITS];
  double fMCPhoton2_hit_t[NMAXPIZEROS][NMAXHITS];
  double fMCPhoton1_hit_q[NMAXPIZEROS][NMAXHITS];
  double fMCPhoton2_hit_q[NMAXPIZEROS][NMAXHITS];
  double fMCPhoton1_hit_pos[NMAXPIZEROS][NMAXHITS][3];
  double fMCPhoton2_hit_pos[NMAXPIZEROS][NMAXHITS][3];

  // Reco pi0 variables
  // Showers
  int fShower1ID[NMAXPIZEROS];
  int fShower2ID[NMAXPIZEROS];
  double fShower1StartPosition[NMAXPIZEROS][3];
  double fShower2StartPosition[NMAXPIZEROS][3];
  double fShower1Direction[NMAXPIZEROS][3];
  double fShower2Direction[NMAXPIZEROS][3];
  double fShower1Length[NMAXPIZEROS];
  double fShower2Length[NMAXPIZEROS];
  double fShower1Completeness[NMAXPIZEROS];
  double fShower2Completeness[NMAXPIZEROS];
  double fShower1Purity[NMAXPIZEROS];
  double fShower2Purity[NMAXPIZEROS];
  // Reco hits related to showers
  int fShower1NumHits[NMAXPIZEROS];
  double fShower1_cnn_sc[NMAXPIZEROS][NMAXHITS];
  double fShower2_cnn_sc[NMAXPIZEROS][NMAXHITS];
  int fShower2NumHits[NMAXPIZEROS];
  double fShower1_cal_pos[NMAXPIZEROS][NMAXHITS][3];
  double fShower2_cal_pos[NMAXPIZEROS][NMAXHITS][3];
  double fShower1_cal_E[NMAXPIZEROS];
  double fShower2_cal_E[NMAXPIZEROS];
  double fShower1Energy[NMAXPIZEROS];
  double fShower2Energy[NMAXPIZEROS];
  double fShower1EnergyFromHits[NMAXPIZEROS];
  double fShower2EnergyFromHits[NMAXPIZEROS];
  int fShower1HasBeamParent[NMAXPIZEROS];
  int fShower2HasBeamParent[NMAXPIZEROS];
  int fShower1HasBeamGrandparent[NMAXPIZEROS];
  int fShower2HasBeamGrandparent[NMAXPIZEROS];
  double fShower1_cal_pitch[NMAXPIZEROS][NMAXHITS];
  double fShower2_cal_pitch[NMAXPIZEROS][NMAXHITS];
  double fShower1_cal_dEdx[NMAXPIZEROS][NMAXHITS];
  double fShower2_cal_dEdx[NMAXPIZEROS][NMAXHITS];
  double fShower1_cal_dQdx[NMAXPIZEROS][NMAXHITS];
  double fShower2_cal_dQdx[NMAXPIZEROS][NMAXHITS];

  // Other
  //  geo::GeometryCore const* fGeometryService;
  geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
  spacecharge::SpaceChargeService::provider_type const* fSCEService;
};


protoana::ProtoDUNEPizeroAnaTree::ProtoDUNEPizeroAnaTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  beamlineUtil(p.get<fhicl::ParameterSet>("BeamLineUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fSimulationTag(p.get<std::string>("SimulationTag")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
  fTrackTag(p.get<std::string>("TrackTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fShowerCaloTag(p.get<std::string>("ShowerCalorimetryTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fVerbose(p.get<int>("Verbose")) {
  // Get a pointer to the geometry service provider.
  //  fGeometryService = lar::providerFrom<geo::Geometry>();
  fSCEService = lar::providerFrom<spacecharge::SpaceChargeService>();
}

void protoana::ProtoDUNEPizeroAnaTree::beginJob(){

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
  fPandoraBeam->Branch("primaryCNNScore",               &fprimaryCNNScore,              "primaryCNNScore/D");
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
  fPandoraBeam->Branch("primaryTruth_E",                &fprimaryTruth_E,               "primaryTruth_E/D");
  fPandoraBeam->Branch("primaryTruth_vtx",              &fprimaryTruth_vtx,             "primaryTruth_vtx[3]/D");
  fPandoraBeam->Branch("primaryTruth_pdg",              &fprimaryTruth_pdg,             "primaryTruth_pdg/I");
  fPandoraBeam->Branch("primaryTruth_trkID",            &fprimaryTruth_trkID,           "primaryTruth_trkID/I");
  fPandoraBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
  fPandoraBeam->Branch("primaryShowerCharge",           &fprimaryShowerCharge,          "primaryShowerCharge/D");
  fPandoraBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/D");
  fPandoraBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/D");

  fPandoraBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
  fPandoraBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
  fPandoraBeam->Branch("primarynCal",                   &fprimarynCal,                  "primarynCal/I");
  fPandoraBeam->Branch("primarydQdx",                   &fprimarydQdx,                  "primarydQdx[primarynCal]/D");
  fPandoraBeam->Branch("primary_cal_pos",               &fprimary_cal_pos,              "primary_cal_pos[primarynCal][3]/D");
  fPandoraBeam->Branch("primary_cal_pitch",             &fprimary_cal_pitch,            "primary_cal_pitch[primarynCal]/D");
  fPandoraBeam->Branch("primarydEdx",                   &fprimarydEdx,                  "primarydEdx[primarynCal]/D");
  fPandoraBeam->Branch("primaryResidualRange",          &fprimaryResidualRange,         "primaryResidualRange[primarynCal]/D");
  fPandoraBeam->Branch("primaryT0",                     &fprimaryT0,                    "primaryT0/D");
  fPandoraBeam->Branch("primaryNDaughters",             &fprimaryNDaughters,            "primaryNDaughters/I");

  fPandoraBeam->Branch("primaryDaughterID",               &fprimaryDaughterID,                "primaryDaughterID[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterIstrack",          &fprimaryDaughterIstrack,           "primaryDaughterIstrack[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterIsshower",         &fprimaryDaughterIsshower,          "primaryDaughterIsshower[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterBDTScore",         &fprimaryDaughterBDTScore,          "primaryDaughterBDTScore[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterCNNScore",         &fprimaryDaughterCNNScore,          "primaryDaughterCNNScore[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterNCollHits",        &fprimaryDaughterNCollHits,         "primaryDaughterNCollHits[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterNHits",            &fprimaryDaughterNHits,             "primaryDaughterNHits[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterEnergy",           &fprimaryDaughterEnergy,            "primaryDaughterEnergy[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterEnergyFromHits",   &fprimaryDaughterEnergyFromHits,    "primaryDaughterEnergyFromHits[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterMomentum",         &fprimaryDaughterMomentum,          "primaryDaughterMomentum[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterEndMomentum",      &fprimaryDaughterEndMomentum,       "primaryDaughterEndMomentum[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterLength",           &fprimaryDaughterLength,            "primaryDaughterLength[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterStartPosition",    &fprimaryDaughterStartPosition,     "primaryDaughterStartPosition[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterEndPosition",      &fprimaryDaughterEndPosition,       "primaryDaughterEndPosition[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterStartDirection",   &fprimaryDaughterStartDirection,    "primaryDaughterStartDirection[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterEndDirection",     &fprimaryDaughterEndDirection,      "primaryDaughterEndDirection[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterParentPdg",        &fprimaryDaughterParentPdg,         "primaryDaughterParentPdg[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterGrandparentPdg",   &fprimaryDaughterGrandparentPdg,    "primaryDaughterGrandparentPdg[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterParentID",         &fprimaryDaughterParentID,          "primaryDaughterParentID[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterGrandparentID",    &fprimaryDaughterGrandparentID,     "primaryDaughterGrandparentID[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterParentE",          &fprimaryDaughterParentE,           "primaryDaughterParentE[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterGrandparentE",     &fprimaryDaughterGrandparentE,      "primaryDaughterGrandparentE[primaryNDaughters]/D");
  fPandoraBeam->Branch("primaryDaughterParentStart",      &fprimaryDaughterParentStart,       "primaryDaughterParentStart[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterGrandparentStart", &fprimaryDaughterGrandparentStart,  "primaryDaughterGrandparentStart[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterParentEnd",        &fprimaryDaughterParentEnd,         "primaryDaughterParentEnd[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterGrandparentEnd",   &fprimaryDaughterGrandparentEnd,    "primaryDaughterGrandparentEnd[primaryNDaughters][3]/D");
  fPandoraBeam->Branch("primaryDaughterParentEndProcess", &fprimaryDaughterParentEndProcess,  "primaryDaughterParentEndProcess[primaryNDaughters]/I");
  fPandoraBeam->Branch("primaryDaughterHitPosition",      &fprimaryDaughterHitPosition,       ("primaryDaughterHitPosition[primaryNDaughters][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("primaryDaughterCaloHitPosition",  &fprimaryDaughterCaloHitPosition,   ("primaryDaughterCaloHitPosition[primaryNDaughters][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  fPandoraBeam->Branch("primaryDaughterHitCharge",        &fprimaryDaughterHitCharge,         ("primaryDaughterHitCharge[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("primaryDaughterHitdEdxArea",      &fprimaryDaughterHitdEdxArea,       ("primaryDaughterHitdEdxArea[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("primaryDaughterHitdQdx",          &fprimaryDaughterHitdQdx,           ("primaryDaughterHitdQdx[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fPandoraBeam->Branch("primaryDaughterCaloHitdQdx",      &fprimaryDaughterCaloHitdQdx,       ("primaryDaughterCaloHitdQdx[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fPandoraBeam->Branch("primaryDaughterCaloHitPitch",     &fprimaryDaughterCaloHitPitch,      ("primaryDaughterCaloHitPitch[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fPandoraBeam->Branch("primaryDaughterOwnHitPitch",      &fprimaryDaughterOwnHitPitch,       ("primaryDaughterOwnHitPitch[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());
  fPandoraBeam->Branch("primaryDaughterHitWire",          &fprimaryDaughterHitWire,           ("primaryDaughterHitWire[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/I").c_str());
  fPandoraBeam->Branch("primaryDaughterHitTime",          &fprimaryDaughterHitTime,           ("primaryDaughterHitTime[primaryNDaughters][" + std::to_string(NMAXHITS) + "]/D").c_str());

  // fPandoraBeam->Branch("ObjID",                &fObjID,                ("ObjID[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjIsTrack",           &fObjIsTrack,           ("ObjIsTrack[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjIsShower",          &fObjIsShower,          ("ObjIsShower[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjIsPrimaryDaughter", &fObjIsPrimaryDaughter, ("ObjIsPrimaryDaughter[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjBDTscore",          &fObjBDTscore,          ("ObjBDTscore[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjCNNscore",          &fObjCNNscore,          ("ObjCNNscore[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjMomentum",          &fObjMomentum,          ("ObjMomentum[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjEndMomentum",       &fObjEndMomentum,       ("ObjEndMomentum[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjLength",            &fObjLength,            ("ObjLength[" + std::to_string(NMAXOBJECTS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjStartPosition",     &fObjStartPosition,     ("ObjStartPosition[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("ObjEndPosition",       &fObjEndPosition,       ("ObjEndPosition[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("ObjStartDirection",    &fObjStartDirection,    ("ObjStartDirection[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("ObjEndDirection",      &fObjEndDirection,      ("ObjEndDirection[" + std::to_string(NMAXOBJECTS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("ObjParentPdg",         &fObjMCParentPdg,         ("ObjParentPdg[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjGrandparentPdg",    &fObjMCGrandparentPdg,    ("ObjGrandparentPdg[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjParentID",          &fObjMCParentID,          ("ObjParentID[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjGrandparentID",     &fObjMCGrandparentID,     ("ObjGrandparentID[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // // Object hit information.
  // fPandoraBeam->Branch("ObjNumHits",           &fObjNumHits,           ("ObjNumHits[" + std::to_string(NMAXOBJECTS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjHitPosition",       &fObjHitPosition,       ("ObjHitPosition[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("ObjHitTime",           &fObjHitTime,           ("ObjHitTime[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjHitWire",           &fObjHitWire,           ("ObjHitWire[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/I").c_str());
  // fPandoraBeam->Branch("ObjHitCharge",         &fObjHitCharge,         ("ObjHitCharge[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjHitPitch",          &fObjHitPitch,          ("ObjHitPitch[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjHitdEdx",           &fObjHitdEdx,           ("ObjHitdEdx[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("ObjHitdQdx",           &fObjHitdQdx,           ("ObjHitdQdx[" + std::to_string(NMAXOBJECTS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());

  // MC pi0 variables
  // fPandoraBeam->Branch("MCPi0ID",                       &fMCPi0ID,                      ("MCPi0ID[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("MCPi0Energy",                   &fMCPi0Energy,                  ("MCPi0Energy[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("MCPi0StartPosition",            &fMCPi0StartPosition,           ("MCPi0StartPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("MCPi0EndPosition",              &fMCPi0EndPosition,             ("MCPi0EndPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton1ID",                   &fMCPhoton1ID,                  ("MCPhoton1ID[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("MCPhoton2ID",                   &fMCPhoton2ID,                  ("MCPhoton2ID[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("MCPhoton1Energy",               &fMCPhoton1Energy,              ("MCPhoton1Energy[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton2Energy",               &fMCPhoton2Energy,              ("MCPhoton2Energy[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton1StartPosition",        &fMCPhoton1StartPosition,       ("MCPhoton1StartPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton2StartPosition",        &fMCPhoton2StartPosition,       ("MCPhoton2StartPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton1EndPosition",          &fMCPhoton1EndPosition,         ("MCPhoton1EndPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton2EndPosition",          &fMCPhoton2EndPosition,         ("MCPhoton2EndPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // // Reco hits related to photons
  // fPandoraBeam->Branch("MCPhoton1NumHits",              &fMCPhoton1NumHits,             ("MCPhoton1NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("MCPhoton2NumHits",              &fMCPhoton2NumHits,             ("MCPhoton2NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("MCPhoton1_hit_w",               &fMCPhoton1_hit_w,              ("MCPhoton1_hit_w[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/I").c_str());
  // fPandoraBeam->Branch("MCPhoton2_hit_w",               &fMCPhoton2_hit_w,              ("MCPhoton2_hit_w[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/I").c_str());
  // fPandoraBeam->Branch("MCPhoton1_hit_q",               &fMCPhoton1_hit_q,              ("MCPhoton1_hit_q[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton2_hit_q",               &fMCPhoton2_hit_q,              ("MCPhoton2_hit_q[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton1_hit_t",               &fMCPhoton1_hit_t,              ("MCPhoton1_hit_t[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton2_hit_t",               &fMCPhoton2_hit_t,              ("MCPhoton2_hit_t[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton1_hit_pos",             &fMCPhoton1_hit_pos,            ("MCPhoton1_hit_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("MCPhoton2_hit_pos",             &fMCPhoton2_hit_pos,            ("MCPhoton2_hit_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());

  // // Shower variables
  // fPandoraBeam->Branch("Shower1ID",                     &fShower1ID,                    ("Shower1ID[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower2ID",                     &fShower2ID,                    ("Shower2ID[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower1StartPosition",          &fShower1StartPosition,         ("Shower1StartPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("Shower2StartPosition",          &fShower2StartPosition,         ("Shower2StartPosition[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("Shower1Direction",              &fShower1Direction,             ("Shower1Direction[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("Shower2Direction",              &fShower2Direction,             ("Shower2Direction[" + std::to_string(NMAXPIZEROS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("Shower1Length",                 &fShower1Length,                ("Shower1Length[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2Length",                 &fShower2Length,                ("Shower2Length[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1Completeness",           &fShower1Completeness,          ("Shower1Completeness[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2Completeness",           &fShower2Completeness,          ("Shower2Completeness[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1Purity",                 &fShower1Purity,                ("Shower1Purity[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2Purity",                 &fShower2Purity,                ("Shower2Purity[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // // Reco hits related to showers
  // fPandoraBeam->Branch("Shower1NumHits",                &fShower1NumHits,               ("Shower1NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower2NumHits",                &fShower2NumHits,               ("Shower2NumHits[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower1_cnn_sc",                &fShower1_cnn_sc,               ("Shower1_cnn_sc[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2_cnn_sc",                &fShower2_cnn_sc,               ("Shower2_cnn_sc[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1_cal_pos",               &fShower1_cal_pos,              ("Shower1_cal_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("Shower2_cal_pos",               &fShower2_cal_pos,              ("Shower2_cal_pos[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "][3]/D").c_str());
  // fPandoraBeam->Branch("Shower1_cal_E",                 &fShower1_cal_E,                ("Shower1_cal_E[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2_cal_E",                 &fShower2_cal_E,                ("Shower2_cal_E[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1Energy",                 &fShower1Energy,                ("Shower1Energy[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2Energy",                 &fShower2Energy,                ("Shower2Energy[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1EnergyFromHits",         &fShower1EnergyFromHits,        ("Shower1EnergyFromHits[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2EnergyFromHits",         &fShower2EnergyFromHits,        ("Shower2EnergyFromHits[" + std::to_string(NMAXPIZEROS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1HasBeamParent",          &fShower1HasBeamParent,         ("Shower1HasBeamParent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower2HasBeamParent",          &fShower2HasBeamParent,         ("Shower2HasBeamParent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower1HasBeamGrandparent",     &fShower1HasBeamGrandparent,    ("Shower1HasBeamGrandparent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower2HasBeamGrandparent",     &fShower2HasBeamGrandparent,    ("Shower2HasBeamGrandparent[" + std::to_string(NMAXPIZEROS) + "]/I").c_str());
  // fPandoraBeam->Branch("Shower1_cal_pitch",             &fShower1_cal_pitch,            ("Shower1_cal_pitch[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2_cal_pitch",             &fShower2_cal_pitch,            ("Shower2_cal_pitch[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1_cal_dEdx",              &fShower1_cal_dEdx,             ("Shower1_cal_dEdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2_cal_dEdx",              &fShower2_cal_dEdx,             ("Shower2_cal_dEdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower1_cal_dQdx",              &fShower1_cal_dQdx,             ("Shower1_cal_dQdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
  // fPandoraBeam->Branch("Shower2_cal_dQdx",              &fShower2_cal_dQdx,             ("Shower2_cal_dQdx[" + std::to_string(NMAXPIZEROS) + "][" + std::to_string(NMAXHITS) + "]/D").c_str());
}

void protoana::ProtoDUNEPizeroAnaTree::analyze(art::Event const & evt){
  // Initialise tree parameters
  Initialise();
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

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  //check for reco pandora stuff
  auto recoParticleHandle = evt.getHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  if( !recoParticleHandle ) return;

  // Get all of the PFParticles, by default from the "pandora" product
  // auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  for(const recob::PFParticle* particle : pfParticles){

    FillPrimaryPFParticle(evt, clockData, detProp, particle);
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackTag);
    // std::cout << "Primary PFParticle ID: " << particle->Self() << '\n';

    fprimaryVertex[0] = vtx.X(); fprimaryVertex[1] = vtx.Y(); fprimaryVertex[2] = vtx.Z();
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackTag);

    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  }

  // // Look at all the reconstructed particles that are fully contained in APA3.
  // FillAPA3Object(evt, recoParticleHandle);

  bool beamTriggerEvent = false;
  if(!evt.isRealData()){
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    // Define and fill a handle to point to a vector of the MCParticles
    auto MCParticleHandle = evt.getHandle<std::vector<simb::MCParticle>>(fSimulationTag);
    if (!MCParticleHandle) {
      // Handle no simb::MCParticles.
      throw cet::exception("ProtoDUNEPizeroAnaTree")
          << " No simb::MCParticle objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    std::map<int, const simb::MCParticle*> MCPmap;
    for(const simb::MCParticle& mcp : *MCParticleHandle) {
      MCPmap[mcp.TrackId()] = &mcp;
    }

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
      fbeam_ntrjPoints = std::min(geantGoodParticle->NumberTrajectoryPoints(), (unsigned)1000);
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

      std::cout << "Primary beam particle PDG code: " << fbeamtrackPdg << " and E: "
                << fbeamtrackEnergy << '\n';

      std::vector<const simb::MCParticle*> pi0s;

      // Search for daughter pi0s if the primary particle is a pion.
      if(abs(fbeamtrackPdg) == 211) {
        std::cout << "Found a pion! (ID " << fbeamtrackID << ", PDG "
                  << fbeamtrackPdg << ")\n";

        std::cout << "Pion has pi0 daughters:\n";
        std::vector<const simb::MCParticle*> daughters;
        daughters.push_back(geantGoodParticle);
        while(daughters.size()) {
          const simb::MCParticle* parent = daughters[0];
          for(int di = 0; di < parent->NumberDaughters(); ++di) {
            const int daughter_ID = parent->Daughter(di);
            const simb::MCParticle* daughter = MCPmap[daughter_ID];
            daughters.push_back(daughter);
            if(daughter->PdgCode() == 111) {
              std::cout << "  ID " << daughter_ID << ", PDG " << daughter->PdgCode()
                        << ", P " << daughter->P() << '\n';
              pi0s.push_back(daughter);
            }
          }
          daughters.erase(daughters.begin());
        }
      }

      // Record MC and shower information
      setPiZeroInfo(evt, clockData, detProp, pi0s);

    } //geantGoodParticle
  } // MC
  else { //data
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
        // if TOFChan == -1, then there was not a successful match, if it's
        // 0, 1, 2, or 3, then there was a good match corresponding to the
        // different pair-wise combinations of the upstream and downstream
        // channels
        if(beaminfo[0]->GetTOFChan() != -1){
          ftof = beaminfo[0]->GetTOF();
        }
        // If possible PIDs contain pion, set PDG code to pion.
        std::vector<int> PIDs = beamlineUtil.GetPID(evt, 6);
        if(std::find(PIDs.begin(), PIDs.end(), 211) != PIDs.end() && fprimaryIsshower != 1) {
          fbeamtrackPdg = 211;
        } else if(std::find(PIDs.begin(), PIDs.end(), 2212) != PIDs.end() && fprimaryIsshower != 1) {
          fbeamtrackPdg = 2212;
        } else if(std::find(PIDs.begin(), PIDs.end(), 321) != PIDs.end() && fprimaryIsshower != 1) {
          fbeamtrackPdg = 321;
        }
        std::cout << "Set PDG to " << fbeamtrackPdg << '\n';
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

        // Record pi0s reconstructed from showers
        // Make recob::Showers accessible
        auto rShowerHandle = evt.getHandle<std::vector<recob::Shower>>(fShowerTag);
        if (!rShowerHandle) {
          // Handle no recob::Showers.
          throw cet::exception("PiZeroExtract")
              << " No recob::Shower objects in this event - "
              << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        // // Loop through showers
        // for(const recob::Shower& shower : *rShowerHandle) {
        //   pizero::PiZeroProcess recopzproc(shower, evt, fShowerTag);
        //   // Only record if the showers come somewhat close
        //   if(!recopzproc.allRecoSet()) continue;
        //   const double shdist = pizero::ClosestDistance(recopzproc.shower1(), recopzproc.shower2());
        //   if(shdist > 0 && shdist < 100) {
        //     setPiZeroInfo(evt, recopzproc);
        //     break;
        //   }
        // }
      } //good beam trigger
    } //good beaminfo
  }//for data

  // Fill tree if event came from beam trigger.
  if(beamTriggerEvent)   fPandoraBeam->Fill();

  // Put some space between events.
  std::cout << "\n\n";
}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEPizeroAnaTree::endJob(){

}

// -----------------------------------------------------------------------------
void
protoana::ProtoDUNEPizeroAnaTree::FillPrimaryPFParticle(art::Event const & evt,
                                                        detinfo::DetectorClocksData const& clockData,
                                                        detinfo::DetectorPropertiesData const& detProp,
                                                        const recob::PFParticle* particle){

  // Pandora's BDT beam-cosmic score
  fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);

  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
  const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle, evt,fPFParticleTag,fTrackTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track>>(fTrackTag);
  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower>>(fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromTracks(recoTracks,evt,fTrackTag);
  // Set some common variables that have different access methods.
  std::vector<art::Ptr<recob::Hit>> pfpHits;
  std::vector<anab::Calorimetry> calovector;
  const simb::MCParticle* mcparticle = NULL;
  if(thisTrack != 0x0) {
    pfpHits = findHitsFromTracks.at(thisTrack->ID());
    calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackTag, fCalorimetryTag);
    mcparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, evt, fTrackTag);
  } else if(thisShower != 0x0) {
    pfpHits = findHitsFromShowers.at(thisShower->ID());
    calovector = showerUtil.GetRecoShowerCalorimetry(*thisShower, evt, fShowerTag, fShowerCaloTag);
    mcparticle = truthUtil.GetMCParticleFromRecoShower(clockData, *thisShower, evt, fShowerTag);
  }
  // NHits associated with this pfParticle.
  fprimaryNHits = pfpHits.size();

  if(mcparticle){
    fprimaryTruth_pdg = mcparticle->PdgCode();
    fprimaryTruth_trkID = mcparticle->TrackId();
    fprimaryTruth_vtx[0] = mcparticle->Vx();
    fprimaryTruth_vtx[1] = mcparticle->Vy();
    fprimaryTruth_vtx[2] = mcparticle->Vx();
    fprimaryTruth_E = mcparticle->E();
    if( fbeamtrackID != -999 && fbeamtrackID == fprimaryTruth_trkID ) {
      fprimaryIsBeamparticle = 1;
    }
  }

  // Calorimetry only collection plane
  for(size_t k = 0; k < calovector.size() && k<3; k++){
    int plane = calovector[k].PlaneID().Plane;
    if(plane !=2 ) continue;
    fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
    fprimaryRange[plane]              = calovector[k].Range();
    fprimarynCal = std::min(calovector[k].dEdx().size(), (size_t)NMAXHITS);
    for(size_t l=0; l<calovector[k].dEdx().size() && l<NMAXHITS; ++l){
      fprimarydEdx[l]= calovector[k].dEdx()[l];
      fprimarydQdx[l]= calovector[k].dQdx()[l];
      fprimaryResidualRange[l]= calovector[k].ResidualRange()[l];
      const auto &pos=(calovector[k].XYZ())[l];
      fprimary_cal_pos[l][0] = pos.X();
      fprimary_cal_pos[l][1] = pos.Y();
      fprimary_cal_pos[l][2] = pos.Z();
      fprimary_cal_pitch[l] =calovector[k].TrkPitchVec()[l];
    }
  }

  // Get the CNN score for the primary particle.
  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );
  fprimaryCNNScore = 0;
  for(const art::Ptr<recob::Hit>& hit : pfpHits) {
    std::array<float,4> cnn_out = hitResults.getOutput( hit );
    const double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
    const double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
    fprimaryCNNScore += cnn_score;
  }
  fprimaryCNNScore /= fprimaryNHits;
  std::cout << "Got primary CNN score " << fprimaryCNNScore << ", true PDG code: " << fprimaryTruth_pdg << '\n';
  std::cout << "Pandora track: " << thisTrack << ", shower: " << thisShower << '\n';

  // Get the T0 for this pfParticle
  std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
  if(!pfT0vec.empty())
    fprimaryT0 = pfT0vec[0].Time();

  if(thisTrack != 0x0){

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

  } // end is track
  else if(thisShower != 0x0){

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

  } // end is shower

  // Record children of primary track
  auto recoParticleHandle = evt.getHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  if( !recoParticleHandle ) return;

  fprimaryNDaughters = std::min(particle->NumDaughters(), NMAXDAUGHTERS);
  for(int di = 0; di < fprimaryNDaughters; ++di) {
    fprimaryDaughterID[di] = particle->Daughter(di);
    const recob::PFParticle* daughter = &recoParticleHandle->at(fprimaryDaughterID[di]);
    if(daughter == 0x0) continue;

    const recob::Track* daughterTrack   = pfpUtil.GetPFParticleTrack(*daughter, evt,fPFParticleTag,fTrackTag);
    const recob::Shower* daughterShower = pfpUtil.GetPFParticleShower(*daughter,evt,fPFParticleTag,fShowerTag);

    // Pandora's BDT beam-cosmic score
    fprimaryDaughterBDTScore[di] = (double)pfpUtil.GetBeamCosmicScore(*daughter,evt,fPFParticleTag);
    // Hits and calorimetry associated with this pfParticle
    std::vector<art::Ptr<recob::Hit>> pfpDHits;
    std::vector<anab::Calorimetry> dcalovector;
    if(daughterTrack != 0x0) {
      pfpDHits = findHitsFromTracks.at( daughterTrack->ID() );
      dcalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackTag, fCalorimetryTag);
    } else if(daughterShower != 0x0) {
      pfpDHits = findHitsFromShowers.at( daughterShower->ID() );
      dcalovector = showerUtil.GetRecoShowerCalorimetry(*daughterShower, evt, fShowerTag, fShowerCaloTag);
    }
    if(dcalovector.size() != 3 && fVerbose > 0) {
      std::cerr << "WARNING::Calorimetry vector size for track is = " << dcalovector.size() << std::endl;
    }

    // Loop over calorimetry.
    std::map<int, double> calo_hit_dQdx;
    std::map<int, double> calo_hit_pitch;
    for(size_t k = 0; k < dcalovector.size() && k<3; k++){
      int plane = dcalovector[k].PlaneID().Plane;
      if(plane !=2 ) continue;
      for(size_t l=0; l<dcalovector[k].dEdx().size() && l<NMAXHITS; ++l){
        // fprimaryDaughterHitdEdx[di][l]= dcalovector[k].dEdx()[l];
        // fprimaryDaughterCaloHitdQdx[di][l]= dcalovector[k].dQdx()[l];
        calo_hit_dQdx[dcalovector[k].TpIndices()[l]] = dcalovector[k].dQdx()[l];
        calo_hit_pitch[dcalovector[k].TpIndices()[l]] = dcalovector[k].TrkPitchVec()[l];
        // tmp_dQdx.push_back(fprimaryDaughterCaloHitdQdx[di][l]);
        // fprimaryDaughterResidualRange[l]= dcalovector[k].ResidualRange()[l];
        const geo::Point_t &pos=(dcalovector[k].XYZ())[l];
        if(daughterShower != 0) {
          const TVector3 hitvec(pos.X(), pos.Y(), pos.Z());
        }
        // fprimaryDaughterCaloHitPosition[di][l][0] = pos.X();
        // fprimaryDaughterCaloHitPosition[di][l][1] = pos.Y();
        // fprimaryDaughterCaloHitPosition[di][l][2] = pos.Z();
        // fprimaryDaughterHitPitch[di][l] = dcalovector[k].TrkPitchVec()[l];
      } // Loop over calorimetry hits.
    } // Loop over calorimetry.

    // Loop over hits, record Aidan's CNN shower-like score.
    art::FindManyP<recob::SpacePoint> spFromHits(pfpDHits, evt, fHitTag);
    fprimaryDaughterCNNScore[di] = 0;
    unsigned rec_hiti = 0;
    for( unsigned i = 0; i < pfpDHits.size() && rec_hiti < NMAXHITS; ++i){
      if(pfpDHits[i]->WireID().Plane != 2) continue;

      const geo::WireGeo* pwire = fWireReadoutGeom.WirePtr(pfpDHits[i]->WireID());
      auto xyzWire = pwire->GetCenter();
      fprimaryDaughterHitWire[di][rec_hiti] = pfpDHits[i]->WireID().Wire;
      fprimaryDaughterHitTime[di][rec_hiti] = pfpDHits[i]->PeakTime();
      fprimaryDaughterHitPosition[di][rec_hiti][0] = detProp.ConvertTicksToX(pfpDHits[i]->PeakTime(),pfpDHits[i]->WireID().Plane,pfpDHits[i]->WireID().TPC,0);
      fprimaryDaughterHitPosition[di][rec_hiti][2] = xyzWire.Z();
      std::vector<art::Ptr<recob::SpacePoint>> sps = spFromHits.at(i);
      if(!sps.empty()) {
        fprimaryDaughterHitPosition[di][rec_hiti][1] = sps[0]->XYZ()[1];
      }
      TVector3 hitpos(fprimaryDaughterHitPosition[di][rec_hiti][0],
                      fprimaryDaughterHitPosition[di][rec_hiti][1],
                      fprimaryDaughterHitPosition[di][rec_hiti][2]);

      std::array<float,4> cnn_out = hitResults.getOutput( pfpDHits[i] );
      const double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
      const double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
      fprimaryDaughterCNNScore[di] += cnn_score;
      // Record general hit information.
      fprimaryDaughterHitCharge[di][rec_hiti] = pfpDHits[i]->Integral();
      // Determine own hit pitch with SCE corrections.
      const double wirePitch = fWireReadoutGeom.Plane(pfpDHits[i]->WireID()).WirePitch();
      TVector3 objdir(0,0,0);
      if(daughterTrack != 0x0) {
        objdir = {daughterTrack->Trajectory().StartDirection().X(),
                  daughterTrack->Trajectory().StartDirection().Y(),
                  daughterTrack->Trajectory().StartDirection().Z()};
      } else if(daughterShower != 0x0) {
        objdir = daughterShower->Direction();
      }
      TVector3 hitdir = objdir.Unit();
      if(objdir.Mag2() != 0 && hitpos[1] != -999.0) {
        // Here we have the object's direction and full hit position.
        // Set hit end to point on next wire in object's direction.
        objdir = objdir.Unit();
        TVector3 hitend = hitpos + objdir * (wirePitch/objdir.Z());
        if(fSCEService->EnableCalSpatialSCE()) {
          // Do SCE corrections on hitpos and hitend to get proper hit direction.
          // std::cout << "Doing SCE corrections.\n";
          const geo::Vector_t poscorr = fSCEService->GetCalPosOffsets(geo::Point_t(hitpos), pfpDHits[i]->WireID().TPC);
          const geo::Vector_t endcorr = fSCEService->GetCalPosOffsets(geo::Point_t(hitend), pfpDHits[i]->WireID().TPC);
          hitpos += TVector3(poscorr.X(), poscorr.Y(), poscorr.Z());
          hitend += TVector3(endcorr.X(), endcorr.Y(), endcorr.Z());
        }
        hitdir = (hitend-hitpos).Unit();
      } //else if(objdir.Mag2() != 0) {
      //   std::cout << "Hit position incomplete: " << hitpos[0] << ',' << hitpos[1] << ',' << hitpos[2] << '\n';
      // } else if(hitpos[1] != -999) {
      //   std::cout << "No object direction.\n";
      // } else { std::cout << "??????\n"; }
      fprimaryDaughterOwnHitPitch[di][rec_hiti] = wirePitch / abs(hitdir.Z());
      fprimaryDaughterCaloHitdQdx[di][rec_hiti] = calo_hit_dQdx[pfpDHits[i].key()];
      fprimaryDaughterCaloHitPitch[di][rec_hiti] = calo_hit_pitch[pfpDHits[i].key()];
      // std::cout << "Bare dx: " << wirePitch / abs(objdir.Z())
      //           << ", own dx: " << fprimaryDaughterOwnHitPitch[di][rec_hiti]
      //           << ", calo dx: " << fprimaryDaughterCaloHitPitch[di][rec_hiti] << '\n';
      ++rec_hiti;
    } // Loop over hits.
    fprimaryDaughterNCollHits[di] = rec_hiti;
    fprimaryDaughterNHits[di] = pfpDHits.size();
    fprimaryDaughterCNNScore[di] /= rec_hiti;

    if(daughterTrack != 0x0) {
      fprimaryDaughterIstrack[di] = 1;
      fprimaryDaughterIsshower[di] = 0;

      fprimaryDaughterMomentum[di]           = daughterTrack->StartMomentum();
      fprimaryDaughterEndMomentum[di]        = daughterTrack->EndMomentum();
      fprimaryDaughterLength[di]             = daughterTrack->Length();
      fprimaryDaughterStartPosition[di][0]   = daughterTrack->Trajectory().Start().X();
      fprimaryDaughterStartPosition[di][1]   = daughterTrack->Trajectory().Start().Y();
      fprimaryDaughterStartPosition[di][2]   = daughterTrack->Trajectory().Start().Z();
      fprimaryDaughterEndPosition[di][0]     = daughterTrack->Trajectory().End().X();
      fprimaryDaughterEndPosition[di][1]     = daughterTrack->Trajectory().End().Y();
      fprimaryDaughterEndPosition[di][2]     = daughterTrack->Trajectory().End().Z();
      fprimaryDaughterStartDirection[di][0]  = daughterTrack->Trajectory().StartDirection().X();
      fprimaryDaughterStartDirection[di][1]  = daughterTrack->Trajectory().StartDirection().Y();
      fprimaryDaughterStartDirection[di][2]  = daughterTrack->Trajectory().StartDirection().Z();
      fprimaryDaughterEndDirection[di][0]    = daughterTrack->Trajectory().EndDirection().X();
      fprimaryDaughterEndDirection[di][1]    = daughterTrack->Trajectory().EndDirection().Y();
      fprimaryDaughterEndDirection[di][2]    = daughterTrack->Trajectory().EndDirection().Z();
    } else if(daughterShower != 0x0) {
      fprimaryDaughterIstrack[di] = 0;
      fprimaryDaughterIsshower[di] = 1;

      std::vector<const recob::Hit*> shHits = showerUtil.GetRecoShowerHits(*daughterShower, evt, fShowerTag);
      fprimaryDaughterEnergyFromHits[di]     = showerUtil.EstimateEnergyFromHitCharge(clockData,
                                                                                      detProp,
                                                                                      shHits,
                                                                                      fCalorimetryAlg)[2];
      fprimaryDaughterLength[di]             = daughterShower->Length();
      fprimaryDaughterStartDirection[di][0]   = daughterShower->Direction().X();
      fprimaryDaughterStartDirection[di][1]   = daughterShower->Direction().Y();
      fprimaryDaughterStartDirection[di][2]   = daughterShower->Direction().Z();
      fprimaryDaughterStartPosition[di][0]    = daughterShower->ShowerStart().X();
      fprimaryDaughterStartPosition[di][1]    = daughterShower->ShowerStart().Y();
      fprimaryDaughterStartPosition[di][2]    = daughterShower->ShowerStart().Z();
    }

    // Attempt to get the (grand)parent truth information of this daughter.
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const simb::MCParticle* parent = truthUtil.GetMCParticleFromReco(clockData, *daughter, evt, fPFParticleTag);
    if(parent != 0x0) {
      fprimaryDaughterParentPdg[di]          = parent->PdgCode();
      // std::cout << "True parent PDG: " << fprimaryDaughterParentPdg[di] << '\n';
      fprimaryDaughterParentID[di]          = parent->TrackId();
      // std::cout << "Daughter parent: " << fprimaryDaughterParentID[di] << '\n';
      fprimaryDaughterParentE[di]          = parent->E();
      // std::cout << "True parent E: " << fprimaryDaughterParentE[di] << '\n';
      if(parent->EndProcess() == "conv") {
        fprimaryDaughterParentEndProcess[di] = 0;
      } else if(parent->EndProcess() == "phot") {
        fprimaryDaughterParentEndProcess[di] = 1;
      }
      parent->Position().Vect().GetXYZ(fprimaryDaughterParentStart[di]);
      parent->EndPosition().Vect().GetXYZ(fprimaryDaughterParentEnd[di]);
      if(parent->Mother() != 0) {
        const simb::MCParticle* gparent = pi_serv->TrackIdToParticle_P(parent->Mother());
        fprimaryDaughterGrandparentPdg[di]   =  gparent->PdgCode();
        fprimaryDaughterGrandparentID[di]   =  gparent->TrackId();
        fprimaryDaughterGrandparentE[di]   =  gparent->E();
        // std::cout << "True grandparent PDG: " << fprimaryDaughterGrandparentPdg[di] << '\n';
        gparent->Position().Vect().GetXYZ(fprimaryDaughterGrandparentStart[di]);
        gparent->EndPosition().Vect().GetXYZ(fprimaryDaughterGrandparentEnd[di]);
      }
    } // If MCParticle parent.
  } // for PFParticle daughters

} // FillPrimaryPFParticle

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// void protoana::ProtoDUNEPizeroAnaTree::FillAPA3Object(art::Event const & evt, art::Handle<std::vector<recob::PFParticle>> pfpHandle) {
//   // Put all PFParticles in a map by ID.
//   std::unordered_map<size_t, const recob::PFParticle*> PFPmap;
//   for(const recob::PFParticle& pfp : *pfpHandle) {
//     PFPmap[pfp.Self()] = &pfp;
//   }
//   // Primary particle ID.
//   const size_t beamID = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag).size() > 0? pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag)[0]->Self(): -1;
//   // Keep track of the number of saved objects
//   unsigned obji = 0;
//   for(const recob::PFParticle& pfp : *pfpHandle) {
//     if(obji >= NMAXOBJECTS) break;
//
//     // Determine whether particle is classified as track or shower.
//     const recob::Track* track   = pfpUtil.GetPFParticleTrack(pfp, evt,fPFParticleTag,fTrackTag);
//     const recob::Shower* shower = pfpUtil.GetPFParticleShower(pfp,evt,fPFParticleTag,fShowerTag);
//
//     // Do initial check whether object is fully contained in APA3.
//     if(track != 0x0) {
//       const int TPCstart = fGeometryService->FindTPCAtPosition(track->Start()).deepestIndex();
//       const int TPCend = fGeometryService->FindTPCAtPosition(track->End()).deepestIndex();
//       // APA3 == TPCID 1??
//       if(TPCstart != 1 || TPCend != 1) continue;
//     } else if(shower != 0x0) {
//       const int TPCstart = fGeometryService->FindTPCAtPosition(shower->ShowerStart()).deepestIndex();
//       const int TPCend = fGeometryService->FindTPCAtPosition(shower->ShowerStart() + shower->Length() * shower->Direction()).deepestIndex();
//       // APA3 == TPCID 1??
//       if(TPCstart != 1 || TPCend != 1) continue;
//     } else {
//       continue;
//     }
//
//     // General PFParticle attributes
//     fObjID[obji] = pfp.Self();
//     if(PFPmap[pfp.Parent()] != 0x0 && PFPmap[pfp.Parent()]->Self() == beamID) {
//       fObjIsPrimaryDaughter[obji] = 1;
//       // std::cout << "Yes, ID " << fObjID[obji] << '\n';
//     } else {
//       // std::cout << "No, ID " << fObjID[obji] << '\n';
//       fObjIsPrimaryDaughter[obji] = 0;
//     }
//
//     // Pandora's BDT beam-cosmic score
//     fObjBDTscore[obji] = (double)pfpUtil.GetBeamCosmicScore(pfp,evt,fPFParticleTag);
//
//     // Get associations between PFParticles, clusters and hits.
//     const std::vector<const recob::Cluster*> pfpClusters = pfpUtil.GetPFParticleClusters(pfp,evt,fPFParticleTag);
//     const auto allClusters = evt.getValidHandle<std::vector<recob::Cluster>>(fPFParticleTag);
//     const art::FindManyP<recob::Hit> findHits(allClusters, evt, fPFParticleTag);
//     std::vector<art::Ptr<recob::Hit>> pfpHits;
//     for(auto cluster : pfpClusters){
//       const std::vector<art::Ptr<recob::Hit>> clusterHits = findHits.at(cluster->ID());
//       for(art::Ptr<recob::Hit> hit : clusterHits){
//         pfpHits.push_back(hit);
//       }
//     }
//     // NHits associated with this pfParticle
//     fObjNumHits[obji] = pfpHits.size();
//
//     // Aidan's CNN track-like score: determine the average for this object and save.
//     anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );
//     fObjCNNscore[obji] = 0;
//     // Loop over PFP hits.
//     for( unsigned i = 0; i < pfpHits.size() && i < NMAXHITS; ++i){
//       // Collection plane only.
//       if(pfpHits[i]->WireID().Plane != 2) continue;
//       // Record CNN score per hit.
//       const std::array<float,4> cnn_out = hitResults.getOutput( pfpHits[i] );
//       const double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
//       const double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
//       fObjCNNscore[obji] += cnn_score;
//       // Record other hit information.
//       fObjHitTime[obji][i] = pfpHits[i]->PeakTime();
//       fObjHitWire[obji][i] = pfpHits[i]->WireID().Wire;
//       fObjHitCharge[obji][i] = pfpHits[i]->Integral();
//     }
//     fObjCNNscore[obji] /= fObjNumHits[obji];
//     // std::cout << "CNN score: " << fObjCNNscore[obji] << ". Pandora says:\n";
//
//     // Distinctive track/shower information.
//     const simb::MCParticle* mcparticle;
//     if(track != 0) {
//       // PFParticle is marked as track by Pandora.
//       mcparticle = truthUtil.GetMCParticleFromReco(*track, evt, fTrackTag);
//
//       fObjIsShower[obji]                    = 0;
//       fObjIsTrack[obji]                     = 1;
//       fObjLength[obji]                      = track->Length();
//       fObjStartPosition[obji][0]            = track->Start().X();
//       fObjStartPosition[obji][1]            = track->Start().Y();
//       fObjStartPosition[obji][2]            = track->Start().Z();
//       fObjEndPosition[obji][0]              = track->End().X();
//       fObjEndPosition[obji][1]              = track->End().Y();
//       fObjEndPosition[obji][2]              = track->End().Z();
//       fObjStartDirection[obji][0]           = track->StartDirection().X();
//       fObjStartDirection[obji][1]           = track->StartDirection().Y();
//       fObjStartDirection[obji][2]           = track->StartDirection().Z();
//       fObjEndDirection[obji][0]             = track->EndDirection().X();
//       fObjEndDirection[obji][1]             = track->EndDirection().Y();
//       fObjEndDirection[obji][2]             = track->EndDirection().Z();
//       fObjMomentum[obji]                    = track->StartMomentum();
//       fObjEndMomentum[obji]                 = track->EndMomentum();
//
//       // Calorimetry only collection plane
//       std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*track, evt, fTrackTag, fCalorimetryTag);
//       if(calovector.size() != 3 && fVerbose > 0)
//         std::cerr << "WARNING::Calorimetry vector size for track is = " << calovector.size() << std::endl;
//
//       for(size_t k = 0; k < calovector.size() && k<3; k++){
//         int plane = calovector[k].PlaneID().Plane;
//         if(plane !=2 ) continue;
//         for(size_t l=0; l<calovector[k].dEdx().size() && l<NMAXHITS; ++l){
//            fObjHitdEdx[obji][l]= calovector[k].dEdx()[l];
//            fObjHitdQdx[obji][l]= calovector[k].dQdx()[l];
//            // fObjResidualRange[l]= calovector[k].ResidualRange()[l];
//            const auto &pos=(calovector[k].XYZ())[l];
//            fObjHitPosition[obji][l][0] = pos.X();
//            fObjHitPosition[obji][l][1] = pos.Y();
//            fObjHitPosition[obji][l][2] = pos.Z();
//            fObjHitPitch[obji][l] = calovector[k].TrkPitchVec()[l];
//         }
//       }
//     } else if(shower != 0) {
//       // PFParticle is marked as shower by Pandora.
//       mcparticle = truthUtil.GetMCParticleFromReco(*shower, evt, fShowerTag);
//
//       fObjIsTrack[obji]                     = 0;
//       fObjIsShower[obji]                    = 1;
//       fObjLength[obji]                      = shower->Length();
//       // fObjShowerBestPlane[obji]          = shower->best_plane();
//       // fObjOpeningAngle[obji]             = shower->OpenAngle();
//       fObjStartPosition[obji][0]            = shower->ShowerStart().X();
//       fObjStartPosition[obji][1]            = shower->ShowerStart().Y();
//       fObjStartPosition[obji][2]            = shower->ShowerStart().Z();
//       const TVector3 shEnd = shower->ShowerStart() + shower->Length()*shower->Direction();
//       fObjEndPosition[obji][0]              = shEnd.X();
//       fObjEndPosition[obji][1]              = shEnd.Y();
//       fObjEndPosition[obji][2]              = shEnd.Z();
//       fObjStartDirection[obji][0]           = shower->Direction().X();
//       fObjStartDirection[obji][1]           = shower->Direction().Y();
//       fObjStartDirection[obji][2]           = shower->Direction().Z();
//       fObjEndDirection[obji][0]             = shower->Direction().X();
//       fObjEndDirection[obji][1]             = shower->Direction().Y();
//       fObjEndDirection[obji][2]             = shower->Direction().Z();
//
//       // Calorimetry only collection plane
//       std::vector<anab::Calorimetry> calovector = showerUtil.GetRecoShowerCalorimetry(*shower, evt, fShowerTag, fShowerCaloTag);
//       if(calovector.size() != 3 && fVerbose > 0)
//         std::cerr << "WARNING::Calorimetry vector size for shower is = " << calovector.size() << std::endl;
//
//       for(size_t k = 0; k < calovector.size() && k<3; k++){
//         int plane = calovector[k].PlaneID().Plane;
//         if(plane !=2 ) continue;
//         for(size_t l=0; l<calovector[k].dEdx().size() && l<NMAXHITS; ++l){
//            fObjHitdEdx[obji][l]= calovector[k].dEdx()[l];
//            fObjHitdQdx[obji][l]= calovector[k].dQdx()[l];
//            // fObjResidualRange[l]= calovector[k].ResidualRange()[l];
//            const auto &pos=(calovector[k].XYZ())[l];
//            fObjHitPosition[obji][l][0] = pos.X();
//            fObjHitPosition[obji][l][1] = pos.Y();
//            fObjHitPosition[obji][l][2] = pos.Z();
//            fObjHitPitch[obji][l] = calovector[k].TrkPitchVec()[l];
//         }
//       }
//     } else {
//       continue;
//     }
//     // MCParticle matching
//     if(mcparticle){
//       art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
//       fObjMCParentPdg[obji] = mcparticle->PdgCode();
//       fObjMCParentID[obji] = mcparticle->TrackId();
//       if(mcparticle->Mother() > 0 && pi_serv->TrackIdToParticle_P(mcparticle->Mother()) != 0x0) {
//         fObjMCGrandparentPdg[obji] = pi_serv->TrackIdToParticle_P(mcparticle->Mother())->PdgCode();
//         fObjMCGrandparentID[obji] = pi_serv->TrackIdToParticle_P(mcparticle->Mother())->TrackId();
//       }
//     }
//
//     // int fObjID[NMAXOBJECTS];
//     // int fObjIsTrack[NMAXOBJECTS];
//     // int fObjIsShower[NMAXOBJECTS];
//     // int fObjIsPrimaryDaughter[NMAXOBJECTS];
//     // double fObjBDTscore[NMAXOBJECTS];
//     // double fObjCNNscore[NMAXOBJECTS];
//     // double fObjMomentum[NMAXOBJECTS];
//     // double fObjEndMomentum[NMAXOBJECTS];
//     // double fObjLength[NMAXOBJECTS];
//     // double fObjStartPosition[NMAXOBJECTS][3];
//     // double fObjEndPosition[NMAXOBJECTS][3];
//     // double fObjStartDirection[NMAXOBJECTS][3];
//     // double fObjEndDirection[NMAXOBJECTS][3];
//     // int fObjMCParentPdg[NMAXOBJECTS];
//     // int fObjMCGrandparentPdg[NMAXOBJECTS];
//     // int fObjMCParentID[NMAXOBJECTS];
//     // int fObjMCGrandparentID[NMAXOBJECTS];
//     // // Hit information
//     // int fObjNumHits[NMAXOBJECTS];
//     // double fObjHitPosition[NMAXOBJECTS][NMAXHITS][3];
//     // double fObjHitTime[NMAXOBJECTS][NMAXHITS];
//     // int fObjHitWire[NMAXOBJECTS][NMAXHITS];
//     // double fObjHitCharge[NMAXOBJECTS][NMAXHITS];
//     // double fObjHitPitch[NMAXOBJECTS][NMAXHITS];
//     // double fObjHitdEdx[NMAXOBJECTS][NMAXHITS];
//     // double fObjHitdQdx[NMAXOBJECTS][NMAXHITS];
//
//     // Only increment the counter if the object was contained in APA3.
//     ++obji;
//   }
// }

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void protoana::ProtoDUNEPizeroAnaTree::setPiZeroInfo(art::Event const & evt,
                                                     detinfo::DetectorClocksData const& clockData,
                                                     detinfo::DetectorPropertiesData const& detProp,
                                                     const std::vector<const simb::MCParticle*>& pi0s) {
  // Clean slate
  ResetPi0Vars();

  // Loop over all pi0s in the array by index.
  for(unsigned pi0i = 0; pi0i < NMAXPIZEROS && pi0i < pi0s.size(); ++pi0i) {
    pizero::PiZeroProcess pzproc(*(pi0s[pi0i]), evt, clockData, detProp, fShowerTag);

    // Set MC object variables
    if(pzproc.allMCSet()) {

      const simb::MCParticle* pi0 = pzproc.pi0();
      const simb::MCParticle* photon1 = pzproc.photon1();
      const simb::MCParticle* photon2 = pzproc.photon2();

      // MC pi0 variables
      fMCPi0ID[pi0i] = pi0->TrackId();
      fMCPi0Energy[pi0i] = pi0->E();
      pi0->Position().Vect().GetXYZ(fMCPi0StartPosition[pi0i]);
      pi0->EndPosition().Vect().GetXYZ(fMCPi0EndPosition[pi0i]);
      fMCPhoton1ID[pi0i] = photon1->TrackId();
      fMCPhoton2ID[pi0i] = photon2->TrackId();
      fMCPhoton1Energy[pi0i] = photon1->E();
      fMCPhoton2Energy[pi0i] = photon2->E();
      photon1->Position().Vect().GetXYZ(fMCPhoton1StartPosition[pi0i]);
      photon2->Position().Vect().GetXYZ(fMCPhoton2StartPosition[pi0i]);
      photon1->EndPosition().Vect().GetXYZ(fMCPhoton1EndPosition[pi0i]);
      photon2->EndPosition().Vect().GetXYZ(fMCPhoton2EndPosition[pi0i]);
      // Reco hits related to photons
      std::vector<const recob::Hit*> ph1hits = truthUtil.GetMCParticleHits(clockData, *photon1, evt, fHitTag);
      std::vector<const recob::Hit*> ph2hits = truthUtil.GetMCParticleHits(clockData, *photon2, evt, fHitTag);
      // Photon 1 hits
      art::FindManyP<recob::SpacePoint> ph1sps(ph1hits, evt, fPFParticleTag);
      unsigned phi1 = 0;
      for(unsigned i=0; i<ph1hits.size(); ++i) {
        if(ph1hits[i]->WireID().Plane != 2) continue;
        fMCPhoton1_hit_w[pi0i][phi1] = ph1hits[phi1]->WireID().Wire;
        fMCPhoton1_hit_t[pi0i][phi1] = ph1hits[phi1]->PeakTime();
        fMCPhoton1_hit_q[pi0i][phi1] = ph1hits[phi1]->Integral();
        std::vector<art::Ptr<recob::SpacePoint>> sps = ph1sps.at(i);
        if(!sps.empty()) {
          fMCPhoton1_hit_pos[pi0i][phi1][0] = sps[0]->XYZ()[0];
          fMCPhoton1_hit_pos[pi0i][phi1][1] = sps[0]->XYZ()[1];
          fMCPhoton1_hit_pos[pi0i][phi1][2] = sps[0]->XYZ()[2];
        }
        ++phi1;
      }
      // Photon 2 hits
      art::FindManyP<recob::SpacePoint> ph2sps(ph2hits, evt, fPFParticleTag);
      unsigned phi2 = 0;
      for(unsigned i=0; i<ph2hits.size(); ++i) {
        if(ph2hits[i]->WireID().Plane != 2) continue;
        fMCPhoton2_hit_w[pi0i][phi2] = ph2hits[phi2]->WireID().Wire;
        fMCPhoton2_hit_t[pi0i][phi2] = ph2hits[phi2]->PeakTime();
        fMCPhoton2_hit_q[pi0i][phi2] = ph2hits[phi2]->Integral();
        std::vector<art::Ptr<recob::SpacePoint>> sps = ph2sps.at(i);
        if(!sps.empty()) {
          fMCPhoton2_hit_pos[pi0i][phi2][0] = sps[0]->XYZ()[0];
          fMCPhoton2_hit_pos[pi0i][phi2][1] = sps[0]->XYZ()[1];
          fMCPhoton2_hit_pos[pi0i][phi2][2] = sps[0]->XYZ()[2];
        }
        ++phi2;
      }
      fMCPhoton1NumHits[pi0i] = phi1;
      fMCPhoton2NumHits[pi0i] = phi2;

    } // if MC set

    //This is how we get cnn score for now
    auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
    auto recoTracks = evt.getValidHandle< std::vector< recob::Track > >(fTrackTag);
    art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);
    art::FindManyP<recob::Hit> findHitsFromTracks(recoTracks,evt,fTrackTag);
    anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );

    // First shower
    const recob::Shower* shower1 = pzproc.shower1();
    if(shower1 != 0x0) {
      fShower1ID[pi0i] = showerUtil.GetShowerIndex(*shower1, evt, fShowerTag);
      std::vector< art::Ptr< recob::Hit > > sh1_hits = findHitsFromShowers.at(fShower1ID[pi0i]);
      for( size_t i=0; i<sh1_hits.size(); ++i){
         std::array<float,4> cnn_out = hitResults.getOutput( sh1_hits[i] );
         double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
         double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
         fShower1_cnn_sc[pi0i][i] = cnn_score;
      }
      shower1->ShowerStart().GetXYZ(fShower1StartPosition[pi0i]);
      shower1->Direction().GetXYZ(fShower1Direction[pi0i]);
      fShower1Length[pi0i] = shower1->Length();
      unsigned numMChits = truthUtil.GetMCParticleHits(clockData, *pzproc.photon1(), evt, fHitTag).size();
      unsigned showerHits = showerUtil.GetRecoShowerHits(*shower1, evt, fShowerTag).size();
      unsigned sharedHits = truthUtil.GetSharedHits(clockData, *pzproc.photon1(), *shower1, evt, fShowerTag, true).size();
      fShower1Completeness[pi0i] = (double)sharedHits/numMChits;
      fShower1Purity[pi0i] = (double)sharedHits/showerHits;
      // Find out whether (grand)parent of shower was primary.
      // Go through all PFParticles until the right one is found.
      const recob::PFParticle* shpfp = 0x0;
      std::vector<const recob::PFParticle*> beamPFParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
      auto PFParticleHandle = evt.getHandle<std::vector<recob::PFParticle> >(fPFParticleTag);
      if (PFParticleHandle && beamPFParticles.size() != 0) {
        const recob::PFParticle* beampfp = beamPFParticles[0];
        for(const recob::PFParticle& pfp : *PFParticleHandle) {
          const recob::Shower* thish = pfpUtil.GetPFParticleShower(pfp,evt,fPFParticleTag,fShowerTag);
          if(thish == 0x0) continue;
          // Test whether it's the right shower by comparing lengths.
          if(thish->Length() - shower1->Length() < 1e-5) {
            shpfp = &pfp;
            break;
          }
        }
        // Is PFP's parent or PFP's grandparent the beam particle?
        const recob::PFParticle* parent = 0x0;
        const recob::PFParticle* gparent = 0x0;
        if(shpfp && shpfp->Parent() != 0 && shpfp->Parent() < PFParticleHandle->size())
          parent = &PFParticleHandle->at(shpfp->Parent());
        if(parent && parent->Parent() != 0 && parent->Parent() < PFParticleHandle->size()) gparent = &PFParticleHandle->at(parent->Parent());

        fShower1HasBeamParent[pi0i] = parent && parent->Self() == beampfp->Self()? 1: 0;
        fShower1HasBeamGrandparent[pi0i] = gparent && gparent->Self() == beampfp->Self()? 1: 0;
      } // if get pfparticles and beampfparticles

      // Calorimetry related to showers
      std::vector<anab::Calorimetry> sh1calo = showerUtil.GetRecoShowerCalorimetry(*shower1, evt, fShowerTag, fShowerCaloTag);
      for(unsigned k=0; k<sh1calo.size(); ++k) {
        if(sh1calo[k].PlaneID().Plane != 2) continue;

        fShower1NumHits[pi0i] = std::min((int)sh1calo[k].dEdx().size(), NMAXHITS);
        fShower1_cal_E[pi0i] = sh1calo[k].KineticEnergy();
        std::vector<const recob::Hit*> shHitPtrs = showerUtil.GetRecoShowerHits(
                    *shower1, evt, fShowerTag);
        fShower1EnergyFromHits[pi0i] = showerUtil.EstimateEnergyFromHitCharge(clockData,
                                                                              detProp,
                                                                              shHitPtrs,
                                                                              fCalorimetryAlg)[2];

        for(int i=0; i<fShower1NumHits[pi0i]; ++i) {
          const geo::Point_t& pos = sh1calo[k].XYZ()[i];
          fShower1_cal_pos[pi0i][i][0] = pos.X();
          fShower1_cal_pos[pi0i][i][1] = pos.Y();
          fShower1_cal_pos[pi0i][i][2] = pos.Z();
          fShower1_cal_pitch[pi0i][i] = sh1calo[k].TrkPitchVec()[i];
          fShower1_cal_dEdx[pi0i][i] = sh1calo[k].dEdx()[i];
          fShower1_cal_dQdx[pi0i][i] = sh1calo[k].dQdx()[i];
        }
      }
    } // if shower1 != 0x0
    // Second shower
    const recob::Shower* shower2 = pzproc.shower2();
    if(shower2 != 0x0) {
      // Reco pi0 variables
      fShower2ID[pi0i] = showerUtil.GetShowerIndex(*shower2, evt, fShowerTag);
      std::vector< art::Ptr< recob::Hit > > sh2_hits = findHitsFromShowers.at(fShower2ID[pi0i]);
      for( size_t i=0; i<sh2_hits.size(); ++i){
         std::array<float,4> cnn_out = hitResults.getOutput( sh2_hits[i] );
         double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ];
         double cnn_score = cnn_out[ hitResults.getIndex("em") ]/p_trk_or_sh;
         fShower2_cnn_sc[pi0i][i] = cnn_score;
      }
      shower2->ShowerStart().GetXYZ(fShower2StartPosition[pi0i]);
      shower2->Direction().GetXYZ(fShower2Direction[pi0i]);
      fShower2Length[pi0i] = shower2->Length();
      unsigned numMChits = truthUtil.GetMCParticleHits(clockData, *pzproc.photon2(), evt, fHitTag).size();
      unsigned showerHits = showerUtil.GetRecoShowerHits(*shower2, evt, fShowerTag).size();
      unsigned sharedHits = truthUtil.GetSharedHits(clockData, *pzproc.photon2(), *shower2, evt, fShowerTag, true).size();
      fShower2Completeness[pi0i] = (double)sharedHits/numMChits;
      fShower2Purity[pi0i] = (double)sharedHits/showerHits;
      // Find out whether (grand)parent of shower was primary.
      // Go through all PFParticles until the right one is found.
      const recob::PFParticle* shpfp = 0x0;
      std::vector<const recob::PFParticle*> beamPFParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
      auto PFParticleHandle = evt.getHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
      if (PFParticleHandle && beamPFParticles.size() != 0) {
        const recob::PFParticle* beampfp = beamPFParticles[0];
        for(const recob::PFParticle& pfp : *PFParticleHandle) {
          const recob::Shower* thish = pfpUtil.GetPFParticleShower(pfp,evt,fPFParticleTag,fShowerTag);
          if(thish == 0x0) continue;
          // Test whether it's the right shower by comparing lengths.
          if(thish->Length() - shower2->Length() < 1e-5) {
            shpfp = &pfp;
            break;
          }
        }
        // Is PFP's parent or PFP's grandparent the beam particle?
        const recob::PFParticle* parent = 0x0;
        const recob::PFParticle* gparent = 0x0;
        if(shpfp && shpfp->Parent() != 0 && shpfp->Parent() < PFParticleHandle->size()) parent = &PFParticleHandle->at(shpfp->Parent());
        if(parent && parent->Parent() != 0 && parent->Parent() < PFParticleHandle->size()) gparent = &PFParticleHandle->at(parent->Parent());

        fShower2HasBeamParent[pi0i] = parent && parent->Self() == beampfp->Self()? 1: 0;
        fShower2HasBeamGrandparent[pi0i] = gparent && gparent->Self() == beampfp->Self()? 1: 0;
      } // if get pfparticles and beampfparticles

      // Calorimetry related to showers
      std::vector<anab::Calorimetry> sh2calo = showerUtil.GetRecoShowerCalorimetry(*shower2, evt, fShowerTag, fShowerCaloTag);
      for(unsigned k=0; k<sh2calo.size(); ++k) {
        if(sh2calo[k].PlaneID().Plane != 2) continue;

        fShower2NumHits[pi0i] = std::min((int)sh2calo[k].dEdx().size(), NMAXHITS);
        fShower2_cal_E[pi0i] = sh2calo[k].KineticEnergy();
        std::vector<const recob::Hit*> shHitPtrs = showerUtil.GetRecoShowerHits(
                    *shower2, evt, fShowerTag);
        fShower2EnergyFromHits[pi0i] = showerUtil.EstimateEnergyFromHitCharge(clockData, detProp,
                                                                              shHitPtrs, fCalorimetryAlg)[2];
        for(int i=0; i<fShower2NumHits[pi0i]; ++i) {
          const geo::Point_t& pos = sh2calo[k].XYZ()[i];
          fShower2_cal_pos[pi0i][i][0] = pos.X();
          fShower2_cal_pos[pi0i][i][1] = pos.Y();
          fShower2_cal_pos[pi0i][i][2] = pos.Z();
          fShower2_cal_pitch[pi0i][i] = sh2calo[k].TrkPitchVec()[i];
          fShower2_cal_dEdx[pi0i][i] = sh2calo[k].dEdx()[i];
          fShower2_cal_dQdx[pi0i][i] = sh2calo[k].dQdx()[i];
        }
      }
    } // if shower2 != 0x0
  } // Loop over pi0s

} // setPiZeroInfo


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void protoana::ProtoDUNEPizeroAnaTree::Initialise(){
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
  for( int m=0; m<NMAXHITS; m ++){
     fprimarydEdx[m]= -999.0;
     fprimarydQdx[m]= -999.0;
     fprimary_cal_pos[m][0] = -999.0;
     fprimary_cal_pos[m][1] = -999.0;
     fprimary_cal_pos[m][2] = -999.0;
     fprimary_cal_pitch[m]=-999.0;
     fprimaryResidualRange[m] = -999.0;

     fprimaryShower_hit_w[m] =-999.0;
     fprimaryShower_hit_q[m] =-999.0;
     fprimaryShower_hit_t[m] =-999.0;
     fprimaryShower_hit_pos[m][0] =-999.0;
     fprimaryShower_hit_pos[m][1] =-999.0;
     fprimaryShower_hit_pos[m][2] =-999.0;
  }
  fprimaryTruth_trkID =-999;
  fprimaryTruth_pdg = -999;
  fprimaryTruth_E = -999;
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
  fprimaryCNNScore = -999.0;
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
  fprimaryShowerMIPEnergy = -999.0;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryT0 = -999.0;
  fprimaryNDaughters = 0;

  for(unsigned k = 0; k < NMAXDAUGHTERS; ++k) {
    fprimaryDaughterID[k] = -999;
    fprimaryDaughterIstrack[k] = -999;
    fprimaryDaughterIsshower[k] = -999;
    fprimaryDaughterBDTScore[k] = -999.0;
    fprimaryDaughterCNNScore[k] = -999.0;
    fprimaryDaughterNCollHits[k] = 0;
    fprimaryDaughterNHits[k] = 0;
    fprimaryDaughterEnergy[k] = -999.0;
    fprimaryDaughterEnergyFromHits[k] = -999.0;
    fprimaryDaughterMomentum[k] = -999.0;
    fprimaryDaughterEndMomentum[k] = -999.0;
    fprimaryDaughterLength[k] = -999.0;
    for(unsigned l = 0; l < 3; ++l) {
      fprimaryDaughterStartPosition[k][l] = -999.0;
      fprimaryDaughterEndPosition[k][l] = -999.0;
      fprimaryDaughterStartDirection[k][l] = -999.0;
      fprimaryDaughterEndDirection[k][l] = -999.0;
    }
    fprimaryDaughterParentPdg[k] = -999;
    fprimaryDaughterGrandparentPdg[k] = -999;
    fprimaryDaughterParentID[k] = -999;
    fprimaryDaughterGrandparentID[k] = -999;
    fprimaryDaughterParentE[k] = -999.0;
    fprimaryDaughterGrandparentE[k] = -999.0;
    fprimaryDaughterParentStart[k][0] = -999.0;
    fprimaryDaughterParentStart[k][1] = -999.0;
    fprimaryDaughterParentStart[k][2] = -999.0;
    fprimaryDaughterGrandparentStart[k][0] = -999.0;
    fprimaryDaughterGrandparentStart[k][1] = -999.0;
    fprimaryDaughterGrandparentStart[k][2] = -999.0;
    fprimaryDaughterParentEnd[k][0] = -999.0;
    fprimaryDaughterParentEnd[k][1] = -999.0;
    fprimaryDaughterParentEnd[k][2] = -999.0;
    fprimaryDaughterGrandparentEnd[k][0] = -999.0;
    fprimaryDaughterGrandparentEnd[k][1] = -999.0;
    fprimaryDaughterGrandparentEnd[k][2] = -999.0;
    fprimaryDaughterParentEndProcess[k] = -999.0;
    for(unsigned l = 0; l < NMAXHITS; ++l) {
      // fprimaryDaughterHitdQdx[k][l] = -999.0;
      fprimaryDaughterCaloHitdQdx[k][l] = -999.0;
      fprimaryDaughterCaloHitPitch[k][l] = -999.0;
      fprimaryDaughterOwnHitPitch[k][l] = -999.0;
      fprimaryDaughterHitCharge[k][l] = -999.0;
      // fprimaryDaughterHitdEdxArea[k][l] = -999.0;
      fprimaryDaughterHitWire[k][l] = -999;
      fprimaryDaughterHitTime[k][l] = -999.0;
      fprimaryDaughterHitPosition[k][l][0] = -999.0;
      fprimaryDaughterHitPosition[k][l][1] = -999.0;
      fprimaryDaughterHitPosition[k][l][2] = -999.0;
      // fprimaryDaughterCaloHitPosition[k][l][0] = -999.0;
      // fprimaryDaughterCaloHitPosition[k][l][1] = -999.0;
      // fprimaryDaughterCaloHitPosition[k][l][2] = -999.0;
    }
  }

  // // Instead of primary daughters, look at all objects in APA3
  // for(int obji = 0; obji < NMAXOBJECTS; ++obji) {
  //   fObjID[obji] = -999;
  //   fObjIsTrack[obji] = -999;
  //   fObjIsShower[obji] = -999;
  //   fObjIsPrimaryDaughter[obji] = -999;
  //   fObjBDTscore[obji] = -999.0;
  //   fObjCNNscore[obji] = -999.0;
  //   fObjMomentum[obji] = -999.0;
  //   fObjEndMomentum[obji] = -999.0;
  //   fObjLength[obji] = -999.0;
  //   for(int pi = 0; pi < 3; ++pi) {
  //     fObjStartPosition[obji][pi] = -999.0;
  //     fObjEndPosition[obji][pi] = -999.0;
  //     fObjStartDirection[obji][pi] = -999.0;
  //     fObjEndDirection[obji][pi] = -999.0;
  //   }
  //   fObjMCParentPdg[obji] = -999;
  //   fObjMCGrandparentPdg[obji] = -999;
  //   fObjMCParentID[obji] = -999;
  //   fObjMCGrandparentID[obji] = -999;
  //   // Hit information
  //   fObjNumHits[obji] = -999;
  //   for(int hiti = 0; hiti < NMAXHITS; ++hiti) {
  //     fObjHitPosition[obji][hiti][0] = -999.0;
  //     fObjHitPosition[obji][hiti][1] = -999.0;
  //     fObjHitPosition[obji][hiti][2] = -999.0;
  //     fObjHitTime[obji][hiti] = -999.0;
  //     fObjHitWire[obji][hiti] = -999;
  //     fObjHitCharge[obji][hiti] = -999.0;
  //     fObjHitPitch[obji][hiti] = -999.0;
  //     fObjHitdEdx[obji][hiti] = -999.0;
  //     fObjHitdQdx[obji][hiti] = -999.0;
  //   }
  // }

  ResetPi0Vars();
} // Initialise

void protoana::ProtoDUNEPizeroAnaTree::ResetPi0Vars() {
  for(int pi0i = 0; pi0i < NMAXPIZEROS; ++pi0i) {
    // MC pi0 variables
    fMCPi0ID[pi0i] = -999;
    fMCPi0Energy[pi0i] = -999.0;
    for(int i=0; i<3; ++i) {
      fMCPi0StartPosition[pi0i][i] = -999.0;
      fMCPi0EndPosition[pi0i][i] = -999.0;
    }
    fMCPhoton1ID[pi0i] = -999;
    fMCPhoton2ID[pi0i] = -999;
    fMCPhoton1Energy[pi0i] = -999.0;
    fMCPhoton2Energy[pi0i] = -999.0;
    for(int i=0; i<3; ++i) {
      fMCPhoton1StartPosition[pi0i][i] = -999.0;
      fMCPhoton2StartPosition[pi0i][i] = -999.0;
      fMCPhoton1EndPosition[pi0i][i] = -999.0;
      fMCPhoton2EndPosition[pi0i][i] = -999.0;
    }
    // Reco hits related to photons
    fMCPhoton1NumHits[pi0i] = -999;
    fMCPhoton2NumHits[pi0i] = -999;
    for(int i=0; i < NMAXHITS; ++i) {
      fMCPhoton1_hit_w[pi0i][i] = -999;
      fMCPhoton2_hit_w[pi0i][i] = -999;
      fMCPhoton1_hit_t[pi0i][i] = -999.0;
      fMCPhoton2_hit_t[pi0i][i] = -999.0;
      fMCPhoton1_hit_q[pi0i][i] = -999.0;
      fMCPhoton2_hit_q[pi0i][i] = -999.0;
      fMCPhoton1_hit_pos[pi0i][i][0] = -999.0;
      fMCPhoton2_hit_pos[pi0i][i][0] = -999.0;
      fMCPhoton1_hit_pos[pi0i][i][1] = -999.0;
      fMCPhoton2_hit_pos[pi0i][i][1] = -999.0;
      fMCPhoton1_hit_pos[pi0i][i][2] = -999.0;
      fMCPhoton2_hit_pos[pi0i][i][2] = -999.0;
    }

    // Reco pi0 variables
    fShower1ID[pi0i] = -999;
    fShower2ID[pi0i] = -999;
    for(int i=0; i<3; ++i) {
      fShower1StartPosition[pi0i][i] = -999.0;
      fShower2StartPosition[pi0i][i] = -999.0;
      fShower1Direction[pi0i][i] = -999.0;
      fShower2Direction[pi0i][i] = -999.0;
    }
    fShower1Length[pi0i] = -999.0;
    fShower2Length[pi0i] = -999.0;
    fShower1Completeness[pi0i] = -999.0;
    fShower2Completeness[pi0i] = -999.0;
    fShower1Purity[pi0i] = -999.0;
    fShower2Purity[pi0i] = -999.0;
    // Reco hits related to showers
    fShower1NumHits[pi0i] = -999;
    fShower2NumHits[pi0i] = -999;
    fShower1_cal_E[pi0i] = -999.0;
    fShower2_cal_E[pi0i] = -999.0;
    fShower1Energy[pi0i] = -999.0;
    fShower2Energy[pi0i] = -999.0;
    fShower1EnergyFromHits[pi0i] = -999.0;
    fShower2EnergyFromHits[pi0i] = -999.0;
    fShower1HasBeamParent[pi0i] = -999;
    fShower2HasBeamParent[pi0i] = -999;
    fShower1HasBeamGrandparent[pi0i] = -999;
    fShower2HasBeamGrandparent[pi0i] = -999;
    for(int i=0; i<NMAXHITS; ++i) {
      fShower1_cnn_sc[pi0i][i] = -999.0;
      fShower2_cnn_sc[pi0i][i] = -999.0;
      fShower1_cal_pos[pi0i][i][0] = -999.0;
      fShower2_cal_pos[pi0i][i][0] = -999.0;
      fShower1_cal_pos[pi0i][i][1] = -999.0;
      fShower2_cal_pos[pi0i][i][1] = -999.0;
      fShower1_cal_pos[pi0i][i][2] = -999.0;
      fShower2_cal_pos[pi0i][i][2] = -999.0;
      fShower1_cal_pitch[pi0i][i] = -999.0;
      fShower2_cal_pitch[pi0i][i] = -999.0;
      fShower1_cal_dEdx[pi0i][i] = -999.0;
      fShower2_cal_dEdx[pi0i][i] = -999.0;
      fShower1_cal_dQdx[pi0i][i] = -999.0;
      fShower2_cal_dQdx[pi0i][i] = -999.0;
    }
  } // Loop over NMAXPIZEROS
} // ResetPi0Vars

DEFINE_ART_MODULE(protoana::ProtoDUNEPizeroAnaTree)
