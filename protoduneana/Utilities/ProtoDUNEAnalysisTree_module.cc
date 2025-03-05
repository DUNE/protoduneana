////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEAnalysisTree
// File:        ProtoDUNEAnalysisTree_module.cc
//
// Extract protoDUNE useful information, do a first pre-selection and save output to a flat tree
// 
// Some parts are copied from the beam module example
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
//
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

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
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
const int NMAXDAUGTHERS = 15;
const int NMAXTRUTHDAUGTHERS = 30;

namespace protoana {
  class ProtoDUNEAnalysisTree;
}

// -----------------------------------------------------------------------------
class protoana::ProtoDUNEAnalysisTree : public art::EDAnalyzer {
public:

  explicit ProtoDUNEAnalysisTree(fhicl::ParameterSet const & p);

  ProtoDUNEAnalysisTree(ProtoDUNEAnalysisTree const &) = delete;
  ProtoDUNEAnalysisTree(ProtoDUNEAnalysisTree &&) = delete;
  ProtoDUNEAnalysisTree & operator = (ProtoDUNEAnalysisTree const &) = delete;
  ProtoDUNEAnalysisTree & operator = (ProtoDUNEAnalysisTree &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  void beginRun(const art::Run& run) override;

  // Required functions.
  void analyze(art::Event const & evt) override;

private:

  // Helper utility functions
  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEBeamlineUtils beamlineUtil;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_service;

  // Track momentum algorithm calculates momentum based on track range
  trkf::TrackMomentumCalculator trmom;

  // Geometry
  geo::GeometryCore const * fGeometry = art::ServiceHandle<geo::Geometry>().get();
  geo::WireReadoutGeom const* fWireReadoutGeom = &art::ServiceHandle<geo::WireReadout>()->Get();

  // Initialise tree variables
  void Initialise();

  // Fill configuration tree
  void FillConfigTree();

  bool FillPrimaryBeamParticle(art::Event const & evt);
  void FillPrimaryPFParticle(art::Event const & evt,
                             detinfo::DetectorClocksData const& clockData,
                             detinfo::DetectorPropertiesData const& detProp,
                             const recob::PFParticle* particle);
  void FillPrimaryDaughterPFParticle(art::Event const & evt,
                                     detinfo::DetectorClocksData const& clockData,
                                     detinfo::DetectorPropertiesData const& detProp,
                                     const recob::PFParticle* daughterParticle,
                                     int daughterID);
  void FillPrimaryGrandDaughterPFParticle(art::Event const & evt,
                                          detinfo::DetectorClocksData const& clockData,
                                          detinfo::DetectorPropertiesData const& detProp,
                                          const recob::PFParticle* gdaughterParticle,
                                          int daughterID, int gdaughterID);

  // fcl parameters
  const art::InputTag fBeamModuleLabel;
  std::string fCalorimetryTag;
  std::string fParticleIDTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fSimulationTag;
  int fVerbose;

  // Define configuration tree
  TTree *fConfigTree;

  int fNAPAs;
  int fNChansPerAPA;
  int fNCryostats;
  int fNTPCs;
  int fNChannels;
  int fNPlanes;
  double fActiveTPCBoundsX[2];
  double fActiveTPCBoundsY[2];
  double fActiveTPCBoundsZ[2];

  // Beam track tree
  TTree *fBeamTrack;

  // Primary reco track/showers tree
  TTree *fPandoraBeam;

  // Reco daughters tree
  TTree* fDaughterTree;

  // Reco granddaughters tree
  TTree* fGrandDaughterTree;

  int fRun;
  int fSubRun;
  int fevent;
  double fTimeStamp;
  int fNactivefembs[6];

  // Beam track
  int fbeamNTracks;
  int fbeamtrigger;
  double ftof;
  double ftof_expElec;
  double ftof_expMuon;
  double ftof_expPion;
  double ftof_expKaon;
  double ftof_expProt;
  double ftof_expDeut;
  int fcerenkovStatus[2];
  double fcerenkovTime[2];
  double fcerenkovPressure[2];
  double fbeamtrackMomentum;
  double fbeamtrackEnergy;
  double fbeamtrackPos[4];
  double fbeamtrackDir[3];
  int fbeamtrackPdg;
  int fbeamtrackID;
  int fBIAndTimingMatched;

  // Check what the beam track will do in the detector. The beam track couldn't be the same with the primary beam PFParticle.
  int fbeamtrack_truth_Origin;
  double fbeamtrack_truth_EndPos[4];
  double fbeamtrack_truth_Momentum[4];
  double fbeamtrack_truth_KinEnergy;
  double fbeamtrack_truth_Pt;
  double fbeamtrack_truth_Mass;
  double fbeamtrack_truth_Theta;
  double fbeamtrack_truth_Phi;
  double fbeamtrack_truth_TotalLength;
  int fbeamtrack_truth_Process;
  int fbeamtrack_truth_EndProcess;
  int fbeamtrack_truth_Pdg;

  double fbeamtrack_truth_Pos_InTPCActive[4];
  double fbeamtrack_truth_Momentum_InTPCActive[4];
  double fbeamtrack_truth_P_InTPCActive;
  double fbeamtrack_truth_Pt_InTPCActive;
  double fbeamtrack_truth_Theta_InTPCActive;
  double fbeamtrack_truth_Phi_InTPCActive;
  double fbeamtrack_truth_TotalLength_InTPCActive;
  double fbeamtrack_truth_KinEnergy_InTPCActive;

  int fbeamtrack_truth_NDAUGTHERS;
  int fbeamtrack_truthdaughter_TrackId[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdaughter_Pdg[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdaughter_Mother[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdaughter_StartPosition[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdaughter_EndPosition[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdaughter_P[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdaughter_Momentum[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdaughter_EndMomentum[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdaughter_Pt[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdaughter_Mass[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdaughter_Theta[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdaughter_Phi[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdaughter_TotalLength[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdaughter_Process[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdaughter_EndProcess[NMAXTRUTHDAUGTHERS];

  int fbeamtrack_truth_NDECAYDAUGTHERS;
  int fbeamtrack_truthdecaydaughter_TrackId[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdecaydaughter_Pdg[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdecaydaughter_Mother[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdecaydaughter_StartPosition[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdecaydaughter_EndPosition[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdecaydaughter_P[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdecaydaughter_Momentum[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdecaydaughter_EndMomentum[NMAXTRUTHDAUGTHERS][4];
  double fbeamtrack_truthdecaydaughter_Pt[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdecaydaughter_Mass[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdecaydaughter_Theta[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdecaydaughter_Phi[NMAXTRUTHDAUGTHERS];
  double fbeamtrack_truthdecaydaughter_TotalLength[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdecaydaughter_Process[NMAXTRUTHDAUGTHERS];
  int fbeamtrack_truthdecaydaughter_EndProcess[NMAXTRUTHDAUGTHERS];

  // Reconstructed tracks/showers
  int fNpfParticles;
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
  double fprimaryShowerEnergy[3];
  double fprimaryShowerMIPEnergy[3];
  double fprimaryShowerdEdx[3];
  double fprimaryMomentumByRangeProton;
  double fprimaryMomentumByRangeMuon;
  double fprimaryKineticEnergy[3];
  double fprimaryRange[3];
  double fprimaryTrkPitchC[3];
  int fprimaryID;
  double fprimaryT0;
  double fprimaryEarliestHitPeakTime;

  int fprimaryPID_Pdg[3];
  int fprimaryPID_Ndf[3];
  double fprimaryPID_MinChi2[3];
  double fprimaryPID_DeltaChi2[3];
  double fprimaryPID_Chi2Proton[3];
  double fprimaryPID_Chi2Kaon[3];
  double fprimaryPID_Chi2Pion[3];
  double fprimaryPID_Chi2Muon[3];
  double fprimaryPID_MissingE[3];
  double fprimaryPID_MissingEavg[3];
  double fprimaryPID_PIDA[3];

  int fprimary_truth_Origin;
  int fprimary_truth_TrackId;
  int fprimary_truth_Pdg;
  int fprimary_truth_Mother;
  double fprimary_truth_StartPosition[4];
  double fprimary_truth_EndPosition[4];
  double fprimary_truth_P;
  double fprimary_truth_Momentum[4];
  double fprimary_truth_EndMomentum[4];
  double fprimary_truth_Pt;
  double fprimary_truth_Mass;
  double fprimary_truth_Theta;
  double fprimary_truth_Phi;
  double fprimary_truth_TotalLength;
  int fprimary_truth_Process;
  int fprimary_truth_EndProcess;
  int fprimary_truth_Isbeammatched;
  double fprimary_truth_EkinAtVertex;
  double fprimary_truth_EkinAtVertex_notcorrected;
  double fprimary_truth_Pos_InTPCActive[4];
  double fprimary_truth_Momentum_InTPCActive[4];
  double fprimary_truth_P_InTPCActive;
  double fprimary_truth_Pt_InTPCActive;
  double fprimary_truth_Theta_InTPCActive;
  double fprimary_truth_Phi_InTPCActive;
  double fprimary_truth_TotalLength_InTPCActive;
  double fprimary_truth_KinEnergy_InTPCActive;

  int fprimary_truth_NDAUGTHERS;
  int fprimary_truthdaughter_TrackId[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdaughter_Pdg[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdaughter_Mother[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdaughter_StartPosition[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdaughter_EndPosition[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdaughter_P[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdaughter_Momentum[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdaughter_EndMomentum[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdaughter_Pt[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdaughter_Mass[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdaughter_Theta[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdaughter_Phi[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdaughter_TotalLength[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdaughter_Process[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdaughter_EndProcess[NMAXTRUTHDAUGTHERS];

  int fprimary_truth_NDECAYDAUGTHERS;
  int fprimary_truthdecaydaughter_TrackId[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdecaydaughter_Pdg[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdecaydaughter_Mother[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdecaydaughter_StartPosition[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdecaydaughter_EndPosition[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdecaydaughter_P[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdecaydaughter_Momentum[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdecaydaughter_EndMomentum[NMAXTRUTHDAUGTHERS][4];
  double fprimary_truthdecaydaughter_Pt[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdecaydaughter_Mass[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdecaydaughter_Theta[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdecaydaughter_Phi[NMAXTRUTHDAUGTHERS];
  double fprimary_truthdecaydaughter_TotalLength[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdecaydaughter_Process[NMAXTRUTHDAUGTHERS];
  int fprimary_truthdecaydaughter_EndProcess[NMAXTRUTHDAUGTHERS];

  // Daughters from primary
  int fNDAUGHTERS;
  double fdaughterVertex[3];
  int fdaughterIstrack[NMAXDAUGTHERS];
  int fdaughterIsshower[NMAXDAUGTHERS];
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
  double fdaughterShowerEnergy[NMAXDAUGTHERS][3];
  double fdaughterShowerMIPEnergy[NMAXDAUGTHERS][3];
  double fdaughterShowerdEdx[NMAXDAUGTHERS][3];
  int fdaughterShowerBestPlane[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeProton[NMAXDAUGTHERS];
  double fdaughterMomentumByRangeMuon[NMAXDAUGTHERS];
  double fdaughterKineticEnergy[NMAXDAUGTHERS][3];
  double fdaughterRange[NMAXDAUGTHERS][3];
  double fdaughterTrkPitchC[NMAXDAUGTHERS][3];
  int fdaughterID[NMAXDAUGTHERS];
  double fdaughterT0[NMAXDAUGTHERS];

  int fdaughterPID_Pdg[NMAXDAUGTHERS][3];
  int fdaughterPID_Ndf[NMAXDAUGTHERS][3];
  double fdaughterPID_MinChi2[NMAXDAUGTHERS][3];
  double fdaughterPID_DeltaChi2[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Proton[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Kaon[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Pion[NMAXDAUGTHERS][3];
  double fdaughterPID_Chi2Muon[NMAXDAUGTHERS][3];
  double fdaughterPID_MissingE[NMAXDAUGTHERS][3];
  double fdaughterPID_MissingEavg[NMAXDAUGTHERS][3];
  double fdaughterPID_PIDA[NMAXDAUGTHERS][3];

  int fdaughter_truth_TrackId[NMAXDAUGTHERS];
  int fdaughter_truth_Pdg[NMAXDAUGTHERS];
  int fdaughter_truth_Mother[NMAXDAUGTHERS];
  double fdaughter_truth_StartPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndPosition[NMAXDAUGTHERS][4];
  double fdaughter_truth_P[NMAXDAUGTHERS];
  double fdaughter_truth_Momentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_EndMomentum[NMAXDAUGTHERS][4];
  double fdaughter_truth_Pt[NMAXDAUGTHERS];
  double fdaughter_truth_Mass[NMAXDAUGTHERS];
  double fdaughter_truth_Theta[NMAXDAUGTHERS];
  double fdaughter_truth_Phi[NMAXDAUGTHERS];
  double fdaughter_truth_TotalLength[NMAXDAUGTHERS];
  int fdaughter_truth_Process[NMAXDAUGTHERS];
  int fdaughter_truth_EndProcess[NMAXDAUGTHERS];

  int fNGRANDDAUGHTERS;
  double fgranddaughterVertex[NMAXDAUGTHERS][3];
  int fgranddaughterIstrack[NMAXDAUGTHERS];
  int fgranddaughterIsshower[NMAXDAUGTHERS];
  int fgranddaughterNHits[NMAXDAUGTHERS];
  double fgranddaughterTheta[NMAXDAUGTHERS];
  double fgranddaughterPhi[NMAXDAUGTHERS];
  double fgranddaughterLength[NMAXDAUGTHERS];
  double fgranddaughterMomentum[NMAXDAUGTHERS];
  double fgranddaughterEndMomentum[NMAXDAUGTHERS];
  double fgranddaughterEndPosition[NMAXDAUGTHERS][3];
  double fgranddaughterStartPosition[NMAXDAUGTHERS][3];
  double fgranddaughterEndDirection[NMAXDAUGTHERS][3];
  double fgranddaughterStartDirection[NMAXDAUGTHERS][3];
  double fgranddaughterOpeningAngle[NMAXDAUGTHERS];
  double fgranddaughterShowerEnergy[NMAXDAUGTHERS][3];
  double fgranddaughterShowerMIPEnergy[NMAXDAUGTHERS][3];
  double fgranddaughterShowerdEdx[NMAXDAUGTHERS][3];
  int fgranddaughterShowerBestPlane[NMAXDAUGTHERS];
  double fgranddaughterMomentumByRangeProton[NMAXDAUGTHERS];
  double fgranddaughterMomentumByRangeMuon[NMAXDAUGTHERS];
  double fgranddaughterKineticEnergy[NMAXDAUGTHERS][3];
  double fgranddaughterRange[NMAXDAUGTHERS][3];
  double fgranddaughterTrkPitchC[NMAXDAUGTHERS][3];
  int fgranddaughterID[NMAXDAUGTHERS];
  int fgranddaughterMotherID[NMAXDAUGTHERS];
  double fgranddaughterT0[NMAXDAUGTHERS];

  int fgranddaughterPID_Pdg[NMAXDAUGTHERS][3];
  int fgranddaughterPID_Ndf[NMAXDAUGTHERS][3];
  double fgranddaughterPID_MinChi2[NMAXDAUGTHERS][3];
  double fgranddaughterPID_DeltaChi2[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Proton[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Kaon[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Pion[NMAXDAUGTHERS][3];
  double fgranddaughterPID_Chi2Muon[NMAXDAUGTHERS][3];
  double fgranddaughterPID_MissingE[NMAXDAUGTHERS][3];
  double fgranddaughterPID_MissingEavg[NMAXDAUGTHERS][3];
  double fgranddaughterPID_PIDA[NMAXDAUGTHERS][3];

  int fgranddaughter_truth_TrackId[NMAXDAUGTHERS];
  int fgranddaughter_truth_Pdg[NMAXDAUGTHERS];
  int fgranddaughter_truth_Mother[NMAXDAUGTHERS];
  double fgranddaughter_truth_StartPosition[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_EndPosition[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_P[NMAXDAUGTHERS];
  double fgranddaughter_truth_Momentum[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_EndMomentum[NMAXDAUGTHERS][4];
  double fgranddaughter_truth_Pt[NMAXDAUGTHERS];
  double fgranddaughter_truth_Mass[NMAXDAUGTHERS];
  double fgranddaughter_truth_Theta[NMAXDAUGTHERS];
  double fgranddaughter_truth_Phi[NMAXDAUGTHERS];
  double fgranddaughter_truth_TotalLength[NMAXDAUGTHERS];
  int fgranddaughter_truth_Process[NMAXDAUGTHERS];
  int fgranddaughter_truth_EndProcess[NMAXDAUGTHERS];

};

// -----------------------------------------------------------------------------
protoana::ProtoDUNEAnalysisTree::ProtoDUNEAnalysisTree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  beamlineUtil(p.get<fhicl::ParameterSet>("BeamLineUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fParticleIDTag(p.get<std::string>("ParticleIDTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fSimulationTag(p.get<std::string>("SimulationTag")),
  fVerbose(p.get<int>("Verbose"))
{

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::beginRun(const art::Run&){

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::analyze(art::Event const & evt){

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

  // Fill beam track info
  bool isbeamTriggerEvent = FillPrimaryBeamParticle(evt);
  // If this is not a beam trigger then do not save
  if(!isbeamTriggerEvent) return;

  /*
  // Now we want to access the output from Pandora. This comes in the form of particle flow objects (recob::PFParticle).
  // The primary PFParticles are those we want to consider and these PFParticles then have a hierarchy of daughters that
  // describe the whole interaction of a given primary particle
  //
  //                     / daughter track    / grand-daughter track
  //                    /                   /
  //  primary track    /                   /
  //  ---------------- ---- daughter track ----- grand-daughter track
  //                   \                   
  //                   /\-
  //                   /\\-- daughter shower
  //
  // The above primary PFParticle will have links to three daughter particles, two track-like and one shower-like
  */

  //trmom.SetMinLength(100);

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  //std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag) << std::endl;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  fNpfParticles = pfParticles.size();

  // We can now look at these particles
  for(const recob::PFParticle* particle : pfParticles){

    FillPrimaryPFParticle(evt, clockData, detProp, particle);

    // Find the particle vertex. We need the tracker tag here because we need to do a bit of
    // additional work if the PFParticle is track-like to find the vertex. 
    const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fprimaryVertex[0] = vtx.X(); fprimaryVertex[1] = vtx.Y(); fprimaryVertex[2] = vtx.Z();

    // Now we can look for the interaction point of the particle if one exists, i.e where the particle
    // scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
    // check this by making sure we get a sensible position
    const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
    fdaughterVertex[0] = interactionVtx.X(); fdaughterVertex[1] = interactionVtx.Y(); fdaughterVertex[2] = interactionVtx.Z();

    // Maximum number of daugthers to be processed
    if(particle->NumDaughters() > NMAXDAUGTHERS)
      std::cout << "INFO::Number of daughters is " << particle->NumDaughters() << ". Only the first NMAXDAUGTHERS are processed." << std::endl;

    // Let's get the daughter PFParticles... we can do this simply without the utility
    for(const int daughterID : particle->Daughters()){
      // Daughter ID is the element of the original recoParticle vector
      const recob::PFParticle *daughterParticle      = &(recoParticles->at(daughterID));
      //std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;

      // Fill tree with daughter info
      FillPrimaryDaughterPFParticle(evt, clockData, detProp, daughterParticle, daughterID);

      // Get the secondary vertex from daughter interactions
      const TVector3 secinteractionVtx = pfpUtil.GetPFParticleSecondaryVertex(*daughterParticle,evt,fPFParticleTag,fTrackerTag);
      fgranddaughterVertex[fNDAUGHTERS][0] = secinteractionVtx.X(); fgranddaughterVertex[fNDAUGHTERS][1] = secinteractionVtx.Y(); fgranddaughterVertex[fNDAUGHTERS][2] = secinteractionVtx.Z();

      for(const int gdaughterID : daughterParticle->Daughters()){
        FillPrimaryGrandDaughterPFParticle(evt, clockData, detProp, daughterParticle, gdaughterID, daughterID);
	
	// Only process NMAXDAUGTHERS
	if(fNGRANDDAUGHTERS > NMAXDAUGTHERS) break;
      }

      // Only process NMAXDAUGTHERS
      if(fNDAUGHTERS > NMAXDAUGTHERS) break;
      
    }

    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  } 

  // Fill trees
  fPandoraBeam->Fill();
  fBeamTrack->Fill();
  fDaughterTree->Fill();
  fGrandDaughterTree->Fill();

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::endJob(){

}

// -----------------------------------------------------------------------------
bool protoana::ProtoDUNEAnalysisTree::FillPrimaryBeamParticle(art::Event const & evt){
  
  bool beamTriggerEvent = false;
  // If this event is MC then we can check what the true beam particle is
  if(!evt.isRealData()){
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    // simulation has fGeneratorTag = "generator"
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    fbeamtrack_truth_Origin = (*mcTruths)[0].Origin();
    fbeamNTracks = 1;

    // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
    // of these, so we pass the first element into the function to get the good particle
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);

    if(geantGoodParticle != 0x0){

      beamTriggerEvent             = true;
      fbeamtrigger                 = 12;
      fbeamtrackPos[0]             = geantGoodParticle->Vx();
      fbeamtrackPos[1]             = geantGoodParticle->Vy();
      fbeamtrackPos[2]             = geantGoodParticle->Vz();
      fbeamtrackPos[3]             = geantGoodParticle->T();
      fbeamtrackMomentum           = geantGoodParticle->P();
      fbeamtrackEnergy             = geantGoodParticle->E();
      fbeamtrackPdg                = geantGoodParticle->PdgCode();
      fbeamtrackID                 = geantGoodParticle->TrackId();

      fbeamtrack_truth_EndPos[0]   = geantGoodParticle->EndX();
      fbeamtrack_truth_EndPos[1]   = geantGoodParticle->EndY();
      fbeamtrack_truth_EndPos[2]   = geantGoodParticle->EndZ();
      fbeamtrack_truth_EndPos[3]   = geantGoodParticle->EndT();
      fbeamtrack_truth_Momentum[0] = geantGoodParticle->Px();
      fbeamtrack_truth_Momentum[1] = geantGoodParticle->Py();
      fbeamtrack_truth_Momentum[2] = geantGoodParticle->Pz();
      fbeamtrack_truth_Momentum[3] = geantGoodParticle->E();
      fbeamtrack_truth_Pt          = geantGoodParticle->Pt();
      fbeamtrack_truth_Mass        = geantGoodParticle->Mass();
      fbeamtrack_truth_Theta       = geantGoodParticle->Momentum().Theta();
      fbeamtrack_truth_Phi         = geantGoodParticle->Momentum().Phi();
      fbeamtrack_truth_TotalLength = geantGoodParticle->Trajectory().TotalLength();
      fbeamtrack_truth_KinEnergy   = std::sqrt(geantGoodParticle->P()*geantGoodParticle->P() + geantGoodParticle->Mass()*geantGoodParticle->Mass()) - geantGoodParticle->Mass();
      //fbeamtrack_truth_Process     = int(geantGoodParticle->Trajectory().ProcessToKey(geantGoodParticle->Process()));
      //fbeamtrack_truth_EndProcess  = int(geantGoodParticle->Trajectory().ProcessToKey(geantGoodParticle->EndProcess()));
      fbeamtrack_truth_Process     = truthUtil.GetProcessKey(geantGoodParticle->Process());
      fbeamtrack_truth_EndProcess  = truthUtil.GetProcessKey(geantGoodParticle->EndProcess());
      // Truth Pdg may not be the same as the beam PID.
      fbeamtrack_truth_Pdg         = geantGoodParticle->PdgCode();
      // Get the position and kinematics in the TPC active volume
      int trajintpcactivevol       = truthUtil.GetFirstTrajectoryPointInTPCActiveVolume(*geantGoodParticle, fActiveTPCBoundsX[0], fActiveTPCBoundsX[1], fActiveTPCBoundsY[0], fActiveTPCBoundsY[1], fActiveTPCBoundsZ[0], fActiveTPCBoundsZ[1]);
      if(trajintpcactivevol > 0){
	fbeamtrack_truth_Pos_InTPCActive[0]      = geantGoodParticle->Vx(trajintpcactivevol);
	fbeamtrack_truth_Pos_InTPCActive[1]      = geantGoodParticle->Vy(trajintpcactivevol);
	fbeamtrack_truth_Pos_InTPCActive[2]      = geantGoodParticle->Vz(trajintpcactivevol);
	fbeamtrack_truth_Pos_InTPCActive[3]      = geantGoodParticle->T(trajintpcactivevol);
	fbeamtrack_truth_Momentum_InTPCActive[0] = geantGoodParticle->Px(trajintpcactivevol);
	fbeamtrack_truth_Momentum_InTPCActive[1] = geantGoodParticle->Py(trajintpcactivevol);
	fbeamtrack_truth_Momentum_InTPCActive[2] = geantGoodParticle->Pz(trajintpcactivevol);
	fbeamtrack_truth_Momentum_InTPCActive[3] = geantGoodParticle->E(trajintpcactivevol);
	fbeamtrack_truth_P_InTPCActive           = geantGoodParticle->P(trajintpcactivevol);
	fbeamtrack_truth_Pt_InTPCActive          = geantGoodParticle->Pt(trajintpcactivevol);
	fbeamtrack_truth_Theta_InTPCActive       = geantGoodParticle->Momentum(trajintpcactivevol).Theta();
	fbeamtrack_truth_Phi_InTPCActive         = geantGoodParticle->Momentum(trajintpcactivevol).Phi();
	fbeamtrack_truth_TotalLength_InTPCActive = truthUtil.GetMCParticleLengthInTPCActiveVolume(*geantGoodParticle, fActiveTPCBoundsX[0], fActiveTPCBoundsX[1], fActiveTPCBoundsY[0], fActiveTPCBoundsY[1], fActiveTPCBoundsZ[0], fActiveTPCBoundsZ[1]);
	fbeamtrack_truth_KinEnergy_InTPCActive   = std::sqrt(geantGoodParticle->P(trajintpcactivevol)*geantGoodParticle->P(trajintpcactivevol) + geantGoodParticle->Mass()*geantGoodParticle->Mass()) - geantGoodParticle->Mass();
	//unsigned int ntrajpoints = geantGoodParticle->NumberTrajectoryPoints();
      }

      // If this beam particle has interacted in the detector, find its daugthers.
      std::vector<art::Ptr<simb::MCParticle> > mcbeampartilelist;
      auto mcbeamHandle = evt.getHandle<std::vector<simb::MCParticle> >(fSimulationTag);
      if (mcbeamHandle) art::fill_ptr_vector(mcbeampartilelist, mcbeamHandle);

      for(size_t p=0; p < mcbeampartilelist.size(); p++){
	auto & beamprim_mcparticle =  mcbeampartilelist[p];
	if(beamprim_mcparticle->Mother() == geantGoodParticle->TrackId()){
	  fbeamtrack_truthdaughter_TrackId[fbeamtrack_truth_NDAUGTHERS]            = beamprim_mcparticle->TrackId();
	  fbeamtrack_truthdaughter_Pdg[fbeamtrack_truth_NDAUGTHERS]                = beamprim_mcparticle->PdgCode();
	  fbeamtrack_truthdaughter_Mother[fbeamtrack_truth_NDAUGTHERS]             = beamprim_mcparticle->Mother();
	  fbeamtrack_truthdaughter_StartPosition[fbeamtrack_truth_NDAUGTHERS][0]   = beamprim_mcparticle->Vx();
	  fbeamtrack_truthdaughter_StartPosition[fbeamtrack_truth_NDAUGTHERS][1]   = beamprim_mcparticle->Vy();
	  fbeamtrack_truthdaughter_StartPosition[fbeamtrack_truth_NDAUGTHERS][2]   = beamprim_mcparticle->Vz();
	  fbeamtrack_truthdaughter_StartPosition[fbeamtrack_truth_NDAUGTHERS][3]   = beamprim_mcparticle->T();
	  fbeamtrack_truthdaughter_EndPosition[fbeamtrack_truth_NDAUGTHERS][0]     = beamprim_mcparticle->EndX();
	  fbeamtrack_truthdaughter_EndPosition[fbeamtrack_truth_NDAUGTHERS][1]     = beamprim_mcparticle->EndY();
	  fbeamtrack_truthdaughter_EndPosition[fbeamtrack_truth_NDAUGTHERS][2]     = beamprim_mcparticle->EndZ();
	  fbeamtrack_truthdaughter_EndPosition[fbeamtrack_truth_NDAUGTHERS][3]     = beamprim_mcparticle->EndT();
	  fbeamtrack_truthdaughter_P[fbeamtrack_truth_NDAUGTHERS]                  = beamprim_mcparticle->P();
	  fbeamtrack_truthdaughter_Momentum[fbeamtrack_truth_NDAUGTHERS][0]        = beamprim_mcparticle->Px();
	  fbeamtrack_truthdaughter_Momentum[fbeamtrack_truth_NDAUGTHERS][1]        = beamprim_mcparticle->Py();
	  fbeamtrack_truthdaughter_Momentum[fbeamtrack_truth_NDAUGTHERS][2]        = beamprim_mcparticle->Pz();
	  fbeamtrack_truthdaughter_Momentum[fbeamtrack_truth_NDAUGTHERS][3]        = beamprim_mcparticle->E();
	  fbeamtrack_truthdaughter_Pt[fbeamtrack_truth_NDAUGTHERS]                 = beamprim_mcparticle->Pt();
	  fbeamtrack_truthdaughter_Mass[fbeamtrack_truth_NDAUGTHERS]               = beamprim_mcparticle->Mass();
	  fbeamtrack_truthdaughter_EndMomentum[fbeamtrack_truth_NDAUGTHERS][0]     = beamprim_mcparticle->EndPx();
	  fbeamtrack_truthdaughter_EndMomentum[fbeamtrack_truth_NDAUGTHERS][1]     = beamprim_mcparticle->EndPy();
	  fbeamtrack_truthdaughter_EndMomentum[fbeamtrack_truth_NDAUGTHERS][2]     = beamprim_mcparticle->EndPz();
	  fbeamtrack_truthdaughter_EndMomentum[fbeamtrack_truth_NDAUGTHERS][3]     = beamprim_mcparticle->EndE();
	  fbeamtrack_truthdaughter_Theta[fbeamtrack_truth_NDAUGTHERS]              = beamprim_mcparticle->Momentum().Theta();
	  fbeamtrack_truthdaughter_Phi[fbeamtrack_truth_NDAUGTHERS]                = beamprim_mcparticle->Momentum().Phi();
	  fbeamtrack_truthdaughter_TotalLength[fbeamtrack_truth_NDAUGTHERS]        = beamprim_mcparticle->Trajectory().TotalLength();
	  fbeamtrack_truthdaughter_Process[fbeamtrack_truth_NDAUGTHERS]            = truthUtil.GetProcessKey(beamprim_mcparticle->Process());
	  fbeamtrack_truthdaughter_EndProcess[fbeamtrack_truth_NDAUGTHERS]         = truthUtil.GetProcessKey(beamprim_mcparticle->EndProcess());

	  fbeamtrack_truth_NDAUGTHERS++;
	  if(fbeamtrack_truth_NDAUGTHERS > NMAXTRUTHDAUGTHERS) break;

	  if(truthUtil.GetProcessKey(beamprim_mcparticle->EndProcess()) == 30){ // Decay particle
	    for(size_t pp=0; pp < mcbeampartilelist.size(); pp++){
	      auto & decay_beamprim_mcparticle =  mcbeampartilelist[pp];
	      if(decay_beamprim_mcparticle->Mother() == beamprim_mcparticle->TrackId()){
		fbeamtrack_truthdecaydaughter_TrackId[fbeamtrack_truth_NDECAYDAUGTHERS]            = decay_beamprim_mcparticle->TrackId();
		fbeamtrack_truthdecaydaughter_Pdg[fbeamtrack_truth_NDECAYDAUGTHERS]                = decay_beamprim_mcparticle->PdgCode();
		fbeamtrack_truthdecaydaughter_Mother[fbeamtrack_truth_NDECAYDAUGTHERS]             = decay_beamprim_mcparticle->Mother();
		fbeamtrack_truthdecaydaughter_StartPosition[fbeamtrack_truth_NDECAYDAUGTHERS][0]   = decay_beamprim_mcparticle->Vx();
		fbeamtrack_truthdecaydaughter_StartPosition[fbeamtrack_truth_NDECAYDAUGTHERS][1]   = decay_beamprim_mcparticle->Vy();
		fbeamtrack_truthdecaydaughter_StartPosition[fbeamtrack_truth_NDECAYDAUGTHERS][2]   = decay_beamprim_mcparticle->Vz();
		fbeamtrack_truthdecaydaughter_StartPosition[fbeamtrack_truth_NDECAYDAUGTHERS][3]   = decay_beamprim_mcparticle->T();
		fbeamtrack_truthdecaydaughter_EndPosition[fbeamtrack_truth_NDECAYDAUGTHERS][0]     = decay_beamprim_mcparticle->EndX();
		fbeamtrack_truthdecaydaughter_EndPosition[fbeamtrack_truth_NDECAYDAUGTHERS][1]     = decay_beamprim_mcparticle->EndY();
		fbeamtrack_truthdecaydaughter_EndPosition[fbeamtrack_truth_NDECAYDAUGTHERS][2]     = decay_beamprim_mcparticle->EndZ();
		fbeamtrack_truthdecaydaughter_EndPosition[fbeamtrack_truth_NDECAYDAUGTHERS][3]     = decay_beamprim_mcparticle->EndT();
		fbeamtrack_truthdecaydaughter_P[fbeamtrack_truth_NDECAYDAUGTHERS]                  = decay_beamprim_mcparticle->P();
		fbeamtrack_truthdecaydaughter_Momentum[fbeamtrack_truth_NDECAYDAUGTHERS][0]        = decay_beamprim_mcparticle->Px();
		fbeamtrack_truthdecaydaughter_Momentum[fbeamtrack_truth_NDECAYDAUGTHERS][1]        = decay_beamprim_mcparticle->Py();
		fbeamtrack_truthdecaydaughter_Momentum[fbeamtrack_truth_NDECAYDAUGTHERS][2]        = decay_beamprim_mcparticle->Pz();
		fbeamtrack_truthdecaydaughter_Momentum[fbeamtrack_truth_NDECAYDAUGTHERS][3]        = decay_beamprim_mcparticle->E();
		fbeamtrack_truthdecaydaughter_Pt[fbeamtrack_truth_NDECAYDAUGTHERS]                 = decay_beamprim_mcparticle->Pt();
		fbeamtrack_truthdecaydaughter_Mass[fbeamtrack_truth_NDECAYDAUGTHERS]               = decay_beamprim_mcparticle->Mass();
		fbeamtrack_truthdecaydaughter_EndMomentum[fbeamtrack_truth_NDECAYDAUGTHERS][0]     = decay_beamprim_mcparticle->EndPx();
		fbeamtrack_truthdecaydaughter_EndMomentum[fbeamtrack_truth_NDECAYDAUGTHERS][1]     = decay_beamprim_mcparticle->EndPy();
		fbeamtrack_truthdecaydaughter_EndMomentum[fbeamtrack_truth_NDECAYDAUGTHERS][2]     = decay_beamprim_mcparticle->EndPz();
		fbeamtrack_truthdecaydaughter_EndMomentum[fbeamtrack_truth_NDECAYDAUGTHERS][3]     = decay_beamprim_mcparticle->EndE();
		fbeamtrack_truthdecaydaughter_Theta[fbeamtrack_truth_NDECAYDAUGTHERS]              = decay_beamprim_mcparticle->Momentum().Theta();
		fbeamtrack_truthdecaydaughter_Phi[fbeamtrack_truth_NDECAYDAUGTHERS]                = decay_beamprim_mcparticle->Momentum().Phi();
		fbeamtrack_truthdecaydaughter_TotalLength[fbeamtrack_truth_NDECAYDAUGTHERS]        = decay_beamprim_mcparticle->Trajectory().TotalLength();
		fbeamtrack_truthdecaydaughter_Process[fbeamtrack_truth_NDECAYDAUGTHERS]            = truthUtil.GetProcessKey(decay_beamprim_mcparticle->Process());
		fbeamtrack_truthdecaydaughter_EndProcess[fbeamtrack_truth_NDECAYDAUGTHERS]         = truthUtil.GetProcessKey(decay_beamprim_mcparticle->EndProcess());

		fbeamtrack_truth_NDECAYDAUGTHERS++;
		if(fbeamtrack_truth_NDECAYDAUGTHERS > NMAXTRUTHDAUGTHERS) break;
	      }
	    }
	  } // decay particle
	  
	}

      } // mc particle vector

    } // beam geant particle  found
  } // is mc
  else{ // data
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);

    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    auto pdbeamHandle = evt.getHandle< std::vector<beam::ProtoDUNEBeamEvent> >(fBeamModuleLabel);
    if (pdbeamHandle) art::fill_ptr_vector(beaminfo, pdbeamHandle);

    fbeamNTracks = (int)beaminfo.size();
    // For now only consider events with one beam track only
    if(fbeamNTracks != 1) return false;
  
    for(unsigned int i = 0; i < beaminfo.size(); ++i){
      fBIAndTimingMatched = beaminfo[i]->CheckIsMatched();
      fbeamtrigger = beaminfo[i]->GetTimingTrigger();
      fbeamtrackPos[3] = (double)beaminfo[i]->GetRDTimestamp();

      // If ToF is 0-3 there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
      if(beaminfo[i]->GetTOFChan() >= 0)
	ftof =  beaminfo[i]->GetTOF();

      // Get Cerenkov
      if(beaminfo[i]->GetBITrigger() == 1){
	fcerenkovStatus[0]   = beaminfo[i]->GetCKov0Status();
	fcerenkovStatus[1]   = beaminfo[i]->GetCKov1Status();
	fcerenkovTime[0]     = beaminfo[i]->GetCKov0Time();
	fcerenkovTime[1]     = beaminfo[i]->GetCKov1Time();
	fcerenkovPressure[0] = beaminfo[i]->GetCKov0Pressure();
	fcerenkovPressure[1] = beaminfo[i]->GetCKov1Pressure();
      }

      // Beam particle could have more than one tracks - for now take the first one, need to do this properly
      auto & tracks = beaminfo[i]->GetBeamTracks();
      if(!tracks.empty()){
	fbeamtrackPos[0] = tracks[0].End().X();
	fbeamtrackPos[1] = tracks[0].End().Y();
	fbeamtrackPos[2] = tracks[0].End().Z();
	fbeamtrackDir[0] = tracks[0].StartDirection().X();
	fbeamtrackDir[1] = tracks[0].StartDirection().Y();
	fbeamtrackDir[2] = tracks[0].StartDirection().Z();
      }
  
      // Beam momentum
      auto & beammom = beaminfo[i]->GetRecoBeamMomenta();
      if(!beammom.empty()){
	fbeamtrackMomentum = beammom[0];
	ftof_expElec =  beamlineUtil.ComputeTOF(11,         beammom[0]);
	ftof_expMuon =  beamlineUtil.ComputeTOF(13,         beammom[0]);
	ftof_expPion =  beamlineUtil.ComputeTOF(211,        beammom[0]);
	ftof_expKaon =  beamlineUtil.ComputeTOF(321,        beammom[0]);
	ftof_expProt =  beamlineUtil.ComputeTOF(2212,       beammom[0]);
	ftof_expDeut =  beamlineUtil.ComputeTOF(1000010020, beammom[0]);
      }

      // For now only take the first beam particle - need to add some criteria if more than one are found
      break;
 
    }
  }

  return beamTriggerEvent;

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::FillPrimaryPFParticle(art::Event const & evt,
                                                            detinfo::DetectorClocksData const& clockData,
                                                            detinfo::DetectorPropertiesData const& detProp,
                                                            const recob::PFParticle* particle){

  // Pandora's BDT beam-cosmic score
  fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);
  
  // NHits associated with this pfParticle
  fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();

  // Earliest hit peak time
  fprimaryEarliestHitPeakTime = pfpUtil.GetPFParticleEarliestHitPeakTime(*particle,evt,fPFParticleTag);
  
  // Get the T0 for this pfParticle
  std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
  if(!pfT0vec.empty())
    fprimaryT0 = pfT0vec[0].Time();

  //std::vector<art::Ptr<simb::MCParticle> > mcpartilelist;
  //mf::LogWarning("ProtoDUNEAnalysisTree") << " Getting MC Particle failed. Ignore if you run on data. " << std::endl;
  //auto mcHandle = evt.getHandle<std::vector<simb::MCParticle> >(fSimulationTag);
  //if (mcHandle)
  //art::fill_ptr_vector(mcpartilelist, mcHandle);

  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
  const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle, evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);

  const simb::MCParticle* mcparticle = NULL; 

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
    
    fprimaryMomentumByRangeMuon        = trmom.GetTrackMomentum(thisTrack->Length(),13);
    fprimaryMomentumByRangeProton      = trmom.GetTrackMomentum(thisTrack->Length(),2212);      
    
    // Calorimetry
    std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
    if(calovector.size() != 3 && fVerbose > 0)
      std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;
    
    for(size_t k = 0; k < calovector.size() && k<3; k++){
      int plane = calovector[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fprimaryKineticEnergy[plane]      = calovector[k].KineticEnergy();
      fprimaryRange[plane]              = calovector[k].Range();
      fprimaryTrkPitchC[plane]          = calovector[k].TrkPitchC();
    }

    // PID
    std::vector<anab::ParticleID> pids = trackUtil.GetRecoTrackPID(*thisTrack, evt, fTrackerTag, fParticleIDTag);

    // Set dummy values
    double pidpdg[3] = {-1,-1,-1};
    double pidchi[3] = {99999.,99999.,99999.};
    double pid2ndchi[3] = {99999.,99999.,99999.};
    int pidndf[3] = {-1,-1,-1};

    for (size_t ipid=0; ipid<pids.size(); ipid++){ 
      std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid].ParticleIDAlgScores();

      // Loop though AlgScoresVec and find the variables we want

      for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
        anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
        /*std::cout << "\n ParticleIDAlg " << AlgScore.fAlgName
          << "\n -- Variable type: " << AlgScore.fVariableType
          << "\n -- Track direction: " << AlgScore.fTrackDir
          << "\n -- Assuming PDG: " << AlgScore.fAssumedPdg
          << "\n -- Number of degrees of freedom: " << AlgScore.fNdf
          << "\n -- Value: " << AlgScore.fValue
          << "\n -- Using planeMask: " << AlgScore.fPlaneMask << std::endl;*/
        if (AlgScore.fPlaneMask.none() || AlgScore.fPlaneMask.count() > 1) continue;

        int planenum = -1;
        if (AlgScore.fPlaneMask.test(0)) planenum = 0;
        if (AlgScore.fPlaneMask.test(1)) planenum = 1;
        if (AlgScore.fPlaneMask.test(2)) planenum = 2;

        if (planenum<0 || planenum>2) continue;
        if (AlgScore.fAlgName == "Chi2"){
          if (TMath::Abs(AlgScore.fAssumedPdg) == 13){ // chi2mu
            fprimaryPID_Chi2Muon[planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 2212){ // chi2pr
            fprimaryPID_Chi2Proton[planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 211){ // chi2pi
            fprimaryPID_Chi2Pion[planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 321){ // chi2ka
            fprimaryPID_Chi2Kaon[planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
        }
        else if (AlgScore.fVariableType==anab::kPIDA){
          fprimaryPID_PIDA[planenum] = AlgScore.fValue;
        }
      } // end loop though AlgScoresVec
    } // end loop over pid[ipid]

      // Finally, set min chi2
    for (size_t plane=0; plane<3; plane++){
      fprimaryPID_MinChi2[plane] = pidchi[plane];
      fprimaryPID_Pdg[plane] = pidpdg[plane];
      fprimaryPID_Ndf[plane] = pidndf[plane];
      fprimaryPID_DeltaChi2[plane] = pid2ndchi[plane]-pidchi[plane];
    }

    // NOT SET with new PID object
    //fprimaryPID_MissingE[plane]       = pids[plane].MissingE();
    //fprimaryPID_MissingEavg[plane]    = pids[plane].MissingEavg();


    // Get the true mc particle
    //const simb::MCParticle* mcparticle = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackerTag);
    mcparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, evt, fTrackerTag);
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
    for(size_t k = 0; k < (thisShower->Energy()).size() && k<3; k++)
      fprimaryShowerEnergy[k] = thisShower->Energy()[k];
    for(size_t k = 0; k < (thisShower->MIPEnergy()).size() && k<3; k++)
      fprimaryShowerMIPEnergy[k] = thisShower->MIPEnergy()[k];
    for(size_t k = 0; k < (thisShower->dEdx()).size() && k<3; k++)
      fprimaryShowerdEdx[k] = thisShower->dEdx()[k];

    //const simb::MCParticle* mcparticle = truthUtil.GetMCParticleFromRecoShower(*thisShower, evt, fShowerTag);
    mcparticle = truthUtil.GetMCParticleFromRecoShower(clockData, *thisShower, evt, fShowerTag);
  } // end is shower
  else{
    if(fVerbose > 0){
      std::cout << "INFO::Primary pfParticle is not track or shower!" << std::endl;
    }
    return;
  }

  if(mcparticle != 0x0){
    fprimary_truth_Origin             = pi_service->TrackIdToMCTruth_P(mcparticle->TrackId())->Origin();
    fprimary_truth_TrackId            = mcparticle->TrackId();
    fprimary_truth_Pdg                = mcparticle->PdgCode();
    fprimary_truth_Mother             = mcparticle->Mother();
    fprimary_truth_StartPosition[0]   = mcparticle->Vx();
    fprimary_truth_StartPosition[1]   = mcparticle->Vy();
    fprimary_truth_StartPosition[2]   = mcparticle->Vz();
    fprimary_truth_StartPosition[3]   = mcparticle->T();
    fprimary_truth_EndPosition[0]     = mcparticle->EndX();
    fprimary_truth_EndPosition[1]     = mcparticle->EndY();
    fprimary_truth_EndPosition[2]     = mcparticle->EndZ();
    fprimary_truth_EndPosition[3]     = mcparticle->EndT();
    fprimary_truth_P                  = mcparticle->P();
    fprimary_truth_Momentum[0]        = mcparticle->Px();
    fprimary_truth_Momentum[1]        = mcparticle->Py();
    fprimary_truth_Momentum[2]        = mcparticle->Pz();
    fprimary_truth_Momentum[3]        = mcparticle->E();
    fprimary_truth_Pt                 = mcparticle->Pt();
    fprimary_truth_Mass               = mcparticle->Mass();
    fprimary_truth_EndMomentum[0]     = mcparticle->EndPx();
    fprimary_truth_EndMomentum[1]     = mcparticle->EndPy();
    fprimary_truth_EndMomentum[2]     = mcparticle->EndPz();
    fprimary_truth_EndMomentum[3]     = mcparticle->EndE();
    fprimary_truth_Theta              = mcparticle->Momentum().Theta();
    fprimary_truth_Phi                = mcparticle->Momentum().Phi();
    //fprimary_truth_NDaughters         = mcparticle->NumberDaughters();
    fprimary_truth_TotalLength        = mcparticle->Trajectory().TotalLength();
    fprimary_truth_Process            = truthUtil.GetProcessKey(mcparticle->Process());
    fprimary_truth_EndProcess         = truthUtil.GetProcessKey(mcparticle->EndProcess());
    fprimary_truth_EkinAtVertex_notcorrected = truthUtil.GetKinEnergyAtVertex(*mcparticle);
    double lastpointEloss = truthUtil.GetDepEnergyAtLastTrajPoint(*mcparticle,fActiveTPCBoundsX[0], fActiveTPCBoundsX[1], fActiveTPCBoundsY[0], fActiveTPCBoundsY[1], fActiveTPCBoundsZ[0], fActiveTPCBoundsZ[1]);
    fprimary_truth_EkinAtVertex = truthUtil.GetKinEnergyAtVertex(*mcparticle,lastpointEloss);
    
    // Get the position and kinematics in the TPC active volume
    int trajintpcactivevol       = truthUtil.GetFirstTrajectoryPointInTPCActiveVolume(*mcparticle, fActiveTPCBoundsX[0], fActiveTPCBoundsX[1], fActiveTPCBoundsY[0], fActiveTPCBoundsY[1], fActiveTPCBoundsZ[0], fActiveTPCBoundsZ[1]);
    if(trajintpcactivevol > 0){
      fprimary_truth_Pos_InTPCActive[0]      = mcparticle->Vx(trajintpcactivevol);
      fprimary_truth_Pos_InTPCActive[1]      = mcparticle->Vy(trajintpcactivevol);
      fprimary_truth_Pos_InTPCActive[2]      = mcparticle->Vz(trajintpcactivevol);
      fprimary_truth_Pos_InTPCActive[3]      = mcparticle->T(trajintpcactivevol);
      fprimary_truth_Momentum_InTPCActive[0] = mcparticle->Px(trajintpcactivevol);
      fprimary_truth_Momentum_InTPCActive[1] = mcparticle->Py(trajintpcactivevol);
      fprimary_truth_Momentum_InTPCActive[2] = mcparticle->Pz(trajintpcactivevol);
      fprimary_truth_Momentum_InTPCActive[3] = mcparticle->E(trajintpcactivevol);
      fprimary_truth_P_InTPCActive           = mcparticle->P(trajintpcactivevol);
      fprimary_truth_Pt_InTPCActive          = mcparticle->Pt(trajintpcactivevol);
      fprimary_truth_Theta_InTPCActive       = mcparticle->Momentum(trajintpcactivevol).Theta();
      fprimary_truth_Phi_InTPCActive         = mcparticle->Momentum(trajintpcactivevol).Phi();
      fprimary_truth_TotalLength_InTPCActive = truthUtil.GetMCParticleLengthInTPCActiveVolume(*mcparticle, fActiveTPCBoundsX[0], fActiveTPCBoundsX[1], fActiveTPCBoundsY[0], fActiveTPCBoundsY[1], fActiveTPCBoundsZ[0], fActiveTPCBoundsZ[1]);
      fprimary_truth_KinEnergy_InTPCActive   = std::sqrt(mcparticle->P(trajintpcactivevol)*mcparticle->P(trajintpcactivevol) + mcparticle->Mass()*mcparticle->Mass()) - mcparticle->Mass();
    }
    
    if(fbeamtrackID != -999 && fprimary_truth_TrackId == fbeamtrackID)
      fprimary_truth_Isbeammatched    = 1;
    else
      fprimary_truth_Isbeammatched    = 0;
    
    // If this beam particle has interacted in the detector, find its daugthers.
    std::vector<art::Ptr<simb::MCParticle> > mcprimarypartilelist;
    auto mcprimaryHandle = evt.getHandle<std::vector<simb::MCParticle> >(fSimulationTag);
    if (mcprimaryHandle) art::fill_ptr_vector(mcprimarypartilelist, mcprimaryHandle);
    
    for(size_t p=0; p < mcprimarypartilelist.size(); p++){
      auto & primary_mcparticle =  mcprimarypartilelist[p];
      
      if(primary_mcparticle->Mother() == mcparticle->TrackId()){
	fprimary_truthdaughter_TrackId[fprimary_truth_NDAUGTHERS]            = primary_mcparticle->TrackId();
	fprimary_truthdaughter_Pdg[fprimary_truth_NDAUGTHERS]                = primary_mcparticle->PdgCode();
	fprimary_truthdaughter_Mother[fprimary_truth_NDAUGTHERS]             = primary_mcparticle->Mother();
	fprimary_truthdaughter_StartPosition[fprimary_truth_NDAUGTHERS][0]   = primary_mcparticle->Vx();
	fprimary_truthdaughter_StartPosition[fprimary_truth_NDAUGTHERS][1]   = primary_mcparticle->Vy();
	fprimary_truthdaughter_StartPosition[fprimary_truth_NDAUGTHERS][2]   = primary_mcparticle->Vz();
	fprimary_truthdaughter_StartPosition[fprimary_truth_NDAUGTHERS][3]   = primary_mcparticle->T();
	fprimary_truthdaughter_EndPosition[fprimary_truth_NDAUGTHERS][0]     = primary_mcparticle->EndX();
	fprimary_truthdaughter_EndPosition[fprimary_truth_NDAUGTHERS][1]     = primary_mcparticle->EndY();
	fprimary_truthdaughter_EndPosition[fprimary_truth_NDAUGTHERS][2]     = primary_mcparticle->EndZ();
	fprimary_truthdaughter_EndPosition[fprimary_truth_NDAUGTHERS][3]     = primary_mcparticle->EndT();
	fprimary_truthdaughter_P[fprimary_truth_NDAUGTHERS]                  = primary_mcparticle->P();
	fprimary_truthdaughter_Momentum[fprimary_truth_NDAUGTHERS][0]        = primary_mcparticle->Px();
	fprimary_truthdaughter_Momentum[fprimary_truth_NDAUGTHERS][1]        = primary_mcparticle->Py();
	fprimary_truthdaughter_Momentum[fprimary_truth_NDAUGTHERS][2]        = primary_mcparticle->Pz();
	fprimary_truthdaughter_Momentum[fprimary_truth_NDAUGTHERS][3]        = primary_mcparticle->E();
	fprimary_truthdaughter_Pt[fprimary_truth_NDAUGTHERS]                 = primary_mcparticle->Pt();
	fprimary_truthdaughter_Mass[fprimary_truth_NDAUGTHERS]               = primary_mcparticle->Mass();
	fprimary_truthdaughter_EndMomentum[fprimary_truth_NDAUGTHERS][0]     = primary_mcparticle->EndPx();
	fprimary_truthdaughter_EndMomentum[fprimary_truth_NDAUGTHERS][1]     = primary_mcparticle->EndPy();
	fprimary_truthdaughter_EndMomentum[fprimary_truth_NDAUGTHERS][2]     = primary_mcparticle->EndPz();
	fprimary_truthdaughter_EndMomentum[fprimary_truth_NDAUGTHERS][3]     = primary_mcparticle->EndE();
	fprimary_truthdaughter_Theta[fprimary_truth_NDAUGTHERS]              = primary_mcparticle->Momentum().Theta();
	fprimary_truthdaughter_Phi[fprimary_truth_NDAUGTHERS]                = primary_mcparticle->Momentum().Phi();
	fprimary_truthdaughter_TotalLength[fprimary_truth_NDAUGTHERS]        = primary_mcparticle->Trajectory().TotalLength();
	fprimary_truthdaughter_Process[fprimary_truth_NDAUGTHERS]            = truthUtil.GetProcessKey(primary_mcparticle->Process());
	fprimary_truthdaughter_EndProcess[fprimary_truth_NDAUGTHERS]         = truthUtil.GetProcessKey(primary_mcparticle->EndProcess());
	
	fprimary_truth_NDAUGTHERS++;
	if(fprimary_truth_NDAUGTHERS > NMAXTRUTHDAUGTHERS) break;
	
	if(truthUtil.GetProcessKey(primary_mcparticle->EndProcess()) == 30){ // Decay particle
	  for(size_t pp=0; pp < mcprimarypartilelist.size(); pp++){
	    auto & decay_primary_mcparticle =  mcprimarypartilelist[pp];
	    
	    if(decay_primary_mcparticle->Mother() == primary_mcparticle->TrackId()){
	      fprimary_truthdecaydaughter_TrackId[fprimary_truth_NDECAYDAUGTHERS]            = decay_primary_mcparticle->TrackId();
	      fprimary_truthdecaydaughter_Pdg[fprimary_truth_NDECAYDAUGTHERS]                = decay_primary_mcparticle->PdgCode();
	      fprimary_truthdecaydaughter_Mother[fprimary_truth_NDECAYDAUGTHERS]             = decay_primary_mcparticle->Mother();
	      fprimary_truthdecaydaughter_StartPosition[fprimary_truth_NDECAYDAUGTHERS][0]   = decay_primary_mcparticle->Vx();
	      fprimary_truthdecaydaughter_StartPosition[fprimary_truth_NDECAYDAUGTHERS][1]   = decay_primary_mcparticle->Vy();
	      fprimary_truthdecaydaughter_StartPosition[fprimary_truth_NDECAYDAUGTHERS][2]   = decay_primary_mcparticle->Vz();
	      fprimary_truthdecaydaughter_StartPosition[fprimary_truth_NDECAYDAUGTHERS][3]   = decay_primary_mcparticle->T();
	      fprimary_truthdecaydaughter_EndPosition[fprimary_truth_NDECAYDAUGTHERS][0]     = decay_primary_mcparticle->EndX();
	      fprimary_truthdecaydaughter_EndPosition[fprimary_truth_NDECAYDAUGTHERS][1]     = decay_primary_mcparticle->EndY();
	      fprimary_truthdecaydaughter_EndPosition[fprimary_truth_NDECAYDAUGTHERS][2]     = decay_primary_mcparticle->EndZ();
	      fprimary_truthdecaydaughter_EndPosition[fprimary_truth_NDECAYDAUGTHERS][3]     = decay_primary_mcparticle->EndT();
	      fprimary_truthdecaydaughter_P[fprimary_truth_NDECAYDAUGTHERS]                  = decay_primary_mcparticle->P();
	      fprimary_truthdecaydaughter_Momentum[fprimary_truth_NDECAYDAUGTHERS][0]        = decay_primary_mcparticle->Px();
	      fprimary_truthdecaydaughter_Momentum[fprimary_truth_NDECAYDAUGTHERS][1]        = decay_primary_mcparticle->Py();
	      fprimary_truthdecaydaughter_Momentum[fprimary_truth_NDECAYDAUGTHERS][2]        = decay_primary_mcparticle->Pz();
	      fprimary_truthdecaydaughter_Momentum[fprimary_truth_NDECAYDAUGTHERS][3]        = decay_primary_mcparticle->E();
	      fprimary_truthdecaydaughter_Pt[fprimary_truth_NDECAYDAUGTHERS]                 = decay_primary_mcparticle->Pt();
	      fprimary_truthdecaydaughter_Mass[fprimary_truth_NDECAYDAUGTHERS]               = decay_primary_mcparticle->Mass();
	      fprimary_truthdecaydaughter_EndMomentum[fprimary_truth_NDECAYDAUGTHERS][0]     = decay_primary_mcparticle->EndPx();
	      fprimary_truthdecaydaughter_EndMomentum[fprimary_truth_NDECAYDAUGTHERS][1]     = decay_primary_mcparticle->EndPy();
	      fprimary_truthdecaydaughter_EndMomentum[fprimary_truth_NDECAYDAUGTHERS][2]     = decay_primary_mcparticle->EndPz();
	      fprimary_truthdecaydaughter_EndMomentum[fprimary_truth_NDECAYDAUGTHERS][3]     = decay_primary_mcparticle->EndE();
	      fprimary_truthdecaydaughter_Theta[fprimary_truth_NDECAYDAUGTHERS]              = decay_primary_mcparticle->Momentum().Theta();
	      fprimary_truthdecaydaughter_Phi[fprimary_truth_NDECAYDAUGTHERS]                = decay_primary_mcparticle->Momentum().Phi();
	      fprimary_truthdecaydaughter_TotalLength[fprimary_truth_NDECAYDAUGTHERS]        = decay_primary_mcparticle->Trajectory().TotalLength();
	      fprimary_truthdecaydaughter_Process[fprimary_truth_NDECAYDAUGTHERS]            = truthUtil.GetProcessKey(decay_primary_mcparticle->Process());
	      fprimary_truthdecaydaughter_EndProcess[fprimary_truth_NDECAYDAUGTHERS]         = truthUtil.GetProcessKey(decay_primary_mcparticle->EndProcess());
	      
	      fprimary_truth_NDECAYDAUGTHERS++;
	      if(fprimary_truth_NDECAYDAUGTHERS > NMAXTRUTHDAUGTHERS) break;
	    }
	  }
	} // decay particle
      }
    }
    
    if(fVerbose > 1){
      std::cout << "INFO::Process = " << (mcparticle->Process()).c_str() << " , End process = " << (mcparticle->EndProcess()).c_str()
		<< " , track ID = " << mcparticle->TrackId()
		<< std::endl;
    }
  }
  else{
    if(fVerbose > 0){
      std::cout << "INFO::No MCParticle for primary found!" << std::endl;
    }
  }

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::FillPrimaryDaughterPFParticle(art::Event const & evt,
                                                                    detinfo::DetectorClocksData const& clockData,
                                                                    detinfo::DetectorPropertiesData const& detProp,
                                                                    const recob::PFParticle* daughterParticle,
                                                                    int daughterID){

  const recob::Track* daughterTrack              = pfpUtil.GetPFParticleTrack(*daughterParticle,evt, fPFParticleTag,fTrackerTag);
  const recob::Shower* daughterShower            = pfpUtil.GetPFParticleShower(*daughterParticle,evt,fPFParticleTag,fShowerTag);
  const simb::MCParticle* mcdaughterparticle     = NULL;
  
  if(daughterTrack != 0x0){
    fdaughterIstrack[fNDAUGHTERS]                = 1;
    fdaughterIsshower[fNDAUGHTERS]               = 0;
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
    
    // Calorimetry
    std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
    if(daughtercalovector.size() != 3 && fVerbose > 0)
      std::cerr << "WARNING::Calorimetry vector size for daughter is = " << daughtercalovector.size() << std::endl;
    
    for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
      int plane = daughtercalovector[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fdaughterKineticEnergy[fNDAUGHTERS][plane] = daughtercalovector[k].KineticEnergy();
      fdaughterRange[fNDAUGHTERS][plane]         = daughtercalovector[k].Range();
      fdaughterTrkPitchC[fNDAUGHTERS][plane]     = daughtercalovector[k].TrkPitchC();
    }

    // PID
    std::vector<anab::ParticleID> daughterpids = trackUtil.GetRecoTrackPID(*daughterTrack, evt, fTrackerTag, fParticleIDTag);

    // Set dummy values
    double pidpdg[3] = {-1,-1,-1};
    double pidchi[3] = {99999.,99999.,99999.};
    double pid2ndchi[3] = {99999.,99999.,99999.};
    int pidndf[3] = {-1,-1,-1};

    for (size_t ipid=0; ipid<daughterpids.size(); ipid++){ 
      std::vector<anab::sParticleIDAlgScores> AlgScoresVec = daughterpids[ipid].ParticleIDAlgScores();

      // Loop though AlgScoresVec and find the variables we want

      for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
        anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
        /*std::cout << "\n ParticleIDAlg " << AlgScore.fAlgName
          << "\n -- Variable type: " << AlgScore.fVariableType
          << "\n -- Track direction: " << AlgScore.fTrackDir
          << "\n -- Assuming PDG: " << AlgScore.fAssumedPdg
          << "\n -- Number of degrees of freedom: " << AlgScore.fNdf
          << "\n -- Value: " << AlgScore.fValue
          << "\n -- Using planeMask: " << AlgScore.fPlaneMask << std::endl;*/
        if (AlgScore.fPlaneMask.none() || AlgScore.fPlaneMask.count() > 1) continue;

        int planenum = -1;
        if (AlgScore.fPlaneMask.test(0)) planenum = 0;
        if (AlgScore.fPlaneMask.test(1)) planenum = 1;
        if (AlgScore.fPlaneMask.test(2)) planenum = 2;

        if (planenum<0 || planenum>2) continue;
        if (AlgScore.fAlgName == "Chi2"){
          if (TMath::Abs(AlgScore.fAssumedPdg) == 13){ // chi2mu
            fdaughterPID_Chi2Muon[fNDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 2212){ // chi2pr
            fdaughterPID_Chi2Proton[fNDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 211){ // chi2pi
            fdaughterPID_Chi2Pion[fNDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 321){ // chi2ka
            fdaughterPID_Chi2Kaon[fNDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
        }
        else if (AlgScore.fVariableType==anab::kPIDA){
          fdaughterPID_PIDA[fNDAUGHTERS][planenum] = AlgScore.fValue;
        }
      } // end loop though AlgScoresVec
    } // end loop over pid[ipid]

      // Finally, set min chi2
    for (size_t plane=0; plane<3; plane++){
      fdaughterPID_MinChi2[fNDAUGHTERS][plane] = pidchi[plane];
      fdaughterPID_Pdg[fNDAUGHTERS][plane] = pidpdg[plane];
      fdaughterPID_Ndf[fNDAUGHTERS][plane] = pidndf[plane];
      fdaughterPID_DeltaChi2[fNDAUGHTERS][plane] = pid2ndchi[plane]-pidchi[plane];
    }

    // NOT SET with new PID object
    //fdaughterPID_MissingE[plane]       = daughterpids[plane].MissingE();
    //fdaughterPID_MissingEavg[plane]    = daughterpids[plane].MissingEavg();


    // Get the true mc particle
    //const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(*daughterTrack, evt, fTrackerTag);
    mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *daughterTrack, evt, fTrackerTag);
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
    for(size_t k = 0; k < (daughterShower->Energy()).size() && k<3; k++)
      fdaughterShowerEnergy[fNDAUGHTERS][k] = daughterShower->Energy()[k];
    for(size_t k = 0; k < (daughterShower->MIPEnergy()).size() && k<3; k++)
      fdaughterShowerMIPEnergy[fNDAUGHTERS][k] = daughterShower->MIPEnergy()[k];
    for(size_t k = 0; k < (daughterShower->dEdx()).size() && k<3; k++)
      fdaughterShowerdEdx[fNDAUGHTERS][k] = daughterShower->dEdx()[k];

    // Get the true mc particle
    //const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoShower(*daughterShower, evt, fShowerTag);
    mcdaughterparticle = truthUtil.GetMCParticleFromRecoShower(clockData, *daughterShower, evt, fShowerTag);
  }
  else{
    if(fVerbose > 0){
      std::cout << "INFO::Daughter pfParticle is not track or shower!" << std::endl;
    }
    //return;
  }
  
  fdaughterID[fNDAUGHTERS]                          = daughterID;
  // NHits associated with this pfParticle
  fdaughterNHits[fNDAUGHTERS]                       = (pfpUtil.GetPFParticleHits(*daughterParticle,evt,fPFParticleTag)).size();
  // T0
  std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*daughterParticle,evt,fPFParticleTag);
  if(!pfdaughterT0vec.empty())
    fdaughterT0[fNDAUGHTERS] = pfdaughterT0vec[0].Time();

  if(mcdaughterparticle != 0x0){
    fdaughter_truth_TrackId[fNDAUGHTERS]          = mcdaughterparticle->TrackId();
    fdaughter_truth_Pdg[fNDAUGHTERS]              = mcdaughterparticle->PdgCode();
    fdaughter_truth_Mother[fNDAUGHTERS]           = mcdaughterparticle->Mother();
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
    fdaughter_truth_TotalLength[fNDAUGHTERS]      = mcdaughterparticle->Trajectory().TotalLength();
    fdaughter_truth_Process[fNDAUGHTERS]          = truthUtil.GetProcessKey(mcdaughterparticle->Process());
    fdaughter_truth_EndProcess[fNDAUGHTERS]       = truthUtil.GetProcessKey(mcdaughterparticle->EndProcess());
  }
  else{
    if(fVerbose > 0){
      std::cout << "INFO::No MCParticle for daughter found!" << std::endl;
    }
  }
  
  // Increment counter
  fNDAUGHTERS++;

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::FillPrimaryGrandDaughterPFParticle(art::Event const & evt,
                                                                         detinfo::DetectorClocksData const& clockData,
                                                                         detinfo::DetectorPropertiesData const& detProp,
                                                                         const recob::PFParticle* gdaughterParticle,
                                                                         int daughterID, int gdaughterID){

  const recob::Track* gdaughterTrack                       = pfpUtil.GetPFParticleTrack(*gdaughterParticle,evt, fPFParticleTag,fTrackerTag);
  const recob::Shower* gdaughterShower                     = pfpUtil.GetPFParticleShower(*gdaughterParticle,evt,fPFParticleTag,fShowerTag);
  const simb::MCParticle* mcdaughterparticle               = NULL;

  if(gdaughterTrack != 0x0){
    fgranddaughterIstrack[fNGRANDDAUGHTERS]                = 1;
    fgranddaughterIsshower[fNGRANDDAUGHTERS]               = 0;
    fgranddaughterTheta[fNGRANDDAUGHTERS]                  = gdaughterTrack->Theta();
    fgranddaughterPhi[fNGRANDDAUGHTERS]                    = gdaughterTrack->Phi();
    fgranddaughterLength[fNGRANDDAUGHTERS]                 = gdaughterTrack->Length();
    fgranddaughterMomentum[fNGRANDDAUGHTERS]               = gdaughterTrack->StartMomentum();
    fgranddaughterEndMomentum[fNGRANDDAUGHTERS]            = gdaughterTrack->EndMomentum();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][0]       = gdaughterTrack->Trajectory().Start().X();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][1]       = gdaughterTrack->Trajectory().Start().Y();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][2]       = gdaughterTrack->Trajectory().Start().Z();
    fgranddaughterEndPosition[fNGRANDDAUGHTERS][0]         = gdaughterTrack->Trajectory().End().X();
    fgranddaughterEndPosition[fNGRANDDAUGHTERS][1]         = gdaughterTrack->Trajectory().End().Y();
    fgranddaughterEndPosition[fNGRANDDAUGHTERS][2]         = gdaughterTrack->Trajectory().End().Z();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][0]      = gdaughterTrack->Trajectory().StartDirection().X();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][1]      = gdaughterTrack->Trajectory().StartDirection().Y();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][2]      = gdaughterTrack->Trajectory().StartDirection().Z();
    fgranddaughterEndDirection[fNGRANDDAUGHTERS][0]        = gdaughterTrack->Trajectory().EndDirection().X();
    fgranddaughterEndDirection[fNGRANDDAUGHTERS][1]        = gdaughterTrack->Trajectory().EndDirection().Y();
    fgranddaughterEndDirection[fNGRANDDAUGHTERS][2]        = gdaughterTrack->Trajectory().EndDirection().Z();
    
    fgranddaughterMomentumByRangeMuon[fNGRANDDAUGHTERS]    = trmom.GetTrackMomentum(gdaughterTrack->Length(),13);
    fgranddaughterMomentumByRangeProton[fNGRANDDAUGHTERS]  = trmom.GetTrackMomentum(gdaughterTrack->Length(),2212);
    
    // Calorimetry
    std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*gdaughterTrack, evt, fTrackerTag, fCalorimetryTag);
    if(daughtercalovector.size() != 3 && fVerbose > 0)
      std::cerr << "WARNING::Calorimetry vector size for grand-daughter is = " << daughtercalovector.size() << std::endl;
    
    for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
      int plane = daughtercalovector[k].PlaneID().Plane;
      if(plane < 0) continue;
      if(plane > 2) continue;
      fgranddaughterKineticEnergy[fNGRANDDAUGHTERS][plane] = daughtercalovector[k].KineticEnergy();
      fgranddaughterRange[fNGRANDDAUGHTERS][plane]         = daughtercalovector[k].Range();
      fgranddaughterTrkPitchC[fNGRANDDAUGHTERS][plane]     = daughtercalovector[k].TrkPitchC();
    }

    // PID
    std::vector<anab::ParticleID> daughterpids = trackUtil.GetRecoTrackPID(*gdaughterTrack, evt, fTrackerTag, fParticleIDTag);

    // Set dummy values
    double pidpdg[3] = {-1,-1,-1};
    double pidchi[3] = {99999.,99999.,99999.};
    double pid2ndchi[3] = {99999.,99999.,99999.};
    int pidndf[3] = {-1,-1,-1};

    for (size_t ipid=0; ipid<daughterpids.size(); ipid++){ 
      std::vector<anab::sParticleIDAlgScores> AlgScoresVec = daughterpids[ipid].ParticleIDAlgScores();

      // Loop though AlgScoresVec and find the variables we want

      for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
        anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
        /*std::cout << "\n ParticleIDAlg " << AlgScore.fAlgName
          << "\n -- Variable type: " << AlgScore.fVariableType
          << "\n -- Track direction: " << AlgScore.fTrackDir
          << "\n -- Assuming PDG: " << AlgScore.fAssumedPdg
          << "\n -- Number of degrees of freedom: " << AlgScore.fNdf
          << "\n -- Value: " << AlgScore.fValue
          << "\n -- Using planeMask: " << AlgScore.fPlaneMask << std::endl;*/
        if (AlgScore.fPlaneMask.none() || AlgScore.fPlaneMask.count() > 1) continue;

        int planenum = -1;
        if (AlgScore.fPlaneMask.test(0)) planenum = 0;
        if (AlgScore.fPlaneMask.test(1)) planenum = 1;
        if (AlgScore.fPlaneMask.test(2)) planenum = 2;

        if (planenum<0 || planenum>2) continue;
        if (AlgScore.fAlgName == "Chi2"){
          if (TMath::Abs(AlgScore.fAssumedPdg) == 13){ // chi2mu
            fgranddaughterPID_Chi2Muon[fNGRANDDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 2212){ // chi2pr
            fgranddaughterPID_Chi2Proton[fNGRANDDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 211){ // chi2pi
            fgranddaughterPID_Chi2Pion[fNGRANDDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
          else if (TMath::Abs(AlgScore.fAssumedPdg) == 321){ // chi2ka
            fgranddaughterPID_Chi2Kaon[fNGRANDDAUGHTERS][planenum] = AlgScore.fValue;
            if (AlgScore.fValue<pidchi[planenum]){
              pidchi[planenum] = AlgScore.fValue;
              pidpdg[planenum] = TMath::Abs(AlgScore.fAssumedPdg);
            }
            else if (AlgScore.fValue<pid2ndchi[planenum]){
              pid2ndchi[planenum] = AlgScore.fValue;
            }
          }
        }
        else if (AlgScore.fVariableType==anab::kPIDA){
          fgranddaughterPID_PIDA[fNGRANDDAUGHTERS][planenum] = AlgScore.fValue;
        }
      } // end loop though AlgScoresVec
    } // end loop over pid[ipid]

      // Finally, set min chi2
    for (size_t plane=0; plane<3; plane++){
      fgranddaughterPID_MinChi2[fNGRANDDAUGHTERS][plane] = pidchi[plane];
      fgranddaughterPID_Pdg[fNGRANDDAUGHTERS][plane] = pidpdg[plane];
      fgranddaughterPID_Ndf[fNGRANDDAUGHTERS][plane] = pidndf[plane];
      fgranddaughterPID_DeltaChi2[fNGRANDDAUGHTERS][plane] = pid2ndchi[plane]-pidchi[plane];
    }

    // NOT SET with new PID object
    //fgranddaughterPID_MissingE[fNGRANDDAUGHTERS][plane]       = daughterpids[plane].MissingE();
    //fgranddaughterPID_MissingEavg[fNGRANDDAUGHTERS][plane]    = daughterpids[plane].MissingEavg();
    
    // Get the true mc particle
    //const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(*gdaughterTrack, evt, fTrackerTag);
    mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *gdaughterTrack, evt, fTrackerTag);
  }
  else if(gdaughterShower != 0x0){
    fgranddaughterIstrack[fNGRANDDAUGHTERS]                   = 0;
    fgranddaughterIsshower[fNGRANDDAUGHTERS]                  = 1;
    fgranddaughterLength[fNGRANDDAUGHTERS]                    = gdaughterShower->Length();
    fgranddaughterShowerBestPlane[fNGRANDDAUGHTERS]           = gdaughterShower->best_plane();
    fgranddaughterOpeningAngle[fNGRANDDAUGHTERS]              = gdaughterShower->OpenAngle();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][0]          = gdaughterShower->ShowerStart().X();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][1]          = gdaughterShower->ShowerStart().Y();
    fgranddaughterStartPosition[fNGRANDDAUGHTERS][2]          = gdaughterShower->ShowerStart().Z();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][0]         = gdaughterShower->Direction().X();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][1]         = gdaughterShower->Direction().Y();
    fgranddaughterStartDirection[fNGRANDDAUGHTERS][2]         = gdaughterShower->Direction().Z();
    for(size_t k = 0; k < (gdaughterShower->Energy()).size() && k<3; k++)
      fgranddaughterShowerEnergy[fNGRANDDAUGHTERS][k] = gdaughterShower->Energy()[k];
    for(size_t k = 0; k < (gdaughterShower->MIPEnergy()).size() && k<3; k++)
      fgranddaughterShowerMIPEnergy[fNGRANDDAUGHTERS][k] = gdaughterShower->MIPEnergy()[k];
    for(size_t k = 0; k < (gdaughterShower->dEdx()).size() && k<3; k++)
      fgranddaughterShowerdEdx[fNGRANDDAUGHTERS][k] = gdaughterShower->dEdx()[k];

    // Get the true mc particle
    //const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoShower(*gdaughterShower, evt, fShowerTag);
    mcdaughterparticle = truthUtil.GetMCParticleFromRecoShower(clockData, *gdaughterShower, evt, fShowerTag);
  }
  else{
    if(fVerbose > 0){
      std::cout << "INFO::GrandDaughter pfParticle is not track or shower!" << std::endl;
    }
    //return;
  }
  
  fgranddaughterID[fNGRANDDAUGHTERS]                          = gdaughterID;
  fgranddaughterMotherID[fNGRANDDAUGHTERS]                    = daughterID;
  // NHits associated with this pfParticle
  fgranddaughterNHits[fNGRANDDAUGHTERS]                       = (pfpUtil.GetPFParticleHits(*gdaughterParticle,evt,fPFParticleTag)).size();
  // T0
  std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*gdaughterParticle,evt,fPFParticleTag);
  if(!pfdaughterT0vec.empty())
    fgranddaughterT0[fNGRANDDAUGHTERS] = pfdaughterT0vec[0].Time();

  if(mcdaughterparticle != 0x0){
    fgranddaughter_truth_TrackId[fNGRANDDAUGHTERS]          = mcdaughterparticle->TrackId();
    fgranddaughter_truth_Pdg[fNGRANDDAUGHTERS]              = mcdaughterparticle->PdgCode();
    fgranddaughter_truth_Mother[fNGRANDDAUGHTERS]           = mcdaughterparticle->Mother();
    fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][0] = mcdaughterparticle->Vx();
    fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][1] = mcdaughterparticle->Vy();
    fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][2] = mcdaughterparticle->Vz();
    fgranddaughter_truth_StartPosition[fNGRANDDAUGHTERS][3] = mcdaughterparticle->T();
    fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][0]   = mcdaughterparticle->EndX();
    fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][1]   = mcdaughterparticle->EndY();
    fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][2]   = mcdaughterparticle->EndZ();
    fgranddaughter_truth_EndPosition[fNGRANDDAUGHTERS][3]   = mcdaughterparticle->EndT();
    fgranddaughter_truth_P[fNGRANDDAUGHTERS]                = mcdaughterparticle->P();
    fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][0]      = mcdaughterparticle->Px();
    fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][1]      = mcdaughterparticle->Py();
    fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][2]      = mcdaughterparticle->Pz();
    fgranddaughter_truth_Momentum[fNGRANDDAUGHTERS][3]      = mcdaughterparticle->E();
    fgranddaughter_truth_Pt[fNGRANDDAUGHTERS]               = mcdaughterparticle->Pt();
    fgranddaughter_truth_Mass[fNGRANDDAUGHTERS]             = mcdaughterparticle->Mass();
    fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][0]   = mcdaughterparticle->EndPx();
    fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][1]   = mcdaughterparticle->EndPy();
    fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][2]   = mcdaughterparticle->EndPz();
    fgranddaughter_truth_EndMomentum[fNGRANDDAUGHTERS][3]   = mcdaughterparticle->EndE();
    fgranddaughter_truth_Theta[fNGRANDDAUGHTERS]            = mcdaughterparticle->Momentum().Theta();
    fgranddaughter_truth_Phi[fNGRANDDAUGHTERS]              = mcdaughterparticle->Momentum().Phi();
    fgranddaughter_truth_TotalLength[fNGRANDDAUGHTERS]      = mcdaughterparticle->Trajectory().TotalLength();
    fgranddaughter_truth_Process[fNGRANDDAUGHTERS]          = truthUtil.GetProcessKey(mcdaughterparticle->Process());
    fgranddaughter_truth_EndProcess[fNGRANDDAUGHTERS]       = truthUtil.GetProcessKey(mcdaughterparticle->EndProcess());
  }
  else{
    if(fVerbose > 0){
      std::cout << "INFO::No MCParticle for granddaughter found!" << std::endl;
    }
  }
  
  // Increment counter
  fNGRANDDAUGHTERS++;

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::FillConfigTree(){

  fNCryostats     = fGeometry->Ncryostats();
  fNTPCs          = fGeometry->NTPC();
  fNChannels      = fWireReadoutGeom->Nchannels();
  fNPlanes        = fWireReadoutGeom->Nplanes();
  fNAPAs          = fGeometry->NTPC()*fGeometry->Ncryostats()/2;
  if(fNAPAs > 0)
    fNChansPerAPA = fWireReadoutGeom->Nchannels()/fNAPAs;

  fActiveTPCBoundsX[0] = fActiveTPCBoundsY[0] = fActiveTPCBoundsZ[0] = 1000000.0;
  fActiveTPCBoundsX[1] = fActiveTPCBoundsY[1] = fActiveTPCBoundsZ[1] = -1000000.0;

  if(fVerbose > 0){
    std::cout << "Detector Name: " << fGeometry->DetectorName() << std::endl;
    std::cout << "GDML file: " << fGeometry->GDMLFile() << std::endl;
  }

  for (geo::TPCGeo const& TPC: fGeometry->Iterate<geo::TPCGeo>()) {
    // get center in world coordinates
    auto const center = TPC.GetCenter();
    double tpc[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
    if(center.X() - tpc[0] < fActiveTPCBoundsX[0]) fActiveTPCBoundsX[0] = center.X() - tpc[0];
    if(center.X() + tpc[0] > fActiveTPCBoundsX[1]) fActiveTPCBoundsX[1] = center.X() + tpc[0];
    if(center.Y() - tpc[1] < fActiveTPCBoundsY[0]) fActiveTPCBoundsY[0] = center.Y() - tpc[1];
    if(center.Y() + tpc[1] > fActiveTPCBoundsY[1]) fActiveTPCBoundsY[1] = center.Y() + tpc[1];
    if(center.Z() - tpc[2] < fActiveTPCBoundsZ[0]) fActiveTPCBoundsZ[0] = center.Z() - tpc[2];
    if(center.Z() + tpc[2] > fActiveTPCBoundsZ[1]) fActiveTPCBoundsZ[1] = center.Z() + tpc[2];
  }

  fConfigTree->Fill();

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  // Configuration Tree
  fConfigTree = tfs->make<TTree>("Config", "Configuration Tree");
  fConfigTree->Branch("NAPAs",                                         &fNAPAs,                                      "NAPAs/I");
  fConfigTree->Branch("NChansPerAPA",                                  &fNChansPerAPA,                               "NChansPerAPA/I");
  fConfigTree->Branch("NCryostats",                                    &fNCryostats,                                 "NCryostats/I");
  fConfigTree->Branch("NTPCs",                                         &fNTPCs,                                      "NTPCs/I");
  fConfigTree->Branch("NChannels",                                     &fNChannels,                                  "NChannels/I");
  fConfigTree->Branch("NPlanes",                                       &fNPlanes,                                    "NPlanes/I");
  fConfigTree->Branch("ActiveTPCBoundsX",                              &fActiveTPCBoundsX,                           "ActiveTPCBoundsX[2]/D");
  fConfigTree->Branch("ActiveTPCBoundsY",                              &fActiveTPCBoundsY,                           "ActiveTPCBoundsY[2]/D");
  fConfigTree->Branch("ActiveTPCBoundsZ",                              &fActiveTPCBoundsZ,                           "ActiveTPCBoundsZ[2]/D");

  // Fill configuration tree
  FillConfigTree();

  fBeamTrack = tfs->make<TTree>("BeamTrack", "Beam track tree");
  fBeamTrack->Branch("run",                                            &fRun,                                        "run/I");
  fBeamTrack->Branch("subrun",                                         &fSubRun,                                     "subrun/I");
  fBeamTrack->Branch("event",                                          &fevent,                                      "event/I");
  fBeamTrack->Branch("beamNTracks",                                    &fbeamNTracks,                                "beamNTracks/I");
  fBeamTrack->Branch("beamtrigger",                                    &fbeamtrigger,                                "beamtrigger/I");
  fBeamTrack->Branch("tof",                                            &ftof,                                        "tof/D");
  fBeamTrack->Branch("tof_expElec",                                    &ftof_expElec,                                "tof_expElec/D");
  fBeamTrack->Branch("tof_expMuon",                                    &ftof_expMuon,                                "tof_expMuon/D");
  fBeamTrack->Branch("tof_expPion",                                    &ftof_expPion,                                "tof_expPion/D");
  fBeamTrack->Branch("tof_expKaon",                                    &ftof_expKaon,                                "tof_expKaon/D");
  fBeamTrack->Branch("tof_expProt",                                    &ftof_expProt,                                "tof_expProt/D");
  fBeamTrack->Branch("tof_expDeut",                                    &ftof_expDeut,                                "tof_expDeut/D");
  fBeamTrack->Branch("cerenkovStatus",                                 &fcerenkovStatus,                             "cerenkovStatus[2]/I");
  fBeamTrack->Branch("cerenkovTime",                                   &fcerenkovTime,                               "cerenkovTime[2]/D");
  fBeamTrack->Branch("cerenkovPressure",                               &fcerenkovPressure,                           "cerenkovPressure[2]/D");
  fBeamTrack->Branch("beamtrackMomentum",                              &fbeamtrackMomentum,                          "beamtrackMomentum/D");
  fBeamTrack->Branch("beamtrackEnergy",                                &fbeamtrackEnergy,                            "beamtrackEnergy/D");
  fBeamTrack->Branch("beamtrackPos",                                   &fbeamtrackPos,                               "beamtrackPos[4]/D");
  fBeamTrack->Branch("beamtrackDir",                                   &fbeamtrackDir,                               "beamtrackDir[3]/D");
  fBeamTrack->Branch("BIAndTimingMatched",                             &fBIAndTimingMatched,                         "BIAndTimingMatched/I");
  
  fBeamTrack->Branch("beamtrack_truth_Origin",                         &fbeamtrack_truth_Origin,                     "beamtrack_truth_Origin/I");
  fBeamTrack->Branch("beamtrackPdg",                                   &fbeamtrackPdg,                               "beamtrackPdg/I");
  fBeamTrack->Branch("beamtrackID",                                    &fbeamtrackID,                                "beamtrackID/I");
  fBeamTrack->Branch("beamtrack_truth_EndPos",                         &fbeamtrack_truth_EndPos,                     "beamtrack_truth_EndPos[4]/D");
  fBeamTrack->Branch("beamtrack_truth_Momentum",                       &fbeamtrack_truth_Momentum,                   "beamtrack_truth_Momentum[4]/D");
  fBeamTrack->Branch("beamtrack_truth_KinEnergy",                      &fbeamtrack_truth_KinEnergy,                  "beamtrack_truth_KinEnergy/D");
  fBeamTrack->Branch("beamtrack_truth_Pt",                             &fbeamtrack_truth_Pt,                         "beamtrack_truth_Pt/D");
  fBeamTrack->Branch("beamtrack_truth_Mass",                           &fbeamtrack_truth_Mass,                       "beamtrack_truth_Mass/D");
  fBeamTrack->Branch("beamtrack_truth_Theta",                          &fbeamtrack_truth_Theta,                      "beamtrack_truth_Theta/D");
  fBeamTrack->Branch("beamtrack_truth_Phi",                            &fbeamtrack_truth_Phi,                        "beamtrack_truth_Phi/D");
  fBeamTrack->Branch("beamtrack_truth_TotalLength",                    &fbeamtrack_truth_TotalLength,                "beamtrack_truth_TotalLength/D");
  fBeamTrack->Branch("beamtrack_truth_Process",                        &fbeamtrack_truth_Process,                    "beamtrack_truth_Process/I");
  fBeamTrack->Branch("beamtrack_truth_EndProcess",                     &fbeamtrack_truth_EndProcess,                 "beamtrack_truth_EndProcess/I");
  fBeamTrack->Branch("beamtrack_truth_Pdg",                            &fbeamtrack_truth_Pdg,                        "beamtrack_truth_Pdg/I");
  fBeamTrack->Branch("beamtrack_truth_NDAUGTHERS",                     &fbeamtrack_truth_NDAUGTHERS,                 "beamtrack_truth_NDAUGTHERS/I");
  fBeamTrack->Branch("beamtrack_truth_NDECAYDAUGTHERS",                &fbeamtrack_truth_NDECAYDAUGTHERS,            "beamtrack_truth_NDECAYDAUGTHERS/I");
  fBeamTrack->Branch("beamtrack_truth_Pos_InTPCActive",                &fbeamtrack_truth_Pos_InTPCActive,            "beamtrack_truth_Pos_InTPCActive[4]/D");
  fBeamTrack->Branch("beamtrack_truth_Momentum_InTPCActive",           &fbeamtrack_truth_Momentum_InTPCActive,       "beamtrack_truth_Momentum_InTPCActive[4]/D");
  fBeamTrack->Branch("beamtrack_truth_P_InTPCActive",                  &fbeamtrack_truth_P_InTPCActive,              "beamtrack_truth_P_InTPCActive/D");
  fBeamTrack->Branch("beamtrack_truth_Pt_InTPCActive",                 &fbeamtrack_truth_Pt_InTPCActive,             "beamtrack_truth_Pt_InTPCActive/D");
  fBeamTrack->Branch("beamtrack_truth_Theta_InTPCActive",              &fbeamtrack_truth_Theta_InTPCActive,          "beamtrack_truth_Theta_InTPCActive/D");
  fBeamTrack->Branch("beamtrack_truth_Phi_InTPCActive",                &fbeamtrack_truth_Phi_InTPCActive,            "beamtrack_truth_Phi_InTPCActive/D");
  fBeamTrack->Branch("beamtrack_truth_TotalLength_InTPCActive",        &fbeamtrack_truth_TotalLength_InTPCActive,    "beamtrack_truth_TotalLength_InTPCActive/D");
  fBeamTrack->Branch("beamtrack_truth_KinEnergy_InTPCActive",          &fbeamtrack_truth_KinEnergy_InTPCActive,      "beamtrack_truth_KinEnergy_InTPCActive/D");

  fBeamTrack->Branch("beamtrack_truthdaughter_TrackId",                &fbeamtrack_truthdaughter_TrackId,            "beamtrack_truthdaughter_TrackId[beamtrack_truth_NDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdaughter_Pdg",                    &fbeamtrack_truthdaughter_Pdg,                "beamtrack_truthdaughter_Pdg[beamtrack_truth_NDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdaughter_Mother",                 &fbeamtrack_truthdaughter_Mother,             "beamtrack_truthdaughter_Mother[beamtrack_truth_NDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdaughter_StartPosition",          &fbeamtrack_truthdaughter_StartPosition,      "beamtrack_truthdaughter_StartPosition[beamtrack_truth_NDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_EndPosition",            &fbeamtrack_truthdaughter_EndPosition,        "beamtrack_truthdaughter_EndPosition[beamtrack_truth_NDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_P",                      &fbeamtrack_truthdaughter_P,                  "beamtrack_truthdaughter_P[beamtrack_truth_NDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_Momentum",               &fbeamtrack_truthdaughter_Momentum,           "beamtrack_truthdaughter_Momentum[beamtrack_truth_NDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_EndMomentum",            &fbeamtrack_truthdaughter_EndMomentum,        "beamtrack_truthdaughter_EndMomentum[beamtrack_truth_NDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_Pt",                     &fbeamtrack_truthdaughter_Pt,                 "beamtrack_truthdaughter_Pt[beamtrack_truth_NDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_Mass",                   &fbeamtrack_truthdaughter_Mass,               "beamtrack_truthdaughter_Mass[beamtrack_truth_NDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_Theta",                  &fbeamtrack_truthdaughter_Theta,              "beamtrack_truthdaughter_Theta[beamtrack_truth_NDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_Phi",                    &fbeamtrack_truthdaughter_Phi,                "beamtrack_truthdaughter_Phi[beamtrack_truth_NDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_TotalLength",            &fbeamtrack_truthdaughter_TotalLength,        "beamtrack_truthdaughter_TotalLength[beamtrack_truth_NDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdaughter_Process",                &fbeamtrack_truthdaughter_Process,            "beamtrack_truthdaughter_Process[beamtrack_truth_NDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdaughter_EndProcess",             &fbeamtrack_truthdaughter_EndProcess,         "beamtrack_truthdaughter_EndProcess[beamtrack_truth_NDAUGTHERS]/I");

  fBeamTrack->Branch("beamtrack_truthdecaydaughter_TrackId",           &fbeamtrack_truthdecaydaughter_TrackId,       "beamtrack_truthdecaydaughter_TrackId[beamtrack_truth_NDECAYDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Pdg",               &fbeamtrack_truthdecaydaughter_Pdg,           "beamtrack_truthdecaydaughter_Pdg[beamtrack_truth_NDECAYDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Mother",            &fbeamtrack_truthdecaydaughter_Mother,        "beamtrack_truthdecaydaughter_Mother[beamtrack_truth_NDECAYDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_StartPosition",     &fbeamtrack_truthdecaydaughter_StartPosition, "beamtrack_truthdecaydaughter_StartPosition[beamtrack_truth_NDECAYDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_EndPosition",       &fbeamtrack_truthdecaydaughter_EndPosition,   "beamtrack_truthdecaydaughter_EndPosition[beamtrack_truth_NDECAYDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_P",                 &fbeamtrack_truthdecaydaughter_P,             "beamtrack_truthdecaydaughter_P[beamtrack_truth_NDECAYDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Momentum",          &fbeamtrack_truthdecaydaughter_Momentum,      "beamtrack_truthdecaydaughter_Momentum[beamtrack_truth_NDECAYDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_EndMomentum",       &fbeamtrack_truthdecaydaughter_EndMomentum,   "beamtrack_truthdecaydaughter_EndMomentum[beamtrack_truth_NDECAYDAUGTHERS][4]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Pt",                &fbeamtrack_truthdecaydaughter_Pt,            "beamtrack_truthdecaydaughter_Pt[beamtrack_truth_NDECAYDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Mass",              &fbeamtrack_truthdecaydaughter_Mass,          "beamtrack_truthdecaydaughter_Mass[beamtrack_truth_NDECAYDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Theta",             &fbeamtrack_truthdecaydaughter_Theta,         "beamtrack_truthdecaydaughter_Theta[beamtrack_truth_NDECAYDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Phi",               &fbeamtrack_truthdecaydaughter_Phi,           "beamtrack_truthdecaydaughter_Phi[beamtrack_truth_NDECAYDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_TotalLength",       &fbeamtrack_truthdecaydaughter_TotalLength,   "beamtrack_truthdecaydaughter_TotalLength[beamtrack_truth_NDECAYDAUGTHERS]/D");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_Process",           &fbeamtrack_truthdecaydaughter_Process,       "beamtrack_truthdecaydaughter_Process[beamtrack_truth_NDECAYDAUGTHERS]/I");
  fBeamTrack->Branch("beamtrack_truthdecaydaughter_EndProcess",        &fbeamtrack_truthdecaydaughter_EndProcess,    "beamtrack_truthdecaydaughter_EndProcess[beamtrack_truth_NDECAYDAUGTHERS]/I");

  // Analysis tree
  fPandoraBeam = tfs->make<TTree>("ReconstructedBeamEvents", "Reconstructed Beam events");
  fPandoraBeam->Branch("run",                                          &fRun,                                        "run/I");
  fPandoraBeam->Branch("subrun",                                       &fSubRun,                                     "subrun/I");
  fPandoraBeam->Branch("event",                                        &fevent,                                      "event/I");
  fPandoraBeam->Branch("timestamp",                                    &fTimeStamp,                                  "timestamp/D");
  fPandoraBeam->Branch("Nactivefembs",                                 &fNactivefembs,                               "Nactivefembs[5]/I");
  fPandoraBeam->Branch("NpfParticles",                                 &fNpfParticles,                               "NpfParticles/I");
  fPandoraBeam->Branch("primaryEarliestHitPeakTim",                    &fprimaryEarliestHitPeakTime,                 "primaryEarliestHitPeakTim/D");
  fPandoraBeam->Branch("primaryVertex",                                &fprimaryVertex,                              "primaryVertex[3]/D");
  fPandoraBeam->Branch("primaryIstrack",                               &fprimaryIstrack,                             "primaryIstrack/I");
  fPandoraBeam->Branch("primaryIsshower",                              &fprimaryIsshower,                            "primaryIsshower/I");
  fPandoraBeam->Branch("primaryBDTScore",                              &fprimaryBDTScore,                            "primaryBDTScore/D");
  fPandoraBeam->Branch("primaryNHits",                                 &fprimaryNHits,                               "primaryNHits/I");
  fPandoraBeam->Branch("primaryTheta",                                 &fprimaryTheta,                               "primaryTheta/D");
  fPandoraBeam->Branch("primaryPhi",                                   &fprimaryPhi,                                 "primaryPhi/D");
  fPandoraBeam->Branch("primaryLength",                                &fprimaryLength,                              "primaryLength/D");
  fPandoraBeam->Branch("primaryMomentum",                              &fprimaryMomentum,                            "primaryMomentum/D");
  fPandoraBeam->Branch("primaryEndMomentum",                           &fprimaryEndMomentum,                         "primaryEndMomentum/D");
  fPandoraBeam->Branch("primaryEndPosition",                           &fprimaryEndPosition,                         "primaryEndPosition[3]/D");
  fPandoraBeam->Branch("primaryStartPosition",                         &fprimaryStartPosition,                       "primaryStartPosition[3]/D");
  fPandoraBeam->Branch("primaryEndDirection",                          &fprimaryEndDirection,                        "primaryEndDirection[3]/D");
  fPandoraBeam->Branch("primaryStartDirection",                        &fprimaryStartDirection,                      "primaryStartDirection[3]/D");
  fPandoraBeam->Branch("primaryOpeningAngle",                          &fprimaryOpeningAngle,                        "primaryOpeningAngle/D");
  fPandoraBeam->Branch("primaryID",                                    &fprimaryID,                                  "primaryID/I");
  fPandoraBeam->Branch("primaryShowerBestPlane",                       &fprimaryShowerBestPlane,                     "primaryShowerBestPlane/I");
  fPandoraBeam->Branch("primaryShowerEnergy",                          &fprimaryShowerEnergy,                        "primaryShowerEnergy[3]/I");
  fPandoraBeam->Branch("primaryShowerMIPEnergy",                       &fprimaryShowerMIPEnergy,                     "primaryShowerMIPEnergy[3]/I");
  fPandoraBeam->Branch("primaryShowerdEdx",                            &fprimaryShowerdEdx,                          "primaryShowerdEdx[3]/I");
  fPandoraBeam->Branch("primaryMomentumByRangeProton",                 &fprimaryMomentumByRangeProton,               "primaryMomentumByRangeProton/D");
  fPandoraBeam->Branch("primaryMomentumByRangeMuon",                   &fprimaryMomentumByRangeMuon,                 "primaryMomentumByRangeMuon/D");
  fPandoraBeam->Branch("primaryKineticEnergy",                         &fprimaryKineticEnergy,                       "primaryKineticEnergy[3]/D");
  fPandoraBeam->Branch("primaryRange",                                 &fprimaryRange,                               "primaryRange[3]/D");
  fPandoraBeam->Branch("primaryTrkPitchC",                             &fprimaryTrkPitchC,                           "primaryTrkPitchC[3]/D");
  fPandoraBeam->Branch("primaryT0",                                    &fprimaryT0,                                  "primaryT0/D");
  fPandoraBeam->Branch("primaryPID_Pdg",                               &fprimaryPID_Pdg,                             "primaryPID_Pdg[3]/I");
  fPandoraBeam->Branch("primaryPID_Ndf",                               &fprimaryPID_Ndf,                             "primaryPID_Ndf[3]/I");
  fPandoraBeam->Branch("primaryPID_MinChi2",                           &fprimaryPID_MinChi2,                         "primaryPID_MinChi2[3]/D");
  fPandoraBeam->Branch("primaryPID_DeltaChi2",                         &fprimaryPID_DeltaChi2,                       "primaryPID_DeltaChi2[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Proton",                        &fprimaryPID_Chi2Proton,                      "primaryPID_Chi2Proton[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Kaon",                          &fprimaryPID_Chi2Kaon,                        "primaryPID_Chi2Kaon[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Pion",                          &fprimaryPID_Chi2Pion,                        "primaryPID_Chi2Pion[3]/D");
  fPandoraBeam->Branch("primaryPID_Chi2Muon",                          &fprimaryPID_Chi2Muon,                        "primaryPID_Chi2Muon[3]/D");
  fPandoraBeam->Branch("primaryPID_MissingE",                          &fprimaryPID_MissingE,                        "primaryPID_MissingE[3]/D");
  fPandoraBeam->Branch("primaryPID_MissingEavg",                       &fprimaryPID_MissingEavg,                     "primaryPID_MissingEavg[3]/D");
  fPandoraBeam->Branch("primaryPID_PIDA",                              &fprimaryPID_PIDA,                            "primaryPID_PIDA[3]/D");

  fPandoraBeam->Branch("primary_truth_Origin",                         &fprimary_truth_Origin,                       "primary_truth_Origin/I");
  fPandoraBeam->Branch("primary_truth_TrackId",                        &fprimary_truth_TrackId,                      "primary_truth_TrackId/I");
  fPandoraBeam->Branch("primary_truth_Pdg",                            &fprimary_truth_Pdg,                          "primary_truth_Pdg/I");
  fPandoraBeam->Branch("primary_truth_Mother",                         &fprimary_truth_Mother,                       "primary_truth_Mother/I");
  fPandoraBeam->Branch("primary_truth_StartPosition",                  &fprimary_truth_StartPosition,                "primary_truth_StartPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_EndPosition",                    &fprimary_truth_EndPosition,                  "primary_truth_EndPosition[4]/D");
  fPandoraBeam->Branch("primary_truth_Momentum",                       &fprimary_truth_Momentum,                     "primary_truth_Momentum[4]/D");
  fPandoraBeam->Branch("primary_truth_EndMomentum",                    &fprimary_truth_EndMomentum,                  "primary_truth_EndMomentum[4]/D");
  fPandoraBeam->Branch("primary_truth_P",                              &fprimary_truth_P,                            "primary_truth_P/D");
  fPandoraBeam->Branch("primary_truth_Pt",                             &fprimary_truth_Pt,                           "primary_truth_Pt/D");
  fPandoraBeam->Branch("primary_truth_Mass",                           &fprimary_truth_Mass,                         "primary_truth_Mass/D");
  fPandoraBeam->Branch("primary_truth_Theta",                          &fprimary_truth_Theta,                        "primary_truth_Theta/D");
  fPandoraBeam->Branch("primary_truth_Phi",                            &fprimary_truth_Phi,                          "primary_truth_Phi/D");
  fPandoraBeam->Branch("primary_truth_TotalLength",                    &fprimary_truth_TotalLength,                  "primary_truth_TotalLength/D");
  fPandoraBeam->Branch("primary_truth_Process",                        &fprimary_truth_Process,                      "primary_truth_Process/I");
  fPandoraBeam->Branch("primary_truth_EndProcess",                     &fprimary_truth_EndProcess,                   "primary_truth_EndProcess/I");
  fPandoraBeam->Branch("primary_truth_Isbeammatched",                  &fprimary_truth_Isbeammatched,                "primary_truth_Isbeammatched/I");

  fPandoraBeam->Branch("primary_truth_EkinAtVertex_notcorrected",      &fprimary_truth_EkinAtVertex_notcorrected,    "primary_truth_EkinAtVertex_notcorrected/D");
  fPandoraBeam->Branch("primary_truth_EkinAtVertex",                   &fprimary_truth_EkinAtVertex,                 "primary_truth_EkinAtVertex/D");
  fPandoraBeam->Branch("primary_truth_Pos_InTPCActive",                &fprimary_truth_Pos_InTPCActive,              "primary_truth_Pos_InTPCActive[4]/D");
  fPandoraBeam->Branch("primary_truth_Momentum_InTPCActive",           &fprimary_truth_Momentum_InTPCActive,         "primary_truth_Momentum_InTPCActive[4]/D");
  fPandoraBeam->Branch("primary_truth_P_InTPCActive",                  &fprimary_truth_P_InTPCActive,                "primary_truth_P_InTPCActive/D");
  fPandoraBeam->Branch("primary_truth_Pt_InTPCActive",                 &fprimary_truth_Pt_InTPCActive,               "primary_truth_Pt_InTPCActive/D");
  fPandoraBeam->Branch("primary_truth_Theta_InTPCActive",              &fprimary_truth_Theta_InTPCActive,            "primary_truth_Theta_InTPCActive/D");
  fPandoraBeam->Branch("primary_truth_Phi_InTPCActive",                &fprimary_truth_Phi_InTPCActive,              "primary_truth_Phi_InTPCActive/D");
  fPandoraBeam->Branch("primary_truth_TotalLength_InTPCActive",        &fprimary_truth_TotalLength_InTPCActive,      "primary_truth_TotalLength_InTPCActive/D");
  fPandoraBeam->Branch("primary_truth_KinEnergy_InTPCActive",          &fprimary_truth_KinEnergy_InTPCActive,        "primary_truth_KinEnergy_InTPCActive/D");

  fPandoraBeam->Branch("primary_truth_NDAUGTHERS",                     &fprimary_truth_NDAUGTHERS,                   "primary_truth_NDAUGTHERS/I");
  fPandoraBeam->Branch("primary_truth_NDECAYDAUGTHERS",                &fprimary_truth_NDECAYDAUGTHERS,              "primary_truth_NDECAYDAUGTHERS/I");
  fPandoraBeam->Branch("primary_truthdaughter_TrackId",                &fprimary_truthdaughter_TrackId,              "primary_truthdaughter_TrackId[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdaughter_Pdg",                    &fprimary_truthdaughter_Pdg,                  "primary_truthdaughter_Pdg[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdaughter_Mother",                 &fprimary_truthdaughter_Mother,               "primary_truthdaughter_Mother[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdaughter_StartPosition",          &fprimary_truthdaughter_StartPosition,        "primary_truthdaughter_StartPosition[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdaughter_EndPosition",            &fprimary_truthdaughter_EndPosition,          "primary_truthdaughter_EndPosition[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdaughter_P",                      &fprimary_truthdaughter_P,                    "primary_truthdaughter_P[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdaughter_Momentum",               &fprimary_truthdaughter_Momentum,             "primary_truthdaughter_Momentum[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdaughter_EndMomentum",            &fprimary_truthdaughter_EndMomentum,          "primary_truthdaughter_EndMomentum[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdaughter_Pt",                     &fprimary_truthdaughter_Pt,                   "primary_truthdaughter_Pt[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdaughter_Mass",                   &fprimary_truthdaughter_Mass,                 "primary_truthdaughter_Mass[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdaughter_Theta",                  &fprimary_truthdaughter_Theta,                "primary_truthdaughter_Theta[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdaughter_Phi",                    &fprimary_truthdaughter_Phi,                  "primary_truthdaughter_Phi[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdaughter_TotalLength",            &fprimary_truthdaughter_TotalLength,          "primary_truthdaughter_TotalLength[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdaughter_Process",                &fprimary_truthdaughter_Process,              "primary_truthdaughter_Process[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdaughter_EndProcess",             &fprimary_truthdaughter_EndProcess,           "primary_truthdaughter_EndProcess[primary_truth_NDAUGTHERS]/I");

  fPandoraBeam->Branch("primary_truthdecaydaughter_TrackId",           &fprimary_truthdecaydaughter_TrackId,         "primary_truthdecaydaughter_TrackId[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Pdg",               &fprimary_truthdecaydaughter_Pdg,             "primary_truthdecaydaughter_Pdg[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Mother",            &fprimary_truthdecaydaughter_Mother,          "primary_truthdecaydaughter_Mother[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdecaydaughter_StartPosition",     &fprimary_truthdecaydaughter_StartPosition,   "primary_truthdecaydaughter_StartPosition[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_EndPosition",       &fprimary_truthdecaydaughter_EndPosition,     "primary_truthdecaydaughter_EndPosition[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_P",                 &fprimary_truthdecaydaughter_P,               "primary_truthdecaydaughter_P[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Momentum",          &fprimary_truthdecaydaughter_Momentum,        "primary_truthdecaydaughter_Momentum[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_EndMomentum",       &fprimary_truthdecaydaughter_EndMomentum,     "primary_truthdecaydaughter_EndMomentum[primary_truth_NDAUGTHERS][4]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Pt",                &fprimary_truthdecaydaughter_Pt,              "primary_truthdecaydaughter_Pt[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Mass",              &fprimary_truthdecaydaughter_Mass,            "primary_truthdecaydaughter_Mass[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Theta",             &fprimary_truthdecaydaughter_Theta,           "primary_truthdecaydaughter_Theta[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Phi",               &fprimary_truthdecaydaughter_Phi,             "primary_truthdecaydaughter_Phi[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_TotalLength",       &fprimary_truthdecaydaughter_TotalLength,     "primary_truthdecaydaughter_TotalLength[primary_truth_NDAUGTHERS]/D");
  fPandoraBeam->Branch("primary_truthdecaydaughter_Process",           &fprimary_truthdecaydaughter_Process,         "primary_truthdecaydaughter_Process[primary_truth_NDAUGTHERS]/I");
  fPandoraBeam->Branch("primary_truthdecaydaughter_EndProcess",        &fprimary_truthdecaydaughter_EndProcess,      "primary_truthdecaydaughter_EndProcess[primary_truth_NDAUGTHERS]/I");

  fDaughterTree = tfs->make<TTree>("ReconstructedDaughters", "Reconstructed Daughters"); 
  fDaughterTree->Branch("NDAUGHTERS",                                  &fNDAUGHTERS,                                 "NDAUGHTERS/I");
  fDaughterTree->Branch("daughterVertex",                              &fdaughterVertex,                             "daughterVertex[3]/D");
  fDaughterTree->Branch("daughterIstrack",                             &fdaughterIstrack,                            "daughterIstrack[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughterIsshower",                            &fdaughterIsshower,                           "daughterIsshower[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughterNHits",                               &fdaughterNHits,                              "daughterNHits[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughterTheta",                               &fdaughterTheta,                              "daughterTheta[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterPhi",                                 &fdaughterPhi,                                "daughterPhi[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterLength",                              &fdaughterLength,                             "daughterLength[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterMomentum",                            &fdaughterMomentum,                           "daughterMomentum[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterEndMomentum",                         &fdaughterEndMomentum,                        "daughterEndMomentum[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterEndPosition",                         &fdaughterEndPosition,                        "daughterEndPosition[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterStartPosition",                       &fdaughterStartPosition,                      "daughterStartPosition[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterStartDirection",                      &fdaughterStartDirection,                     "daughterStartDirection[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterEndDirection",                        &fdaughterEndDirection,                       "daughterEndDirection[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterOpeningAngle",                        &fdaughterOpeningAngle,                       "daughterOpeningAngle[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterShowerBestPlane",                     &fdaughterShowerBestPlane,                    "daughterShowerBestPlane[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughterShowerEnergy",                        &fdaughterShowerEnergy,                       "daughterShowerEnergy[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterShowerMIPEnergy",                     &fdaughterShowerMIPEnergy,                    "daughterShowerMIPEnergy[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterShowerdEdx",                          &fdaughterShowerdEdx,                         "daughterShowerdEdx[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterMomentumByRangeProton",               &fdaughterMomentumByRangeProton,              "daughterMomentumByRangeProton[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterMomentumByRangeMuon",                 &fdaughterMomentumByRangeMuon,                "daughterMomentumByRangeMuon[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughterKineticEnergy",                       &fdaughterKineticEnergy,                      "daughterKineticEnergy[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterRange",                               &fdaughterRange,                              "daughterRange[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterTrkPitchC",                           &fdaughterTrkPitchC,                          "daughterTrkPitchC[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterID",                                  &fdaughterID,                                 "daughterID[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughterT0",                                  &fdaughterT0,                                 "daughterT0[NDAUGHTERS]/D");

  fDaughterTree->Branch("daughterPID_Pdg",                             &fdaughterPID_Pdg,                            "daughterPID_Pdg[NDAUGHTERS][3]/I");
  fDaughterTree->Branch("daughterPID_Ndf",                             &fdaughterPID_Ndf,                            "daughterPID_Ndf[NDAUGHTERS][3]/I");
  fDaughterTree->Branch("daughterPID_MinChi2",                         &fdaughterPID_MinChi2,                        "daughterPID_MinChi2[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_DeltaChi2",                       &fdaughterPID_DeltaChi2,                      "daughterPID_DeltaChi2[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_Chi2Proton",                      &fdaughterPID_Chi2Proton,                     "daughterPID_Chi2Proton[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_Chi2Kaon",                        &fdaughterPID_Chi2Kaon,                       "daughterPID_Chi2Kaon[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_Chi2Pion",                        &fdaughterPID_Chi2Pion,                       "daughterPID_Chi2Pion[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_Chi2Muon",                        &fdaughterPID_Chi2Muon,                       "daughterPID_Chi2Muon[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_MissingE",                        &fdaughterPID_MissingE,                       "daughterPID_MissingE[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_MissingEavg",                     &fdaughterPID_MissingEavg,                    "daughterPID_MissingEavg[NDAUGHTERS][3]/D");
  fDaughterTree->Branch("daughterPID_PIDA",                            &fdaughterPID_PIDA,                           "daughterPID_PIDA[NDAUGHTERS][3]/D");

  fDaughterTree->Branch("daughter_truth_TrackId",                      &fdaughter_truth_TrackId,                     "daughter_truth_TrackId[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughter_truth_Pdg",                          &fdaughter_truth_Pdg,                         "daughter_truth_Pdg[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughter_truth_Mother",                       &fdaughter_truth_Mother,                      "daughter_truth_Mother[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughter_truth_StartPosition",                &fdaughter_truth_StartPosition,               "daughter_truth_StartPosition[NDAUGHTERS][4]/D");
  fDaughterTree->Branch("daughter_truth_EndPosition",                  &fdaughter_truth_EndPosition,                 "daughter_truth_EndPosition[NDAUGHTERS][4]/D");
  fDaughterTree->Branch("daughter_truth_Momentum",                     &fdaughter_truth_Momentum,                    "daughter_truth_Momentum[NDAUGHTERS][4]/D");
  fDaughterTree->Branch("daughter_truth_EndMomentum",                  &fdaughter_truth_EndMomentum,                 "daughter_truth_EndMomentum[NDAUGHTERS][4]/D");
  fDaughterTree->Branch("daughter_truth_P",                            &fdaughter_truth_P,                           "daughter_truth_P[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughter_truth_Pt",                           &fdaughter_truth_Pt,                          "daughter_truth_Pt[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughter_truth_Mass",                         &fdaughter_truth_Mass,                        "daughter_truth_Mass[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughter_truth_Theta",                        &fdaughter_truth_Theta,                       "daughter_truth_Theta[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughter_truth_Phi",                          &fdaughter_truth_Phi,                         "daughter_truth_Phi[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughter_truth_TotalLength",                  &fdaughter_truth_TotalLength,                 "daughter_truth_TotalLength[NDAUGHTERS]/D");
  fDaughterTree->Branch("daughter_truth_Process",                      &fdaughter_truth_Process,                     "daughter_truth_Process[NDAUGHTERS]/I");
  fDaughterTree->Branch("daughter_truth_EndProcess",                   &fdaughter_truth_EndProcess,                  "daughter_truth_EndProcess[NDAUGHTERS]/I");

  fGrandDaughterTree = tfs->make<TTree>("ReconstructedGrandDaughters", "Reconstructed GrandDaughters");
  fGrandDaughterTree->Branch("NDAUGHTERS",                             &fNDAUGHTERS,                                 "NDAUGHTERS/I");
  fGrandDaughterTree->Branch("NGRANDDAUGHTERS",                        &fNGRANDDAUGHTERS,                            "NGRANDDAUGHTERS/I");
  fGrandDaughterTree->Branch("granddaughterVertex",                    &fgranddaughterVertex,                        "granddaughterVertex[NDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterIstrack",                   &fgranddaughterIstrack,                       "granddaughterIstrack[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughterIsshower",                  &fgranddaughterIsshower,                      "granddaughterIsshower[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughterNHits",                     &fgranddaughterNHits,                         "granddaughterNHits[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughterTheta",                     &fgranddaughterTheta,                         "granddaughterTheta[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterPhi",                       &fgranddaughterPhi,                           "granddaughterPhi[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterLength",                    &fgranddaughterLength,                        "granddaughterLength[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterMomentum",                  &fgranddaughterMomentum,                      "granddaughterMomentum[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterEndMomentum",               &fgranddaughterEndMomentum,                   "granddaughterEndMomentum[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterEndPosition",               &fgranddaughterEndPosition,                   "granddaughterEndPosition[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterStartPosition",             &fgranddaughterStartPosition,                 "granddaughterStartPosition[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterStartDirection",            &fgranddaughterStartDirection,                "granddaughterStartDirection[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterEndDirection",              &fgranddaughterEndDirection,                  "granddaughterEndDirection[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterOpeningAngle",              &fgranddaughterOpeningAngle,                  "granddaughterOpeningAngle[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterShowerBestPlane",           &fgranddaughterShowerBestPlane,               "granddaughterShowerBestPlane[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughterShowerEnergy",              &fgranddaughterShowerEnergy,                  "granddaughterShowerEnergy[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterShowerMIPEnergy",           &fgranddaughterShowerMIPEnergy,               "granddaughterShowerMIPEnergy[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterShowerdEdx",                &fgranddaughterShowerdEdx,                    "granddaughterShowerdEdx[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterMomentumByRangeProton",     &fgranddaughterMomentumByRangeProton,         "granddaughterMomentumByRangeProton[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterMomentumByRangeMuon",       &fgranddaughterMomentumByRangeMuon,           "granddaughterMomentumByRangeMuon[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughterKineticEnergy",             &fgranddaughterKineticEnergy,                 "granddaughterKineticEnergy[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterRange",                     &fgranddaughterRange,                         "granddaughterRange[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterTrkPitchC",                 &fgranddaughterTrkPitchC,                     "granddaughterTrkPitchC[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterID",                        &fgranddaughterID,                            "granddaughterID[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughterMotherID",                  &fgranddaughterMotherID,                      "granddaughterMotherID[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughterT0",                        &fgranddaughterT0,                            "granddaughterT0[NGRANDDAUGHTERS]/D");

  fGrandDaughterTree->Branch("granddaughterPID_Pdg",                   &fgranddaughterPID_Pdg,                       "granddaughterPID_Pdg[NGRANDDAUGHTERS][3]/I");
  fGrandDaughterTree->Branch("granddaughterPID_Ndf",                   &fgranddaughterPID_Ndf,                       "granddaughterPID_Ndf[NGRANDDAUGHTERS][3]/I");
  fGrandDaughterTree->Branch("granddaughterPID_MinChi2",               &fgranddaughterPID_MinChi2,                   "granddaughterPID_MinChi2[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_DeltaChi2",             &fgranddaughterPID_DeltaChi2,                 "granddaughterPID_DeltaChi2[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_Chi2Proton",            &fgranddaughterPID_Chi2Proton,                "granddaughterPID_Chi2Proton[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_Chi2Kaon",              &fgranddaughterPID_Chi2Kaon,                  "granddaughterPID_Chi2Kaon[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_Chi2Pion",              &fgranddaughterPID_Chi2Pion,                  "granddaughterPID_Chi2Pion[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_Chi2Muon",              &fgranddaughterPID_Chi2Muon,                  "granddaughterPID_Chi2Muon[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_MissingE",              &fgranddaughterPID_MissingE,                  "granddaughterPID_MissingE[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_MissingEavg",           &fgranddaughterPID_MissingEavg,               "granddaughterPID_MissingEavg[NGRANDDAUGHTERS][3]/D");
  fGrandDaughterTree->Branch("granddaughterPID_PIDA",                  &fgranddaughterPID_PIDA,                      "granddaughterPID_PIDA[NGRANDDAUGHTERS][3]/D");

  fGrandDaughterTree->Branch("granddaughter_truth_TrackId",            &fgranddaughter_truth_TrackId,                "granddaughter_truth_TrackId[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughter_truth_Pdg",                &fgranddaughter_truth_Pdg,                    "granddaughter_truth_Pdg[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughter_truth_Mother",             &fgranddaughter_truth_Mother,                 "granddaughter_truth_Mother[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughter_truth_StartPosition",      &fgranddaughter_truth_StartPosition,          "granddaughter_truth_StartPosition[NGRANDDAUGHTERS][4]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_EndPosition",        &fgranddaughter_truth_EndPosition,            "granddaughter_truth_EndPosition[NGRANDDAUGHTERS][4]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_Momentum",           &fgranddaughter_truth_Momentum,               "granddaughter_truth_Momentum[NGRANDDAUGHTERS][4]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_EndMomentum",        &fgranddaughter_truth_EndMomentum,            "granddaughter_truth_EndMomentum[NGRANDDAUGHTERS][4]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_P",                  &fgranddaughter_truth_P,                      "granddaughter_truth_P[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_Pt",                 &fgranddaughter_truth_Pt,                     "granddaughter_truth_Pt[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_Mass",               &fgranddaughter_truth_Mass,                   "granddaughter_truth_Mass[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_Theta",              &fgranddaughter_truth_Theta,                  "granddaughter_truth_Theta[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_Phi",                &fgranddaughter_truth_Phi,                    "granddaughter_truth_Phi[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_TotalLength",        &fgranddaughter_truth_TotalLength,            "granddaughter_truth_TotalLength[NGRANDDAUGHTERS]/D");
  fGrandDaughterTree->Branch("granddaughter_truth_Process",            &fgranddaughter_truth_Process,                "granddaughter_truth_Process[NGRANDDAUGHTERS]/I");
  fGrandDaughterTree->Branch("granddaughter_truth_EndProcess",         &fgranddaughter_truth_EndProcess,             "granddaughter_truth_EndProcess[NGRANDDAUGHTERS]/I");

}

// -----------------------------------------------------------------------------
void protoana::ProtoDUNEAnalysisTree::Initialise(){
  
  fRun = -999;
  fSubRun = -999;
  fevent = -999;
  fTimeStamp = -999.0;
  fNpfParticles = 0;
  for(int k=0; k < 5; k++)
    fNactivefembs[k] = -999;

  for(int k=0; k < 3; k++){
    fprimaryVertex[k] = -999.0;
    fdaughterVertex[k] = -999.0;
    fprimaryEndPosition[k] = -999.0;
    fprimaryStartPosition[k] = -999.0;
    fprimaryEndDirection[k] = -999.0;
    fprimaryStartDirection[k] = -999.0;
    fprimaryKineticEnergy[k] = -999.0;
    fprimaryRange[k] = -999.0;
    fprimaryTrkPitchC[k] = -999.0;

    fprimaryPID_Pdg[k] = -999.0;
    fprimaryPID_Ndf[k] = -999.0;
    fprimaryPID_MinChi2[k] = -999.0;
    fprimaryPID_DeltaChi2[k] = -999.0;
    fprimaryPID_Chi2Proton[k] = -999.0;
    fprimaryPID_Chi2Kaon[k] = -999.0;
    fprimaryPID_Chi2Pion[k] = -999.0;
    fprimaryPID_Chi2Muon[k] = -999.0;
    fprimaryPID_MissingE[k] = -999.0;
    fprimaryPID_MissingEavg[k] = -999.0;
    fprimaryPID_PIDA[k] = -999.0;
  }

  fbeamtrack_truth_Origin = -999;
  fbeamtrigger = -999;
  ftof = -999.0;
  fbeamNTracks = -999;
  ftof_expElec = -999.0;
  ftof_expMuon = -999.0;
  ftof_expPion = -999.0;
  ftof_expKaon = -999.0;
  ftof_expProt = -999.0;
  ftof_expDeut = -999.0;
  for(int k=0; k < 2; k++){
    fcerenkovPressure[k] = -999.0;
    fcerenkovTime[k] = -999.0; 
    fcerenkovStatus[k] = -999;
  }
  fbeamtrackMomentum = -999.0;
  fbeamtrackEnergy = 999.0;
  fbeamtrackPdg = -999;
  fbeamtrackID = -999;
  fbeamtrack_truth_Pt = -999.0;
  fbeamtrack_truth_Mass = -999.0;
  fbeamtrack_truth_Theta = -999.0;
  fbeamtrack_truth_Phi = -999.0;
  fbeamtrack_truth_TotalLength = -999.0;
  fbeamtrack_truth_KinEnergy  = -999.0;
  fbeamtrack_truth_P_InTPCActive  = -999.0;
  fbeamtrack_truth_Pt_InTPCActive  = -999.0;
  fbeamtrack_truth_Theta_InTPCActive  = -999.0;
  fbeamtrack_truth_Phi_InTPCActive  = -999.0;
  fbeamtrack_truth_TotalLength_InTPCActive  = -999.0;
  fbeamtrack_truth_KinEnergy_InTPCActive  = -999.0;

  for(int l=0; l < 3; l++){
    fbeamtrackPos[l] = -999.0;
    fbeamtrackDir[l] = -999.0;
  }
  for(int l=0; l < 4; l++){
    fbeamtrack_truth_EndPos[l] = -999.0;
    fbeamtrack_truth_Momentum[l] = -999.0;
    fbeamtrack_truth_Pos_InTPCActive[l] = -999.0;
    fbeamtrack_truth_Momentum_InTPCActive[l] = -999.0;
  }
  fbeamtrack_truth_NDAUGTHERS = 0;
  fbeamtrack_truth_NDECAYDAUGTHERS = 0;
  for(int k=0; k < NMAXTRUTHDAUGTHERS; k++){
    fbeamtrack_truthdaughter_TrackId[k] = -999;
    fbeamtrack_truthdaughter_Pdg[k] = -999;
    fbeamtrack_truthdaughter_Mother[k] = -999;
    fbeamtrack_truthdaughter_P[k] = -999.0;
    fbeamtrack_truthdaughter_Pt[k] = -999.0;
    fbeamtrack_truthdaughter_Mass[k] = -999.0;
    fbeamtrack_truthdaughter_Theta[k] = -999.0;
    fbeamtrack_truthdaughter_Phi[k] = -999.0;
    fbeamtrack_truthdaughter_TotalLength[k] = -999.0;
    fbeamtrack_truthdaughter_Process[k] = -999;
    fbeamtrack_truthdaughter_EndProcess[k] = -999;

    fbeamtrack_truthdecaydaughter_TrackId[k] = -999;
    fbeamtrack_truthdecaydaughter_Pdg[k] = -999;
    fbeamtrack_truthdecaydaughter_Mother[k] = -999;
    fbeamtrack_truthdecaydaughter_P[k] = -999.0;
    fbeamtrack_truthdecaydaughter_Pt[k] = -999.0;
    fbeamtrack_truthdecaydaughter_Mass[k] = -999.0;
    fbeamtrack_truthdecaydaughter_Theta[k] = -999.0;
    fbeamtrack_truthdecaydaughter_Phi[k] = -999.0;
    fbeamtrack_truthdecaydaughter_TotalLength[k] = -999.0;
    fbeamtrack_truthdecaydaughter_Process[k] = -999;
    fbeamtrack_truthdecaydaughter_EndProcess[k] = -999;
    for(int l=0; l < 4; l++){
      fbeamtrack_truthdaughter_StartPosition[k][l] = -999.0;
      fbeamtrack_truthdaughter_EndPosition[k][l] = -999.0;
      fbeamtrack_truthdaughter_Momentum[k][l] = -999.0;
      fbeamtrack_truthdaughter_EndMomentum[k][l] = -999.0;

      fbeamtrack_truthdecaydaughter_StartPosition[k][l] = -999.0;
      fbeamtrack_truthdecaydaughter_EndPosition[k][l] = -999.0;
      fbeamtrack_truthdecaydaughter_Momentum[k][l] = -999.0;
      fbeamtrack_truthdecaydaughter_EndMomentum[k][l] = -999.0;
    }
  }

  fprimaryIstrack = 0;
  fprimaryIsshower = 0;

  fprimaryBDTScore = -999.0;
  fprimaryNHits = -999;
  fprimaryTheta = -999.0;
  fprimaryPhi = -999.0;
  fprimaryLength = -999.0;
  fprimaryMomentum = -999.0;
  fprimaryEndMomentum = -999.0;
  fprimaryOpeningAngle = -999.0;
  fprimaryShowerBestPlane = -999;
  fprimaryID = -999;
  fprimaryMomentumByRangeProton = -999.0;
  fprimaryMomentumByRangeMuon = -999.0;
  fprimaryT0 = -999.0;
  fprimaryEarliestHitPeakTime = -999.0;

  fprimary_truth_Origin = -999;
  fprimary_truth_TrackId = -999;
  fprimary_truth_Pdg = -999;
  fprimary_truth_Mother = -999;
  fprimary_truth_P = -999.0;
  fprimary_truth_Pt = -999.0;
  fprimary_truth_Mass = -999.0;
  fprimary_truth_Theta = -999.0;
  fprimary_truth_Phi = -999.0;
  fprimary_truth_Process = -999;
  fprimary_truth_Isbeammatched = -999;
  fprimary_truth_TotalLength = -999;
  
  for(int l=0; l < 4; l++){
    fprimary_truth_StartPosition[l] = -999.0;
    fprimary_truth_EndPosition[l] = -999.0;
    fprimary_truth_Momentum[l] = -999.0;
    fprimary_truth_EndMomentum[l] = -999.0;
    fprimary_truth_Pos_InTPCActive[l] = -999.0;
    fprimary_truth_Momentum_InTPCActive[l] = -999.0;
  }
  for(int l=0; l < 3; l++){
    fprimaryShowerEnergy[l] = -999;
    fprimaryShowerMIPEnergy[l] = -999;
    fprimaryShowerdEdx[l] = -999;
  }

  fprimary_truth_NDAUGTHERS = 0;
  fprimary_truth_NDECAYDAUGTHERS = 0;
  for(int k=0; k < NMAXTRUTHDAUGTHERS; k++){
    fprimary_truthdaughter_TrackId[k] = -999;
    fprimary_truthdaughter_Pdg[k] = -999;
    fprimary_truthdaughter_Mother[k] = -999;
    fprimary_truthdaughter_P[k] = -999.0;
    fprimary_truthdaughter_Pt[k] = -999.0;
    fprimary_truthdaughter_Mass[k] = -999.0;
    fprimary_truthdaughter_Theta[k] = -999.0;
    fprimary_truthdaughter_Phi[k] = -999.0;
    fprimary_truthdaughter_TotalLength[k] = -999.0;
    fprimary_truthdaughter_Process[k] = -999;
    fprimary_truthdaughter_EndProcess[k] = -999;

    fprimary_truthdecaydaughter_TrackId[k] = -999;
    fprimary_truthdecaydaughter_Pdg[k] = -999;
    fprimary_truthdecaydaughter_Mother[k] = -999;
    fprimary_truthdecaydaughter_P[k] = -999.0;
    fprimary_truthdecaydaughter_Pt[k] = -999.0;
    fprimary_truthdecaydaughter_Mass[k] = -999.0;
    fprimary_truthdecaydaughter_Theta[k] = -999.0;
    fprimary_truthdecaydaughter_Phi[k] = -999.0;
    fprimary_truthdecaydaughter_TotalLength[k] = -999.0;
    fprimary_truthdecaydaughter_Process[k] = -999;
    fprimary_truthdecaydaughter_EndProcess[k] = -999;
    for(int l=0; l < 4; l++){
      fprimary_truthdaughter_StartPosition[k][l] = -999.0;
      fprimary_truthdaughter_EndPosition[k][l] = -999.0;
      fprimary_truthdaughter_Momentum[k][l] = -999.0;
      fprimary_truthdaughter_EndMomentum[k][l] = -999.0;

      fprimary_truthdecaydaughter_StartPosition[k][l] = -999.0;
      fprimary_truthdecaydaughter_EndPosition[k][l] = -999.0;
      fprimary_truthdecaydaughter_Momentum[k][l] = -999.0;
      fprimary_truthdecaydaughter_EndMomentum[k][l] = -999.0;
    }
  }

  fNDAUGHTERS = 0;
  fNGRANDDAUGHTERS = 0;
  for(int k=0; k < NMAXDAUGTHERS; k++){
    fdaughterIstrack[k] = -999;
    fdaughterIsshower[k] = -999;
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
      fdaughterTrkPitchC[k][l] = -999.0;

      fdaughterPID_Pdg[k][l] = -999.0;
      fdaughterPID_Ndf[k][l] = -999.0;
      fdaughterPID_MinChi2[k][l] = -999.0;
      fdaughterPID_DeltaChi2[k][l] = -999.0;
      fdaughterPID_Chi2Proton[k][l] = -999.0;
      fdaughterPID_Chi2Kaon[k][l] = -999.0;
      fdaughterPID_Chi2Pion[k][l] = -999.0;
      fdaughterPID_Chi2Muon[k][l] = -999.0;
      fdaughterPID_MissingE[k][l] = -999.0;
      fdaughterPID_MissingEavg[k][l] = -999.0;
      fdaughterPID_PIDA[k][l] = -999.0;
      fgranddaughterVertex[k][l] = -999.0;

      fdaughterShowerEnergy[k][l] = -999;
      fdaughterShowerMIPEnergy[k][l] = -999;
      fdaughterShowerdEdx[k][l] = -999;
    }
    fdaughterOpeningAngle[k] = -999.0;
    fdaughterShowerBestPlane[k] = -999;
    fdaughterMomentumByRangeProton[k] = -999.0;
    fdaughterMomentumByRangeMuon[k] = -999.0;
    fdaughterID[k] = -999;
    fdaughterT0[k] = -999;

    fdaughter_truth_TrackId[k] = -999;
    fdaughter_truth_Pdg[k] = -999;
    fdaughter_truth_Mother[k] = -999;
    fdaughter_truth_P[k] = -999.0;
    fdaughter_truth_Pt[k] = -999.0;
    fdaughter_truth_Mass[k] = -999.0;
    fdaughter_truth_Theta[k] = -999.0;
    fdaughter_truth_Phi[k] = -999.0;
    fdaughter_truth_Process[k] = -999;
    fdaughter_truth_TotalLength[k] = -999.0;
    for(int l=0; l < 4; l++){
      fdaughter_truth_StartPosition[k][l] = -999.0;
      fdaughter_truth_EndPosition[k][l] = -999.0;
      fdaughter_truth_Momentum[k][l] = -999.0;
      fdaughter_truth_EndMomentum[k][l] = -999.0;
    }

    // Grand-daugthers
    fgranddaughterIstrack[k] = -999;
    fgranddaughterIsshower[k] = -999;
    fgranddaughterNHits[k] = -999;
    fgranddaughterTheta[k] = -999.0;
    fgranddaughterPhi[k] = -999.0;
    fgranddaughterLength[k] = -999.0;
    fgranddaughterMomentum[k] = -999.0;
    fgranddaughterEndMomentum[k] = -999.0;
    for(int l=0; l < 3; l++){
      fgranddaughterEndPosition[k][l] = -999.0;
      fgranddaughterStartPosition[k][l] = -999.0;
      fgranddaughterEndDirection[k][l] = -999.0;
      fgranddaughterStartDirection[k][l] = -999.0;
      fgranddaughterKineticEnergy[k][l] = -999.0;
      fgranddaughterRange[k][l] = -999.0;
      fgranddaughterTrkPitchC[k][l] = -999.0;

      fgranddaughterPID_Pdg[k][l] = -999.0;
      fgranddaughterPID_Ndf[k][l] = -999.0;
      fgranddaughterPID_MinChi2[k][l] = -999.0;
      fgranddaughterPID_DeltaChi2[k][l] = -999.0;
      fgranddaughterPID_Chi2Proton[k][l] = -999.0;
      fgranddaughterPID_Chi2Kaon[k][l] = -999.0;
      fgranddaughterPID_Chi2Pion[k][l] = -999.0;
      fgranddaughterPID_Chi2Muon[k][l] = -999.0;
      fgranddaughterPID_MissingE[k][l] = -999.0;
      fgranddaughterPID_MissingEavg[k][l] = -999.0;
      fgranddaughterPID_PIDA[k][l] = -999.0;

      fgranddaughterShowerEnergy[k][l] = -999;
      fgranddaughterShowerMIPEnergy[k][l] = -999;
      fgranddaughterShowerdEdx[k][l] = -999;
    }
    fgranddaughterOpeningAngle[k] = -999.0;
    fgranddaughterShowerBestPlane[k] = -999;
    fgranddaughterMomentumByRangeProton[k] = -999.0;
    fgranddaughterMomentumByRangeMuon[k] = -999.0;
    fgranddaughterID[k] = -999;
    fgranddaughterMotherID[k] = -999;
    fgranddaughterT0[k] = -999;

    fgranddaughter_truth_TrackId[k] = -999;
    fgranddaughter_truth_Pdg[k] = -999;
    fgranddaughter_truth_Mother[k] = -999;
    fgranddaughter_truth_P[k] = -999.0;
    fgranddaughter_truth_Pt[k] = -999.0;
    fgranddaughter_truth_Mass[k] = -999.0;
    fgranddaughter_truth_Theta[k] = -999.0;
    fgranddaughter_truth_Phi[k] = -999.0;
    fgranddaughter_truth_Process[k] = -999;
    fgranddaughter_truth_TotalLength[k] = -999.0;
    for(int l=0; l < 4; l++){
      fgranddaughter_truth_StartPosition[k][l] = -999.0;
      fgranddaughter_truth_EndPosition[k][l] = -999.0;
      fgranddaughter_truth_Momentum[k][l] = -999.0;
      fgranddaughter_truth_EndMomentum[k][l] = -999.0;
    }

  }


}

DEFINE_ART_MODULE(protoana::ProtoDUNEAnalysisTree)
