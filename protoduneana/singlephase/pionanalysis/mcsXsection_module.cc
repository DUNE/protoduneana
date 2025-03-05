////////////////////////////////////////////////////////////////////////
// Class:       mcsXsection                                              //
// File:        mcsXsection_module.cc                                    //
// Ajib Paudel MCS cross section module                               //   
// ajib.paudel@phys.ksu.edu                                           //
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "larcorealg/Geometry/Exceptions.h"
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
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
//#include "dunesim/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
//#include "dunesim/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNETrackUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEPFParticleUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "larsim/MCCheater/ParticleInventoryService.h"


//#g4reweight
#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "geant4reweight/src/ReweightBase/G4ReweighterFactory.hh"
#include "geant4reweight/src/ReweightBase/G4Reweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/PropBase/G4ReweightParameterMaker.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "protoduneana/Utilities/G4ReweightUtils.h"


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
#include <tuple>

// Maximum number of beam particles to save
const int NMAXDAUGTHERS = 30;
double  m_pi=0.1395;//GeV/c^2
double prim_energy=0.0;
//double geant_energy=0;


float soln(float x1, float x2, float x3, float x4, float y1, float y2, float y3, float y4, float result[2]){
  if(x3==x4){
    result[0]=x3;
    result[1]=y3;
  }
  else if((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)==0){
    result[0]=x3;
    result[1]=y3;
  }
  else{
    result[0]=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
    result[1]=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
  }
  return *result;
}
double angle2d(double x0, double y0, double x1, double y1, double x2, double y2){
  double theta=((x1-x0)*(x2-x1)+(y1-y0)*(y2-y1))/sqrt(((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))*((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
  return acos(theta);

} 


double theta12(double x1, double x2, double y1, double y2, double z1, double z2,double x1p, double x2p, double y1p, double y2p, double z1p, double z2p){
  double numer=(x2-x1)*(x2p-x1p)+(y2-y1)*(y2p-y1p)+(z2-z1)*(z2p-z1p);
  double den1=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
  double den2=(x2p-x1p)*(x2p-x1p)+(y2p-y1p)*(y2p-y1p)+(z2p-z1p)*(z2p-z1p);
  return 180/3.14*acos(numer/sqrt(den1*den2));
} 


namespace protoana {
  class mcsXsection;
}

using protoana::G4ReweightUtils::CreateRWTraj;
using protoana::G4ReweightUtils::CreateNRWTrajs;
using protoana::G4ReweightUtils::GetNTrajPMSigmaWeights;
using protoana::G4ReweightUtils::GetNTrajWeightFromSetPars;

class protoana::mcsXsection : public art::EDAnalyzer {
public:

  explicit mcsXsection(fhicl::ParameterSet const & p);

  mcsXsection(mcsXsection const &) = delete;
  mcsXsection(mcsXsection &&) = delete;
  mcsXsection & operator = (mcsXsection const &) = delete;
  mcsXsection & operator = (mcsXsection &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & evt) override;

private:
  const art::InputTag fTrackModuleLabel;
  const art::InputTag fBeamModuleLabel;
  // Helper utility functions
  protoana::ProtoDUNEBeamCuts beam_cuts;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNEDataUtils dataUtil;

  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
 

  // Initialise tree variables
  void Initialise();

  // Fill cosmics tree
  void FillCosmicsTree(art::Event const & evt, std::string pfParticleTag);

  // fcl parameters
  //const art::InputTag fNNetModuleLabel; //label of the module used for CNN tagging
 
 
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  //bool fVerbose;

  geo::GeometryCore const * fGeometry;
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
  std::vector<int> beamtrk_z_wire;
  std::vector<int> beamtrk_z_tpc;

  std::vector<float> beamtrk_Px;
  std::vector<float> beamtrk_Py;
  std::vector<float> beamtrk_Pz;
  std::vector<float> beamtrk_Eng;
  std::vector<float> peakT_2;
  std::vector<float> hitz_2;
  // std::vector<float> inelscore;
  //std::vector<float> elscore;
  //std::vector<float> nonescore;
  
  /////////
  //beam info from spectrometer
  std::vector<double> beamMomentum_spec;
  std::vector<double> beamPosx_spec;
  std::vector<double> beamPosy_spec;
  std::vector<double> beamPosz_spec;
  
  std::vector<double> beamDirx_spec;
  std::vector<double> beamDiry_spec;
  std::vector<double> beamDirz_spec;






  std::vector<float> int_2;
  std::vector<float> hitz_1;
  std::vector<float> int_1;
  std::vector<float> hitz_0;
  std::vector<float> int_0;

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
  // double fprimary_truth_Startwiretime[4];
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
  std::vector<double> interactionU;
  std::vector<double> interactionV;
  std::vector<double> interactionW;
  std::vector<double> interactionT;



  std::vector<std::string> interactionProcesslist;
  std::vector<double> interactionAngles;
  std::vector<double> interactionAnglesUT;
  std::vector<double> interactionAnglesVT;
  std::vector<double> interactionAnglesZT;
  std::vector<double> Zintersection;
  std::vector<double> Zintersection1;
  std::vector<double> timeintersection;
  std::vector<double> timeintersection1;
  std::vector<double> deltaZint;
  std::vector<double> deltatimeint;
  //calo info
  std::vector< std::vector<double> > primtrk_dqdx;
  std::vector< std::vector<double> > primtrk_resrange;
  std::vector< std::vector<double> > primtrk_dedx;
  std::vector<double> primtrk_range;
  std::vector< std::vector<double> > primtrk_hitx;
  std::vector< std::vector<double> > primtrk_hity;
  std::vector< std::vector<double> > primtrk_hitz;
  std::vector< std::vector<int> > primtrk_hitz_wire;
  std::vector< std::vector<int> > primtrk_hitz_tpc;

  std::vector< std::vector<double> > primtrk_pitch;
  std::vector<std::vector<double> > primtrk_truth_Z;
  std::vector<std::vector<int> > primtrk_truth_Z_wire;
  std::vector<std::vector<int> > primtrk_truth_Z_tpc;


  std::vector<std::vector<double> > primtrk_truth_Eng;
  std::vector<std::vector<double> > primtrk_truth_trkide;
  std::vector<std::vector<float> > wireno_0;
  std::vector<std::vector<float> > wireno_1;
  std::vector<std::vector<float> > wireno_2;
  std::vector<std::vector<float> > peakTime_0;
  std::vector<std::vector<float> > peakTime_1;
  std::vector<std::vector<float> > peakTime_2;
  std::vector<std::vector<float> > dq_0;
  std::vector<std::vector<float> > dq_1;
  std::vector<std::vector<float> > dq_2;

  ////Michel tagging info

  std::vector< std::vector<int> > endhitssecondary;
  std::vector< std::vector<double> > secondarystartx;
  std::vector< std::vector<double> > secondaryendx;
  std::vector< std::vector<double> > secondarystarty;
  std::vector< std::vector<double> > secondaryendy;
  std::vector< std::vector<double> > secondarystartz;
  std::vector< std::vector<double> > secondaryendz;
  std::vector< std::vector<double> > dQmichel;
  std::vector< std::vector<double> > dQtrackend;
  std::vector< std::vector<double> > dQtrackbegin;
  std::vector< std::vector<double> > primsectheta;
  std::vector< std::vector<double> > tracklengthsecondary;
  std::vector< std::vector<int> > MtrackID;
  ////CNN for michel tagging
  std::vector< std::vector<double> > trackscore;
  std::vector< std::vector<double> > emscore;
  std::vector< std::vector<double> > michelscore;
  std::vector< std::vector<double> > nonescore;


 
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


  int true_beam_PDG;
  int true_beam_ID;

  int n_g4rw;
  double g4rw_react_st;
  double g4rw_d_react;
  double g4rw_elast_st;
  double g4rw_d_elast;
  std::vector<double> g4rw_react; 
  std::vector<double> g4rw_elast; 
  std::vector<double> g4rw_set_weights;
  //Geant4Reweight stuff
  int RW_PDG;
  TFile FracsFile, XSecFile;
  std::vector<fhicl::ParameterSet> ParSet;
  G4ReweightParameterMaker ParMaker;
  G4ReweightManager RWManager;
  G4MultiReweighter MultiRW;
  //G4ReweighterFactory RWFactory;
  //G4Reweighter * theRW;

  //bool fVerbose;


};


protoana::mcsXsection::mcsXsection(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
 
  // fNNetModuleLabel(p.get<art::InputTag>("NNetModuleLabel")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  beam_cuts(p.get<fhicl::ParameterSet>("BeamCuts")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils")),
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  RW_PDG(p.get<int>("RW_PDG")),
  FracsFile( (p.get< std::string >( "FracsFile" )).c_str(), "OPEN" ),
  XSecFile( (p.get< std::string >( "XSecFile" )).c_str(), "OPEN"),
  ParSet(p.get<std::vector<fhicl::ParameterSet>>("ParameterSet")),
  ParMaker(ParSet, RW_PDG),
  RWManager({p.get<fhicl::ParameterSet>("Material")}),
  MultiRW(RW_PDG, FracsFile, ParSet,
	p.get<fhicl::ParameterSet>("Material"),
	&RWManager) {
  std::cout << "done" << std::endl;
}

/* fVerbose(p.get<bool>("Verbose"))
     {

     }
  */
  void protoana::mcsXsection::beginJob(){

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
    fPandoraBeam->Branch("beamtrk_z_wire",&beamtrk_z_wire);
    fPandoraBeam->Branch("beamtrk_z_tpc",&beamtrk_z_tpc);

    fPandoraBeam->Branch("beamtrk_Px",&beamtrk_Px);
    fPandoraBeam->Branch("beamtrk_Py",&beamtrk_Py);
    fPandoraBeam->Branch("beamtrk_Pz",&beamtrk_Pz);
    fPandoraBeam->Branch("beamtrk_Eng",&beamtrk_Eng);
    ////.............spectrometer parameters
    fPandoraBeam->Branch("beamMomentum_spec",&beamMomentum_spec);
    fPandoraBeam->Branch("beamPosx_spec",&beamPosx_spec);
    fPandoraBeam->Branch("beamPosy_spec",&beamPosy_spec);
    fPandoraBeam->Branch("beamPosz_spec",&beamPosz_spec);
    fPandoraBeam->Branch("beamDirx_spec",&beamDirx_spec);
    fPandoraBeam->Branch("beamDiry_spec",&beamDiry_spec);
    fPandoraBeam->Branch("beamDirz_spec",&beamDirz_spec);




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
    // fPandoraBeam->Branch("primary_truth_Startwiretime",   &fprimary_truth_Startwiretime,   "primary_truth_Startwiretime[4]/D");
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
    fPandoraBeam->Branch("interactionU",&interactionU);
    fPandoraBeam->Branch("interactionV",&interactionV);
    fPandoraBeam->Branch("interactionW",&interactionW);
    fPandoraBeam->Branch("interactionT",&interactionT);
    fPandoraBeam->Branch("interactionProcesslist",&interactionProcesslist);
    fPandoraBeam->Branch("interactionAngles",&interactionAngles);
    fPandoraBeam->Branch("interactionAnglesUT",&interactionAnglesUT);
    fPandoraBeam->Branch("interactionAnglesVT",&interactionAnglesVT);
    fPandoraBeam->Branch("interactionAnglesZT",&interactionAnglesZT);

    fPandoraBeam->Branch("Zintersection",&Zintersection);
    fPandoraBeam->Branch("timeintersection",&timeintersection);
    fPandoraBeam->Branch("Zintersection1",&Zintersection1);
    fPandoraBeam->Branch("timeintersection1",&timeintersection1);
    fPandoraBeam->Branch("deltaZint",&deltaZint);
    fPandoraBeam->Branch("deltatimeint",&deltatimeint);
    fPandoraBeam->Branch("peakT_2",&peakT_2);
    fPandoraBeam->Branch("hitz_2",&hitz_2);
    // fPandoraBeam->Branch("inelscore",&inelscore);
    //fPandoraBeam->Branch("elscore",&elscore);
    //fPandoraBeam->Branch("nonescore",&nonescore);

    fPandoraBeam->Branch("int_2",&int_2);
    fPandoraBeam->Branch("hitz_1",&hitz_1);
    fPandoraBeam->Branch("int_1",&int_1);
    fPandoraBeam->Branch("hitz_0",&hitz_0);
    fPandoraBeam->Branch("int_0",&int_0);
    fPandoraBeam->Branch("wireno_0",&wireno_0);
    fPandoraBeam->Branch("wireno_1",&wireno_1);
    fPandoraBeam->Branch("wireno_2",&wireno_2);
    fPandoraBeam->Branch("peakTime_0",&peakTime_0);
    fPandoraBeam->Branch("peakTime_1",&peakTime_1);
    fPandoraBeam->Branch("peakTime_2",&peakTime_2);
    fPandoraBeam->Branch("dq_0",&dq_0);
    fPandoraBeam->Branch("dq_1",&dq_1);
    fPandoraBeam->Branch("dq_2",&dq_2);
    //Michel tagging
    fPandoraBeam->Branch("endhitssecondary",&endhitssecondary);
    fPandoraBeam->Branch("secondarystartx",&secondarystartx);
    fPandoraBeam->Branch("secondarystarty",&secondarystarty);
    fPandoraBeam->Branch("secondarystartz",&secondarystartz);
    fPandoraBeam->Branch("secondaryendx",&secondaryendx);
    fPandoraBeam->Branch("secondaryendy",&secondaryendy);
    fPandoraBeam->Branch("secondaryendz",&secondaryendz);
    fPandoraBeam->Branch("dQmichel",&dQmichel);
    fPandoraBeam->Branch("primsectheta",&primsectheta);
    fPandoraBeam->Branch("dQtrackbegin",&dQtrackbegin);
    fPandoraBeam->Branch("dQtrackend",&dQtrackend);
    fPandoraBeam->Branch("tracklengthsecondary",&tracklengthsecondary);
    fPandoraBeam->Branch("MtrackID",&MtrackID);
    fPandoraBeam->Branch("trackscore",&trackscore);
    fPandoraBeam->Branch("emscore",&emscore);
    fPandoraBeam->Branch("michelscore",&michelscore);
    fPandoraBeam->Branch("nonescore",&nonescore);
    /////////



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
    fPandoraBeam->Branch("primtrk_hitz_wire",&primtrk_hitz_wire);
    fPandoraBeam->Branch("primtrk_hitz_tpc",&primtrk_hitz_tpc);
    fPandoraBeam->Branch("primtrk_pitch",&primtrk_pitch);

    ////g4reweight
    fPandoraBeam->Branch("n_g4rw",  &n_g4rw, "n_g4rw/I");
    fPandoraBeam->Branch("g4rw_react_st",  &g4rw_react_st, "g4rw_react_st/D");
    fPandoraBeam->Branch("g4rw_d_react",  &g4rw_d_react, "g4rw_d_react/D");
    fPandoraBeam->Branch("g4rw_elast_st",  &g4rw_elast_st, "g4rw_elast_st/D");
    fPandoraBeam->Branch("g4rw_d_elast",  &g4rw_d_elast, "g4rw_d_elast/D");
    fPandoraBeam->Branch("g4rw_react", &g4rw_react);
    fPandoraBeam->Branch("g4rw_elast", &g4rw_elast);
    fPandoraBeam->Branch("g4rw_set_weights", &g4rw_set_weights);



    ///////////////////////////////////////////////////
    fPandoraBeam->Branch("primtrk_truth_Z",&primtrk_truth_Z);//primary track true Z positions 
    fPandoraBeam->Branch("primtrk_truth_Z_wire",&primtrk_truth_Z_wire);//primary track true Z positions 
    fPandoraBeam->Branch("primtrk_truth_Z_tpc",&primtrk_truth_Z_tpc);//primary track true Z positions 


    fPandoraBeam->Branch("primtrk_truth_Eng",&primtrk_truth_Eng);//primary track true Energy deposited for each Z position
    fPandoraBeam->Branch("primtrk_truth_trkide",&primtrk_truth_trkide);//primary track true Energy deposited for each Z position
  }//begin job 

  void protoana::mcsXsection::analyze(art::Event const & evt){

    // Initialise tree parameters
    Initialise();

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    fGeometry = art::ServiceHandle<geo::Geometry>().get();
    // const sim::ParticleList& plist=pi_serv->ParticleList();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    // anab::MVAReader<recob::Hit,3> hitResults(evt, fNNetModuleLabel);
    anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel");

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
    std::cout<<"hhecking code "<<std::endl;
    bool beamTriggerEvent = false;
    // If this event is MC then we can check what the true beam particle is
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    if(!evt.isRealData()){
      //Spectrometer values
      //for prod. 3, new implementation to access the beam momentum from spectrometer //////////////////////////////
      auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("generator");
      std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
      if( beamHandle.isValid()){
	art::fill_ptr_vector(beamVec, beamHandle);
      }
      const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one

      //Access momentum
      const std::vector< double > & momenta = beamEvent.GetRecoBeamMomenta();
      std::cout << "Number of reconstructed beam momenta from spec: " << momenta.size() << std::endl;

      if( momenta.size() > 0 ) std::cout << "Measured Beam Momentum from spec: " << momenta.at(0) << std::endl;
		
      //std::cout<<"beam mom size:"<<momenta.size()<<std::endl;
      for (size_t i = 0; i<momenta.size(); ++i){
	beamMomentum_spec.push_back(momenta[i]);
	//std::cout<<"beam mom["<<i<<"]:"<<momenta[i]<<" [GeV]"<<std::endl;
      }

      auto & btracks = beamEvent.GetBeamTracks();
      std::cout<<"beam trk size:"<<btracks.size()<<std::endl;
      for (size_t i = 0; i<btracks.size(); ++i){
	std::cout<<"beamPosx/beamPosy/beamPosz:"<<btracks[i].End().X()<<"/"<<btracks[i].End().Y()<<"/"<<btracks[i].End().Z()<<std::endl;
	std::cout<<"beamDirx/beamDiry/beamDirz:"<<btracks[i].StartDirection().X()<<"/"<<btracks[i].StartDirection().Y()<<"/"<<btracks[i].StartDirection().Z()<<std::endl;

	beamPosx_spec.push_back(btracks[i].End().X());
	beamPosy_spec.push_back(btracks[i].End().Y());
	beamPosz_spec.push_back(btracks[i].End().Z());
	beamDirx_spec.push_back(btracks[i].StartDirection().X());
	beamDiry_spec.push_back(btracks[i].StartDirection().Y());
	beamDirz_spec.push_back(btracks[i].StartDirection().Z());

      }

      //////Spectrometer values
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
	std::vector<int> tmp_wire;
	std::vector<int> tmp_tpc;      

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
	beamid             = geantGoodParticle->TrackId();

	///////////////////////////////////////
	//g4 reweight -------------------------------------------------------------------------------------------------//
	true_beam_PDG = geantGoodParticle->PdgCode();
	true_beam_ID = geantGoodParticle->TrackId();
	//true_beam_len = geantGoodParticle->Trajectory().TotalLength();

	std::cout << "Doing reweight" << std::endl;
	if (true_beam_PDG == RW_PDG) { //if PDG=RW_PDG
	  std::vector<G4ReweightTraj *> trajs = CreateNRWTrajs(*geantGoodParticle, pi_serv->ParticleList(),*fGeometry, fevent, "LAr", true);
	  //G4ReweightTraj theTraj(true_beam_ID, true_beam_PDG, 0, fevent, {0,0});
	  int n_rw_array=200;
	  double react_st=0.0;
	  double d_react=.01;
	  double elast_st=.9;
	  double d_elast=.05;
	  n_g4rw=n_rw_array;
	  g4rw_react_st=react_st; g4rw_d_react=d_react;
	  g4rw_elast_st=elast_st; g4rw_d_elast=d_elast;
	  std::cout<<"G4reweighting loop step 1"<<std::endl;
	  if (ParSet.size()==2) { //if parset size=2 from fichl setting, i.e. dim(el+inel)=2 
	    for (int i = 0; i < n_rw_array; ++i) {
	      for (int j = 0; j < n_rw_array; ++j) {
		std::cout<<"loop i, j "<<i<<"  "<<j<<std::endl;
		double tmp_p1=(react_st + i*d_react); //reac
		double tmp_p2=(elast_st + j*d_elast); //elast
		std::vector<double> input_values = {tmp_p1,tmp_p2}; //HY:the index here based on the sequence of the fichl setting
		//	std::vector<double> input_values = {ParSet.size(), 1.};
		std::cout<<"loop temp1, 2 "<<tmp_p1<<"  "<<tmp_p2<<std::endl;			
		bool set_values = MultiRW.SetAllParameterValues(input_values);
		std::cout<<"setvalues "<<set_values<<std::endl;
		if (!set_values) continue;
		double temp_w = GetNTrajWeightFromSetPars(trajs, MultiRW);
		//double temp_w = MultiRW.GetWeightFromSetParameters(theTraj);
		std::cout<<"weight "<<temp_w<<std::endl;
		g4rw_react.push_back(tmp_p1);
		g4rw_elast.push_back(tmp_p2);
		g4rw_set_weights.push_back(temp_w);
		//std::cout<<"(inel,el,w):("<<tmp_p1<<","<<tmp_p2<<","<<temp_w<<")"<<std::endl;
	      }
	    }
	  } //if parset size=2
	} //if PDG=RW_PDG
	//g4 reweight -------------------------------------------------------------------------------------------------//


	////////////////////////////////////////

	///

	truth_last_process=geantGoodParticle->EndProcess();
	////here is another new parameter added

	fprimary_truth_Process = int(geantGoodParticle->Trajectory().ProcessToKey(geantGoodParticle->Process()));
	prim_energy=0;
	for(size_t i_s=0; i_s < geantGoodParticle->NumberTrajectoryPoints(); i_s++){ //loop over beam tracks
	  //	if(geantGoodParticle->Position(i_s).Z()>0) break;
          auto const pos_true = geo::vect::toPoint(geantGoodParticle->Position(i_s).Vect());
          beamtrk_x.push_back(pos_true.X());
          beamtrk_y.push_back(pos_true.Y());
          beamtrk_z.push_back(pos_true.Z());


	  /*	if(geantGoodParticle->Position(i_s).Z()>=0 && geantGoodParticle->Position(i_s).Z()<=230){
		geo::WireID wireID = fWireReadoutGeom.NearestWireID(pos_true, 2);
		if (!wireID) wireID = fWireReadoutGeom.Plane(2).ClosestWireID(wireID);
		beamtrk_z_wire.push_back(wireID.Wire);
		beamtrk_z_tpc.push_back(wireID.TPC);
		}
		if(geantGoodParticle->Position(i_s).Z()<0 || geantGoodParticle->Position(i_s).Z()>230){
		beamtrk_z_wire.push_back(-99999);
		beamtrk_z_tpc.push_back(-99999);

		}*/
	  geo::TPCID tpc = fGeometry->FindTPCAtPosition(pos_true);
	  if(tpc.isValid){
	    int tpc_no=tpc.TPC;
	    geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
	    geo::WireID wireID;
	    //	geo::WireID wireID = fWireReadoutGeom.NearestWireID(pos_true, planeID);
	    try{
              wireID = fWireReadoutGeom.Plane(planeID).NearestWireID(pos_true);
	    }
	    catch(geo::InvalidWireError const& e) {
	      wireID = e.suggestedWireID(); // pick the closest valid wire
	    }
	    beamtrk_z_wire.push_back(wireID.Wire);
	    beamtrk_z_tpc.push_back(wireID.TPC);
	  }
	  if(!tpc.isValid){
	    beamtrk_z_wire.push_back(-9999);
	    beamtrk_z_tpc.push_back(-9999);

	  }


	


	  /*	geo::PlaneGeo const& plane = fWireReadoutGeom.Plane(2);
		geo::WireID wireID;
		try {
		wireID = plane.NearestWireID(pos_true);
		}
		catch (geo::InvalidWireIDError const& e) {
		if (!e.hasSuggestedWire()) throw;
		wireID = plane.ClosestWireID(e.suggestedWireID());
		}
		beamtrk_z_wire.push_back(wireID.Wire);
		beamtrk_z_tpc.push_back(wireID.TPC);*/
	  ///////trying new stuff
        


	  beamtrk_Px.push_back(geantGoodParticle->Momentum(i_s).X());
	  beamtrk_Py.push_back(geantGoodParticle->Momentum(i_s).Y());
	  beamtrk_Pz.push_back(geantGoodParticle->Momentum(i_s).Z());
	
	  beamtrk_Eng.push_back(geantGoodParticle->Momentum(i_s).E()-geantGoodParticle->Mass());
	  if(geantGoodParticle->Position(i_s).Z()<0) prim_energy=1000*(geantGoodParticle->Momentum(i_s).E()-geantGoodParticle->Mass());
	  if(geantGoodParticle->Position(i_s).Z()<0) fbeamtrackEnergy=prim_energy; //correct energy at the beginning of track
	} //loop over beam trks

	//new section 
	art::ServiceHandle<cheat::BackTrackerService> bt_serv;


	// art::ServiceHandle<geo::Geometry> geom;
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
	    double xval=((truetraj.at(couple1.first)).first).X();
	    double zval=((truetraj.at(couple1.first)).first).Z();
	    unsigned int tpcno=1;
	    if(xval<=0 && zval<232) tpcno=1;
	    if(xval<=0 && zval>232 && zval<464) tpcno=5; 
	    if(xval<=0 && zval>=464) tpcno=9;
	    if(xval>0 && zval<232) tpcno=2; 
	    if(xval>0 && zval>232 && zval<464) tpcno=6; 
	    if(xval>0 && zval>=464) tpcno=10;

            auto const pos = geo::vect::toPoint(truetraj.at(couple1.first).first.Vect());
            geo::TPCID const tpcID{0, tpcno};
            interactionT.push_back(detProp.ConvertXToTicks(pos.X(), geo::PlaneID{tpcID, 2}));
            interactionU.push_back(fWireReadoutGeom.Plane(geo::PlaneID{tpcID, 0}).WireCoordinate(pos));
            interactionV.push_back(fWireReadoutGeom.Plane(geo::PlaneID{tpcID, 1}).WireCoordinate(pos));
            interactionW.push_back(fWireReadoutGeom.Plane(geo::PlaneID{tpcID, 2}).WireCoordinate(pos));
	    interactionProcesslist.push_back(truetraj.KeyToProcess(couple1.second));
	    std::cout<<"number of interactions "<<thisTrajectoryProcessMap1.size()<<std::endl;
            std::cout<<"int X, Y, Z and process "<<pos.X()<<" "<<pos.Y()<<" "<<pos.Z()<<" "<<truetraj.KeyToProcess(couple1.second)<<std::endl;
	    ///get the interaction angle here
	    double interactionAngle = 999999.; // This needs to be changed
	    double ut=999999;
	    double vt=999999;
	    double zt=999999;
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
	      std::cout<<"distance between points values "<<distanceBtwPoint.X()<<" "<<distanceBtwPoint.Z()<<" "<<distanceBtwPointNext.X()<<" "<<distanceBtwPointNext.Z()<<std::endl;

	      interactionAngles.push_back(TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  ));

	   
              geo::TPCID const tpc_0_1{0, 1};
              geo::PlaneID const plane_0{tpc_0_1, 0};
              geo::PlaneID const plane_1{tpc_0_1, 1};
              using geo::vect::toPoint;
              double u0=0.4669*fWireReadoutGeom.Plane(plane_0).WireCoordinate(toPoint(prevInteractionPos3D));
              double v0=0.4669*fWireReadoutGeom.Plane(plane_1).WireCoordinate(toPoint(prevInteractionPos3D));
	      double z0=0.4792*prevInteractionPos3D.Z();
	      double x0=prevInteractionPos3D.X();

              double u1=0.4669*fWireReadoutGeom.Plane(plane_0).WireCoordinate(toPoint(interactionPos3D));
              double v1=0.4669*fWireReadoutGeom.Plane(plane_1).WireCoordinate(toPoint(interactionPos3D));
	      double z1=0.4792*interactionPos3D.Z();
	      double x1=interactionPos3D.Z();

              double u2=0.4669*fWireReadoutGeom.Plane(plane_0).WireCoordinate(toPoint(nextInteractionPos3D));
              double v2=0.4669*fWireReadoutGeom.Plane(plane_1).WireCoordinate(toPoint(nextInteractionPos3D));
	      double z2=0.4792*nextInteractionPos3D.Z();
	      double x2=nextInteractionPos3D.Z();

	      interactionAnglesUT.push_back(angle2d(u0,x0,u1,x1,u2,x2));
	      interactionAnglesVT.push_back(angle2d(v0,x0,v1,x1,v2,x2));
	      interactionAnglesZT.push_back(angle2d(z0,x0,z1,x1,z2,x2));
 

	      continue;
	    }
	    interactionAnglesUT.push_back(ut);
	    interactionAnglesVT.push_back(vt);
	    interactionAnglesZT.push_back(zt);
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
          geo::Point_t const pos_true{currentIde.x,currentIde.y,currentIde.z};
	  /*if(currentIde.z>=0 && currentIde.z<=230){
	    geo::WireID wireID = fWireReadoutGeom.NearestWireID(pos_true, 2);
	    if (!wireID) wireID = fWireReadoutGeom.Plane(2).ClosestWireID(wireID);
	    tmp_wire.push_back(wireID.Wire);
	    tmp_tpc.push_back(wireID.TPC);
	    }
	    if(currentIde.z<0 && currentIde.z>230){
	    tmp_wire.push_back(-99999);
	    tmp_tpc.push_back(-99999);
	    }*/
	  geo::TPCID tpc1 = fGeometry->FindTPCAtPosition(pos_true);
	  if(tpc1.isValid){
	    int tpc_no=tpc1.TPC;
	    geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
	    //	geo::WireID wireID = fWireReadoutGeom.NearestWireID(pos_true, planeID);
	    geo::WireID wireID;
	    try{
              wireID = fWireReadoutGeom.Plane(planeID).NearestWireID(pos_true);
	    }
	    catch(geo::InvalidWireError const& e) {
	      wireID = e.suggestedWireID(); // pick the closest valid wire
	    }
	    tmp_wire.push_back(wireID.Wire);
	    tmp_tpc.push_back(wireID.TPC);
	  }
	  if(!tpc1.isValid){
	    tmp_wire.push_back(-9999);
	    tmp_tpc.push_back(-9999);

	  }


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
	primtrk_truth_Z_wire.push_back(tmp_wire);
	primtrk_truth_Z_tpc.push_back(tmp_tpc);


	tmp_primtrk_truth_Z.clear();
	tmp_primtrk_truth_Eng.clear();
	tmp_primtrk_truth_trkide.clear();
	tmp_wire.clear();
	tmp_tpc.clear();
	//new section 

      }// geantGoodParticle
    }//is not real data loop

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

    std::vector<art::Ptr<recob::Hit> > hitlist;
    auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >("hitpdune"); 
    if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);

    // Implementation of required member function here.
    std::vector<art::Ptr<recob::Track> > tracklist;
    auto trackListHandle = evt.getHandle< std::vector<recob::Track> >("pandoraTrack");
    if (trackListHandle) art::fill_ptr_vector(tracklist, trackListHandle);
    else return;
    art::FindManyP<recob::Track> thass(hitListHandle, evt, "pandoraTrack"); //to associate hit just trying

    std::vector<art::Ptr<recob::PFParticle> > pfplist;
    auto PFPListHandle = evt.getHandle< std::vector<recob::PFParticle> >("pandora");
    if (PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);
  
    std::vector<art::Ptr<recob::Cluster>> clusterlist;
    auto clusterListHandle = evt.getHandle< std::vector<recob::Cluster> >("pandora"); // to get information about the hits
    if (clusterListHandle) art::fill_ptr_vector(clusterlist, clusterListHandle);

    art::FindManyP<recob::Cluster> fmcp(PFPListHandle,evt,"pandora");
    art::FindManyP<recob::Track> pftrack(PFPListHandle,evt,"pandoraTrack");
    art::FindManyP<recob::Hit> clhit(clusterListHandle,evt,"pandora");

    std::cout<<"number of pfp_particles "<<pfplist.size()<<std::endl;
    std::cout<<" size of pfParticles testing size "<<pfParticles.size()<<std::endl;
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt,"pandoraTrack"); // to associate tracks and hits


    std::cout<<"outside pfparticle loop"<<std::endl;
    // We can now look at these particles
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

      /////////////////////Michel tagging stage//////////

      double beamstx=-30;
      double beamendx=-30;
      double beamsty=420;
      double beamendy=420;
      double beamstz=30;
      double beamendz=100;    

      if(thisTrack != 0x0){
	if(!beam_cuts.IsBeamlike(*thisTrack, evt, "1")) return;
	beamstx=thisTrack->Start().X();
	beamsty=thisTrack->Start().Y();
	beamstz=thisTrack->Start().Z();
	beamendx=thisTrack->End().X();
	beamendy=thisTrack->End().Y();
	beamendz=thisTrack->End().Z();
	std::cout<<"beamstx "<<beamstx<<std::endl;
      
	//////Michel tagging here
	std::vector<double> secondarystartx1;
	std::vector<double> secondarystarty1;
	std::vector<double> secondarystartz1;
	std::vector<double> secondaryendx1;
	std::vector<double> secondaryendy1;
	std::vector<double> secondaryendz1;
	std::vector<double> dQmichel1;
	std::vector<double> dQtrackbegin1;
	std::vector<double> dQtrackend1;
	std::vector<double> tracklengthsecondary1;
	std::vector<double> primsectheta1; 
	std::vector<int> endhitssecondary1;
	std::vector<int> MtrackID1;
	std::vector<double> trks1;
	std::vector<double> ems1;
	std::vector<double> michels1;
	std::vector<double> nones1;
	size_t NTracks = tracklist.size();
	std::cout<<"number of tracks "<<NTracks<<std::endl;
	for(size_t i=0;i<NTracks;i++){
	  art::Ptr<recob::Track> ptrack(trackListHandle, i);
	  const recob::Track& track = *ptrack;
	  auto pos = track.Vertex();
	  auto end = track.End();
	  int counter1=0;
	  double startx=pos.X();   
	  double starty=pos.Y();
	  double startz=pos.Z();
	  double endx=end.X();
	  double endy=end.Y();
	  double endz=end.Z();
	  if(track.Length()<5) continue;
	  if(TMath::Max(endy,starty)>500 || TMath::Min(endy, starty)<200 || TMath::Max(startx, endx)>0||TMath::Min(startx,endx)<-200||TMath::Max(startz,endz)<30) continue;
	  //if(TMath::Max(endy,starty)>520 || TMath::Min(endy, starty)<150 || TMath::Max(startx, endx)>20||TMath::Min(startx,endx)<-300||TMath::Max(startz,endz)<230) continue;
      
	  std::vector<int> wirenos;
	  std::vector<float> peakts,dqbuff1;
	  std::vector<float> dQstart,dQend;
	  std::vector<double> micheldq;
	  wirenos.clear();peakts.clear();dqbuff1.clear();
	  float peaktime=-1;
	  int wireno=-99999;
	  int tpcno=-1;
	  float zlast0=-99999;
	  float zlast=-99999;
	  std::vector<std::tuple<double,double,double,double,int,double>> buff_ZYXTWQ;
	  buff_ZYXTWQ.clear();
	  double thetavalue=theta12(beamstx,beamendx,beamsty,beamendy,beamstz,beamendz,startx,endx,starty,endy,startz,endz);
	  if(fmthm.isValid()){
	    auto vhit=fmthm.at(i);
	    auto vmeta=fmthm.data(i);
	    for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
	      bool fBadhit = false;
	      if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
		fBadhit = true;
		//cout<<"fBadHit"<<fBadhit<<endl;
		continue;
	      }
	      if (vmeta[ii]->Index()>=tracklist[i]->NumberTrajectoryPoints()){
		throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[i]->NumberTrajectoryPoints()<<" for track index "<<i<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
	      }
	      if (!tracklist[i]->HasValidPoint(vmeta[ii]->Index())){
		fBadhit = true;
		// cout<<"had valid point "<<fBadhit<<endl;
		continue;
	      }
        
	      auto loc = tracklist[i]->LocationAtPoint(vmeta[ii]->Index());
	      if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
	      if (loc.Z()<-100) continue; //hit not on track
	      if(vhit[ii]->WireID().Plane==2){
		buff_ZYXTWQ.push_back(std::make_tuple(loc.Z(),loc.Y(),loc.X(),vhit[ii]->PeakTime(),vhit[ii]->WireID().Wire,vhit[ii]->Integral()));
		wirenos.push_back(vhit[ii]->WireID().Wire);
		peakts.push_back(vhit[ii]->PeakTime());
		zlast=loc.Z();
		if(zlast>zlast0){
		  zlast0=zlast;
		  wireno=vhit[ii]->WireID().Wire;
		  peaktime=vhit[ii]->PeakTime();
		  tpcno=vhit[ii]->WireID().TPC;
		}        
	      }//planenum 2
	    }//loop over vhit
	  }//fmthm valid
	  //save start and end point of each track
	  //taking care of flipped start and end point
	  if(endz<startz){
	    startx=end.X();   
	    starty=end.Y();
	    startz=end.Z();
	    endx=pos.X();
	    endy=pos.Y();
	    endz=pos.Z();
	  }
	  double trk_score=0.0;
	  double em_score=0;
	  double michel_score=0;
	  double none_score=0;
	  for(size_t hitl=0;hitl<hitlist.size();hitl++){
	    std::array<float,4> cnn_out=hitResults.getOutput(hitlist[hitl]);
	    auto & tracks = thass.at(hitlist[hitl].key());
	    if (!tracks.empty() && tracks[0].key()!=ptrack.key() && tracklist[tracks[0].key()]->Length()>25) continue;
	    // if (!tracks.empty() && tracks[0].key()!=ptrack.key() && tracklist[tracks[0].key()]->Length()>25) continue;
	    bool test=true;
	    float peakth1=hitlist[hitl]->PeakTime();
	    int wireh1=hitlist[hitl]->WireID().Wire;
	    for(size_t m=0;m<wirenos.size();m++){
	      if(wireh1==wirenos[m] && peakth1==peakts[m]){
		test=false;
		break;
	      }
	    }
	    if(!test) continue;
	    int planeid=hitlist[hitl]->WireID().Plane;
	    int tpcid=hitlist[hitl]->WireID().TPC;
	    if(abs(wireh1-wireno)<15 && abs(peakth1-peaktime)<100 && planeid==2 && tpcid==tpcno){
	      // if(abs(wireh1-wireno)<20 && abs(peakth1-peaktime)<150 && planeid==2 && tpcid==tpcno){
	      counter1++;
	      // std::cout<<"wireno, counter "<<wireno<<" "<<counter1<<std::endl;
	      micheldq.push_back(hitlist[hitl]->Integral());
	      trk_score+=cnn_out[hitResults.getIndex("track")];
	      em_score+=cnn_out[hitResults.getIndex("em")];
	      michel_score+=cnn_out[hitResults.getIndex("michel")];
	      none_score+=cnn_out[hitResults.getIndex("none")];
	      // std::cout<<"track, em, michel  none"<<cnn_out[hitResults.getIndex("track")]<<" "<<cnn_out[hitResults.getIndex("em")]<<" "<<cnn_out[hitResults.getIndex("michel")]<<" "<<cnn_out[hitResults.getIndex("none")]<<std::endl;
	      //std::cout<<"wire, peaktime selected hit, wires "<<wireno<<" "<<peaktime<<" "<<hitlist[hitl]->PeakTime()<<" "<<hitlist[hitl]->WireID().Wire<<std::endl;
	    }
	  }//hitlist loop
	  if(buff_ZYXTWQ.size()<10) continue;
	  sort(buff_ZYXTWQ.begin(),buff_ZYXTWQ.end());
	  dQstart.clear(); dQend.clear();
	  int qi11=buff_ZYXTWQ.size();
	  for(int qi=5;qi<TMath::Min(15,qi11);qi++){
	    dQstart.push_back(std::get<5>(buff_ZYXTWQ[qi]));
	  }
    
	  for(int qi=qi11-5;qi<qi11;qi++){
	    dQend.push_back(std::get<5>(buff_ZYXTWQ[qi]));
	  }
	  secondarystartx1.push_back(startx);
	  secondarystarty1.push_back(starty);
	  secondarystartz1.push_back(startz);
	  secondaryendx1.push_back(endx);
	  secondaryendy1.push_back(endy);
	  secondaryendz1.push_back(endz);
	  endhitssecondary1.push_back(counter1);
	  tracklengthsecondary1.push_back(track.Length());
	  dQmichel1.push_back(TMath::Median(micheldq.size(),&micheldq[0]));
	  dQtrackbegin1.push_back(TMath::Median(dQstart.size(),&dQstart[0]));
	  dQtrackend1.push_back(TMath::Median(dQend.size(),&dQend[0]));
	  primsectheta1.push_back(thetavalue);
	  MtrackID1.push_back(track.ID());
	  trks1.push_back(trk_score);
	  ems1.push_back(em_score);
	  michels1.push_back(michel_score);
	  nones1.push_back(none_score);
	  std::cout<<"avg, trackscore, emscore, michelscore, nonescore "<<trk_score<<" "<<em_score<<" "<<michel_score<<" "<<none_score<<std::endl; 
	}//Ntracks
	secondarystartx.push_back(secondarystartx1);
	secondarystarty.push_back(secondarystarty1);
	secondarystartz.push_back(secondarystartz1);
	secondaryendx.push_back(secondaryendx1);
	secondaryendy.push_back(secondaryendy1);
	secondaryendz.push_back(secondaryendz1);
	endhitssecondary.push_back(endhitssecondary1);
	tracklengthsecondary.push_back(tracklengthsecondary1);
	dQmichel.push_back(dQmichel1);
	dQtrackbegin.push_back(dQtrackbegin1);
	dQtrackend.push_back(dQtrackend1);
	primsectheta.push_back(primsectheta1);
	MtrackID.push_back(MtrackID1);

	trackscore.push_back(trks1);
	emscore.push_back(ems1);
	michelscore.push_back(michels1);
	nonescore.push_back(nones1);
      
	endhitssecondary1.clear();
	secondarystartx1.clear();
	secondaryendx1.clear();
	secondarystarty1.clear();
	secondaryendy1.clear();
	secondarystartz1.clear();
	secondaryendz1.clear();
	dQmichel1.clear();
	primsectheta1.clear();
	dQtrackbegin1.clear();
	dQtrackend1.clear();
	tracklengthsecondary1.clear();
	MtrackID1.clear();
	trks1.clear();
	ems1.clear();
	michels1.clear();
	nones1.clear();
      }
      ///////////////End of Michel checking
  
      /////Michel tagging ends here
      if(thisTrack != 0x0){
	if(!beam_cuts.IsBeamlike(*thisTrack, evt, "1")) return;
	// Get the true mc particle
	const simb::MCParticle* mcparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, evt, fTrackerTag);
	if(mcparticle!=0x0){
	  std::cout<<"ftruth pdg "<<mcparticle->PdgCode()<<std::endl;
	  ftruthpdg=mcparticle->PdgCode();
	  truthid=mcparticle->TrackId();
	  fprimary_truth_Pdg= mcparticle->PdgCode();
	  std::cout<<"primary_truth_Pdg "<<fprimary_truth_Pdg<<std::endl;
	  fprimary_truth_Isbeammatched=0;

	  fprimary_truth_TrackId          = mcparticle->TrackId();
	  fprimary_truth_StartPosition[3] = mcparticle->T();
	  fprimary_truth_EndPosition[3]   = mcparticle->EndT();
	  fprimary_truth_P                = mcparticle->P();
	  fprimary_truth_Momentum[3]      = mcparticle->E();
	  fprimary_truth_Pt               = mcparticle->Pt();
	  fprimary_truth_Mass             = mcparticle->Mass();
	  fprimary_truth_EndMomentum[0]   = mcparticle->EndPx();
	  fprimary_truth_EndMomentum[1]   = mcparticle->EndPy();
	  fprimary_truth_EndMomentum[2]   = mcparticle->EndPz();
	  fprimary_truth_EndMomentum[3]   = mcparticle->EndE();
	  fprimary_truth_Theta            = mcparticle->Momentum().Theta();
	  fprimary_truth_Phi              = mcparticle->Momentum().Phi();
	  fprimary_truth_NDaughters       = mcparticle->NumberDaughters();
	  if(beamid==truthid) fprimary_truth_Isbeammatched=1;
	}//mcparticle loop


	fisprimarytrack               = 1;
	fisprimaryshower              = 0;
	fprimaryID                    = thisTrack->ID();
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
	std::vector<int> tmp_zwire;
	std::vector<int> tmp_ztpc;


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

	      /* if(primtrk_pos.Z()>=0 && primtrk_pos.Z()<=230){
		 geo::WireID wireID = fWireReadoutGeom.NearestWireID(pos_true, 2);
		 if (!wireID) wireID = fWireReadoutGeom.Plane(2).ClosestWireID(wireID);
		 tmp_zwire.push_back(wireID.Wire);
		 tmp_ztpc.push_back(wireID.TPC);
		 }
		 if(primtrk_pos.Z()<0 || primtrk_pos.Z()>230){
		 tmp_zwire.push_back(-99999);
		 tmp_ztpc.push_back(-99999);
		 }
	      */
              geo::TPCID tpc2 = fGeometry->FindTPCAtPosition(primtrk_pos);
	      if(tpc2.isValid){
		int tpc_no=tpc2.TPC;
		geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
		geo::WireID wireID;
		//geo::WireID wireID = fWireReadoutGeom.NearestWireID(pos_true, planeID);
		try{
                  wireID = fWireReadoutGeom.Plane(planeID).NearestWireID(primtrk_pos);
		}
		catch(geo::InvalidWireError const& e) {
		  wireID = e.suggestedWireID(); // pick the closest valid wire
		}
		tmp_zwire.push_back(wireID.Wire);
		tmp_ztpc.push_back(wireID.TPC);
	      }
	      if(!tpc2.isValid){
		tmp_zwire.push_back(-9999);
		tmp_ztpc.push_back(-9999);
	      }
	    } //loop over hits
	  } //only collection plane
	}//calovector

	if(tmp_primtrk_dqdx.size()!=0){
	  if(tmp_primtrk_hitz[0]>tmp_primtrk_hitz[tmp_primtrk_hitz.size()-1]){
	    std::reverse(tmp_primtrk_hitz.begin(),tmp_primtrk_hitz.end());
	    std::reverse(tmp_primtrk_hity.begin(),tmp_primtrk_hity.end());
	    std::reverse(tmp_primtrk_hitx.begin(),tmp_primtrk_hitx.end());
	    std::reverse(tmp_primtrk_pitch.begin(),tmp_primtrk_pitch.end());
	    std::reverse(tmp_primtrk_dedx.begin(),tmp_primtrk_dedx.end());
	    std::reverse(tmp_primtrk_dqdx.begin(),tmp_primtrk_dqdx.end());
	    std::reverse(tmp_primtrk_resrange.begin(),tmp_primtrk_resrange.end());

	    std::reverse(tmp_ztpc.begin(),tmp_ztpc.end());
	    std::reverse(tmp_zwire.begin(),tmp_zwire.end());
	  }
	  primtrk_dqdx.push_back(tmp_primtrk_dqdx);
	  primtrk_resrange.push_back(tmp_primtrk_resrange);
	  primtrk_dedx.push_back(tmp_primtrk_dedx);
	  primtrk_hitx.push_back(tmp_primtrk_hitx);
	  primtrk_hity.push_back(tmp_primtrk_hity);
	  primtrk_hitz.push_back(tmp_primtrk_hitz);
	  primtrk_pitch.push_back(tmp_primtrk_pitch);

	  primtrk_hitz_wire.push_back(tmp_zwire);
	  primtrk_hitz_tpc.push_back(tmp_ztpc);

	}
	tmp_primtrk_dqdx.clear();
	tmp_primtrk_resrange.clear();
	tmp_primtrk_dedx.clear();
	tmp_primtrk_hitx.clear();
	tmp_primtrk_hity.clear();
	tmp_primtrk_hitz.clear();
	tmp_primtrk_pitch.clear();
	tmp_ztpc.clear();
	tmp_zwire.clear();
     
	for(size_t k = 0; k < calovector.size() && k<3; k++){
	  fprimaryKineticEnergy[k] = calovector[k].KineticEnergy();
	  fprimaryRange[k] = calovector[k].Range();
	  //const std::vector< double > & dedxvec = calovector[k].dEdx();
	}

     

	int planenum=999;
	float zpos=-999;

	//////////Section for looking into integral dQ vs hitz position//////////// 
	//hits and calorimetry loop
	if(fmthm.isValid()){
	  auto vhit=fmthm.at(fprimaryID);
	  auto vmeta=fmthm.data(fprimaryID);
	  for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
	    bool fBadhit = false;
	    if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
	      fBadhit = true;
	      //cout<<"fBadHit"<<fBadhit<<endl;
	      continue;
	    }
	    if (vmeta[ii]->Index()>=tracklist[fprimaryID]->NumberTrajectoryPoints()){
	      throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[fprimaryID]->NumberTrajectoryPoints()<<" for track index "<<fprimaryID<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
	    }
	    if (!tracklist[fprimaryID]->HasValidPoint(vmeta[ii]->Index())){
	      fBadhit = true;
	      // cout<<"had valid point "<<fBadhit<<endl;
	      continue;
	    }
	
	    auto loc = tracklist[fprimaryID]->LocationAtPoint(vmeta[ii]->Index());
	    // xpos=loc.X();
	    // ypos=loc.Y();
	    zpos=loc.Z();
	    //	cout<<"x, y, z "<<xpos<<"  "<<ypos<<"  "<<zpos<<endl;
	    //	cout<<"BadHit"<<fBadhit<<endl;
	    if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
	    if (zpos<-100) continue; //hit not on track
	    planenum=vhit[ii]->WireID().Plane;
	    if(planenum==2){
	      // std::cout<<"inside loop"<<std::endl;
	      // std::array<float,3> cnn_out=hitResults.getOutput(vhit[ii]);
	      peakT_2.push_back(vhit[ii]->PeakTime());
	      int_2.push_back(vhit[ii]->Integral());
	      hitz_2.push_back(zpos);
	      // inelscore.push_back(cnn_out[hitResults.getIndex("inel")]);
	      // elscore.push_back(cnn_out[hitResults.getIndex("el")]);
	      // nonescore.push_back(cnn_out[hitResults.getIndex("none")]);

	      //  std::cout<<"inel, el, none  score"<<cnn_out[hitResults.getIndex("inel")]<<" "<<cnn_out[hitResults.getIndex("el")]<<" "<<cnn_out[hitResults.getIndex("none")]<<std::endl;
	      // std::cout<<"peaktime "<<vhit[ii]->PeakTime()<<std::endl;	
	    }//planenum 2
	    if(planenum==1){
	      int_1.push_back(vhit[ii]->Integral());
	      hitz_1.push_back(zpos);	
	    }//planenum 1
	    if(planenum==0){
	      int_0.push_back(vhit[ii]->Integral());
	      hitz_0.push_back(zpos);	
	    }//planenum 0
	  }//loop over vhit
	}//fmthm valid
	//hits and calorimetry loop   


	//**************Section for dE/dx correction*****************************************************************************************************************************//
	//**********************************************************************************************************************************************************************//
	std::vector<float> Stw, Endw, Stt, Endt, Stwires, Endwires, Stticks, Endticks, TPCb, TPCcl, primwireb, primtickb, wirebuf, tickbuf, primdqb;
	std::vector<std::vector<float> > primwire;
	std::vector<std::vector<float> > clwire;
	std::vector<std::vector<float> > primtick;
	std::vector<std::vector<float> >  cltick;
	Stw.clear(); Endw.clear(); Stt.clear();  Endt.clear();  Stwires.clear();  Endwires.clear(); Stticks.clear(); Endticks.clear(); TPCb.clear(); TPCcl.clear();
	float den;
	float numw, numt,wire_no,ticks_no;
	for(size_t p1=0;p1<pfplist.size();p1++){
	  std::vector<art::Ptr<recob::Track>> trk=pftrack.at(p1);
	  std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(p1);
	  for(size_t c1=0;c1<allClusters.size();c1++){
	    if(trk.size() && int(trk[0].key())==fprimaryID) std::cout<<"cluster tpc and plane "<<allClusters[c1]->Plane().TPC<<"  "<<allClusters[c1]->Plane().Plane<<" size "<<clhit.at(allClusters[c1].key()).size()<<std::endl;
	    if(allClusters[c1]->Plane().TPC!=1) continue;
	    if(allClusters[c1]->Plane().Plane==0 && trk.size() && int(trk[0].key())==fprimaryID){ 
	      std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	      for(size_t h1=0;h1<allHits.size();h1++){
		primwireb.push_back(allHits[h1]->WireID().Wire);
		primtickb.push_back(allHits[h1]->PeakTime());
		primdqb.push_back(allHits[h1]->Integral());
	      }
	      wireno_0.push_back(primwireb);
	      peakTime_0.push_back(primtickb);
	      dq_0.push_back(primdqb);
	      primwireb.clear();
	      primtickb.clear();
	      primdqb.clear();
	    }
	    if(allClusters[c1]->Plane().Plane==1 && trk.size() && int(trk[0].key())==fprimaryID){ 
	      std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	      for(size_t h1=0;h1<allHits.size();h1++){
		primwireb.push_back(allHits[h1]->WireID().Wire);
		primtickb.push_back(allHits[h1]->PeakTime());
		primdqb.push_back(allHits[h1]->Integral());
	      }
	      wireno_1.push_back(primwireb);
	      peakTime_1.push_back(primtickb);
	      dq_1.push_back(primdqb);
	      primwireb.clear();
	      primtickb.clear();
	      primdqb.clear();
	    }
	  
	    if(allClusters[c1]->Plane().Plane!=2) continue;
	    if(trk.size() && int(trk[0].key())==fprimaryID){
	      std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	      Stw.push_back(allClusters[c1]->StartWire());
	      Endw.push_back(allClusters[c1]->EndWire());
	      Stt.push_back(allClusters[c1]->StartTick());
	      Endt.push_back(allClusters[c1]->EndTick());
	      TPCb.push_back(allClusters[c1]->Plane().TPC);
	      for(size_t h1=0;h1<allHits.size();h1++){
		primwireb.push_back(allHits[h1]->WireID().Wire);
		primtickb.push_back(allHits[h1]->PeakTime());
		primdqb.push_back(allHits[h1]->Integral());
	      }
	      wireno_2.push_back(primwireb);
	      peakTime_2.push_back(primtickb);
	      dq_2.push_back(primdqb);
	      primwire.push_back(primwireb);
	      primtick.push_back(primtickb);
	      primwireb.clear();
	      primtickb.clear();
	      primdqb.clear();
	    }
	    else if(trk.size()){
	      std::vector<art::Ptr<recob::Hit>> allHits=clhit.at(allClusters[c1].key());
	      // std::cout<<"size of hits "<<allHits.size()<<std::endl;
	      Stwires.push_back(allClusters[c1]->StartWire());
	      Endwires.push_back(allClusters[c1]->EndWire());
	      Stticks.push_back(allClusters[c1]->StartTick());
	      Endticks.push_back(allClusters[c1]->EndTick());
	      TPCcl.push_back(allClusters[c1]->Plane().TPC);
	      for(size_t h2=0;h2<allHits.size();h2++){
		wirebuf.push_back(allHits[h2]->WireID().Wire);
		tickbuf.push_back(allHits[h2]->PeakTime());
	      }
	      clwire.push_back(wirebuf);
	      cltick.push_back(tickbuf);
	      wirebuf.clear();
	      tickbuf.clear();
	    }
	  }
	}
	//solving the equations here
	for(size_t clt=0;clt<Stw.size();clt++){
	  for(size_t cl1=0;cl1<Stwires.size();cl1++){
	    if(TPCcl[cl1]!=TPCb[clt]) continue;
	    den=(Stw[clt]-Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]-Endwires[cl1]);
	    bool testv=false;
	    if(Stwires[cl1]==Endwires[cl1]){
	      wire_no=Stwires[cl1];
	      ticks_no=Stticks[cl1];
	      testv=true;
	    }
	    if(den==0  && !testv) continue;
	    if(!testv){
	      numw=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stwires[cl1]-Endwires[cl1])-(Stw[clt]-Endw[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
	      numt=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
	      wire_no=numw/den;
	      ticks_no=numt/den;
	    }
	    if(((Stw[clt]-20<wire_no && Endw[clt]+20>wire_no)||(Stw[clt]+20>wire_no && Endw[clt]<wire_no-20))&&((Stt[clt]-100<ticks_no && Endt[clt]+100>ticks_no)||(Stt[clt]+100>ticks_no && Endt[clt]-100<ticks_no)) && ((Stwires[cl1]-20<wire_no && Endwires[cl1]+20>wire_no)||(Stwires[cl1]+20>wire_no && Endwires[cl1]-20<wire_no)) && ((Stticks[cl1]-100<ticks_no && Endticks[cl1]+100>ticks_no)||(Stticks[cl1]+100>ticks_no && Endticks[cl1]-100<ticks_no))){
	      double xyzStart[3];
	      double xyzEnd[3];
	      unsigned int wireno=std::round(wire_no);
	      // geo::WireID wireid(0,TPCb[clt],2,wireno);
	      if(wireno>=0 && wireno<=479){
                fWireReadoutGeom.WireEndPoints(geo::WireID(0,TPCb[clt],2,wireno), xyzStart, xyzEnd);
		Zintersection.push_back(xyzStart[2]);
		timeintersection.push_back(ticks_no);
		// std::cout<<"intersecting "<<xyzStart[2]<<std::endl;
	      }
	      else{
		Zintersection.push_back(wireno*0.47);
		timeintersection.push_back(ticks_no);
	      }
	      //go for the second solution now
	      float pws=wireno; //p=primary w=wire s=small b=big t=tick c=other clusters
	      float pwb=wireno;
	      float pts=ticks_no;
	      float ptb=ticks_no;
	      float smallest=wireno;
	      float biggest=wireno;
	      float small, big;
	      //finding new variables for beam track
	      for(size_t hit1=0;hit1<primwire[clt].size();hit1++){
		if(primwire[clt][hit1]>=wireno-30 && primwire[clt][hit1]<=wireno+30){
		  if(primwire[clt][hit1]>=wireno-30 && primwire[clt][hit1]<=wireno){
		    small=primwire[clt][hit1];
		    if(small<smallest){
		      pws=primwire[clt][hit1];
		      pts=primtick[clt][hit1];
		    }
		    smallest=primwire[clt][hit1];
		  }
		  if(primwire[clt][hit1]<=wireno+30 && primwire[clt][hit1]>=wireno){
		    big=primwire[clt][hit1];
		    if(big>biggest){
		      pwb=primwire[clt][hit1];
		      ptb=primtick[clt][hit1];
		    }
		    biggest=primwire[clt][hit1];
		  }
		  
		}//if +_30 loop
		
	      }//hit1 loop
	      //finding new variables for remaining clusters
	      float cws=wireno;
	      float cwb=wireno;
	      float cts=ticks_no;
	      float ctb=ticks_no;
	      smallest=wireno;
	      biggest=wireno;
	      float csmall, cbig;
	      for(size_t hit1=0;hit1<clwire[cl1].size();hit1++){
		if(clwire[cl1][hit1]>=wireno-30 && clwire[cl1][hit1]<=wireno+30){
		  if(clwire[cl1][hit1]>=wireno-30 && clwire[cl1][hit1]<=wireno){
		    csmall=clwire[cl1][hit1];
		    if(csmall<smallest){
		      cws=clwire[cl1][hit1];
		      cts=cltick[cl1][hit1];
		    }
		    smallest=clwire[cl1][hit1];
		  }
		  if(clwire[cl1][hit1]<=wireno+30 && clwire[cl1][hit1]>=wireno){
		    cbig=clwire[cl1][hit1];
		    if(cbig>biggest){
		      cwb=clwire[cl1][hit1];
		      ctb=cltick[cl1][hit1];
		    }
		    biggest=clwire[cl1][hit1];
		  }
			
		}//if +_30 loop
	      }//hit1 loop
	      float ans[2];
	      soln(pws, pwb, cws, cwb, pts, ptb, cts, ctb,ans);
	      // std::cout<<"answer values "<<ans[0]<<"  "<<ans[1]<<std::endl;
	      // std::cout<<"wire number and new start end wire "<<wireno<<"  st w "<<pws<<"  end w "<<pwb<<" time and new times "<<ticks_no<<"  "<<pts<<"  "<<ptb<<" same for intersecting cluster "<<cws<<" "<<cwb<<" "<<cts<<"  "<<ctb<<" hits in a cluster "<<clwire[cl1].size()<<"start wire and ticks "<<Stwires[cl1]<<" "<<Endwires[cl1]<<" "<<Stticks[cl1]<<"  "<<Endticks[cl1]<<std::endl;    
	      int wir=round(ans[0]);
	      double xyzStart1[3];
	      double xyzEnd1[3];

	      /* fWireReadoutGeom.WireEndPoints(0,2,2,10, xyzStart1, xyzEnd1);//will remove this line
		 stdd::cout<<"Wire start X in the postive direction "<<xyzStart1[0]<<std::endl;
		 fWireReadoutGeom.WireEndPoints(0,1,2,10, xyzStart1, xyzEnd1);//will remove this line
		 std::cout<<"Wire start X in the negative direction "<<xyzStart1[0]<<std::endl;*/

	      if(((Stw[clt]-5<ans[0] && Endw[clt]+5>ans[0])||(Stw[clt]+5>ans[0] && Endw[clt]<ans[0]-5))&&((Stt[clt]-50<ans[1] && Endt[clt]+50>ans[1])||(Stt[clt]+50>ans[1] && Endt[clt]-50<ans[1])) && ((Stwires[cl1]-5<ans[0] && Endwires[cl1]+5>ans[0])||(Stwires[cl1]+5>ans[0] && Endwires[cl1]-5<ans[0])) && ((Stticks[cl1]-50<ans[1] && Endticks[cl1]+50>ans[1])||(Stticks[cl1]+50>ans[1] && Endticks[cl1]-50<ans[1]))){
		if(wir>=0 && wir<=479){
                  fWireReadoutGeom.WireEndPoints(geo::WireID(0,TPCb[clt],2,wir), xyzStart1, xyzEnd1);
		  Zintersection1.push_back(xyzStart1[2]);
		  timeintersection1.push_back(ans[1]);
		  std::cout<<"wire X position is here  "<<xyzStart1[0]<<std::endl;
		  std::cout<<"wire X position End  is here  "<<xyzEnd1[0]<<std::endl;
		}
		else{
		  Zintersection1.push_back(wire_no*0.47);
		  timeintersection1.push_back(ticks_no);
		  ans[0]=wire_no;
		  ans[1]=ticks_no;
		}
	      }
	      else{
		Zintersection1.push_back(99999);
		timeintersection1.push_back(99999);
	      }
					     
	      //  if((pws<=ans[0] && pwb>=ans[0])&&(cws<=ans[0] && cwb>=ans[0])&&(pts<=ans[1] && ptb>=ans[1])&&(cts<=ans[1] && ctb>=ans[1])){
	      if(((Stw[clt]<=ans[0] && Endw[clt]>=ans[0])||(Stw[clt]>=ans[0] && Endw[clt]<=ans[0]))&&((Stt[clt]<=ans[1] && Endt[clt]>ans[1])||(Stt[clt]>ans[1] && Endt[clt]<ans[1])) && ((Stwires[cl1]<ans[0] && Endwires[cl1]>ans[0])||(Stwires[cl1]>ans[0] && Endwires[cl1]<ans[0])) && ((Stticks[cl1]<ans[1] && Endticks[cl1]>ans[1])||(Stticks[cl1]>ans[1] && Endticks[cl1]<ans[1]))){
		deltaZint.push_back(0);
		deltatimeint.push_back(0);
	      }
	      else if((Stw[clt]<=ans[0] && Endw[clt]>=ans[0])||(Stw[clt]>=ans[0] && Endw[clt]<=ans[0])){
		deltaZint.push_back(abs(std::min(abs(ans[0]-Stwires[cl1]),abs(ans[0]-Endwires[cl1]))));
		deltatimeint.push_back(abs(std::min(abs(ans[1]-Stticks[cl1]),abs(ans[1]-Endticks[cl1]))));
	      }
	      else if((Stwires[cl1]<=ans[0] && Endwires[cl1]>=ans[0])||(Stwires[cl1]>=ans[0] && Endwires[cl1]<=ans[0])){
		deltaZint.push_back(abs(std::min(abs(ans[0]-Stw[clt]),abs(ans[0]-Endw[clt]))));
		deltatimeint.push_back(abs(std::min(abs(ans[1]-Stt[clt]),abs(ans[1]-Endt[clt]))));
	      }
	      else{
		deltaZint.push_back(99999);
		deltatimeint.push_back(99999);
	      }

	    }//if solution lies within wire coordinates
	    //	}//if wire coordinates
	  }//cl1 loop
	}//clt
	primwire.clear();
	clwire.clear();
	primtick.clear();
	cltick.clear();
      
	////*****************section for dE/dx correction*******************************////

	
      }//this track
     
  


      //*********************************************//
      //*********************************************//
      ///////////Showere case and daughter trajectories
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

  void protoana::mcsXsection::endJob(){

  }

  void protoana::mcsXsection::FillCosmicsTree(art::Event const & evt, std::string pfParticleTag){

    // To fill

  }



  void protoana::mcsXsection::Initialise(){
  
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
    beamtrk_z_wire.clear();
    beamtrk_z_tpc.clear();


    beamMomentum_spec.clear();
    beamPosx_spec.clear();
    beamPosy_spec.clear();
    beamPosz_spec.clear();
    beamDirx_spec.clear();
    beamDiry_spec.clear();
    beamDirz_spec.clear();


 
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
    int_2.clear();
    hitz_2.clear();
    // inelscore.clear();
    // elscore.clear();
    // nonescore.clear();
    int_1.clear();
    hitz_1.clear();
    int_0.clear();
    hitz_0.clear();
    peakT_2.clear();
    primtrk_dqdx.clear();
    primtrk_resrange.clear();
    primtrk_dedx.clear();
    primtrk_range.clear();
    primtrk_hitx.clear();
    primtrk_hity.clear();
    primtrk_hitz.clear();
    primtrk_hitz_wire.clear();
    primtrk_hitz_tpc.clear();

    primtrk_pitch.clear();
    primtrk_truth_Z.clear();
    primtrk_truth_Z_wire.clear();
    primtrk_truth_Z_tpc.clear();


    primtrk_truth_Eng.clear();
    primtrk_truth_trkide.clear();
    interactionX.clear();
    interactionY.clear();
    interactionZ.clear();
    interactionU.clear();
    interactionV.clear();
    interactionW.clear();
    interactionT.clear();


    interactionProcesslist.clear();
    interactionAngles.clear();
    interactionAnglesUT.clear();
    interactionAnglesVT.clear();
    interactionAnglesZT.clear();
    Zintersection.clear();
    timeintersection.clear();
    Zintersection1.clear();
    timeintersection1.clear();
    deltaZint.clear();
    deltatimeint.clear();
    wireno_0.clear();
    wireno_1.clear();
    wireno_2.clear();
    dq_0.clear();
    dq_1.clear();
    dq_2.clear();
    peakTime_0.clear();
    peakTime_1.clear();
    peakTime_2.clear();
    //Michel tagging
    endhitssecondary.clear();
    secondarystartx.clear();
    secondaryendx.clear();
    secondarystarty.clear();
    secondaryendy.clear();
    secondarystartz.clear();
    secondaryendz.clear();
    dQmichel.clear();
    primsectheta.clear();
    dQtrackbegin.clear();
    dQtrackend.clear();
    tracklengthsecondary.clear();
    MtrackID.clear();
    trackscore.clear();
    emscore.clear();
    michelscore.clear();
    nonescore.clear();
    //geant4 reweight
    n_g4rw=0;
    g4rw_react_st=1;
    g4rw_d_react=0;
    g4rw_elast_st=1;
    g4rw_d_elast=0;
				
    g4rw_set_weights.clear();
    g4rw_elast.clear();
    g4rw_react.clear();

  }

  DEFINE_ART_MODULE(protoana::mcsXsection)
