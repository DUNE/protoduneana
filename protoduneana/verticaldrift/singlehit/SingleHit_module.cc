////////////////////////////////////////////////////////////////////////
// Class:       SingleHit
// Plugin Type: analyzer (Unknown Unknown)
// File:        SingleHit_module.cc
//
// Generated at Thu Feb  1 04:19:28 2024 by Emile Lavaut using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include <limits>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

//LArSoft
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/ServicePack.h" 
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RawData/RDTimeStamp.h"

// ROOT
#include "Math/ProbFunc.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TGraph.h"
#include "TEllipse.h"

using std::vector;
using std::string;

typedef struct 
{ 
  float y, z; 
  int group, index;  
} point_t, *point;

typedef struct                                                                
{                                                                                 
  float Sumy , Sumz , ECol , EInd1 , EInd2 , PeakTime;                  
  int Npoint, NCol , NInd1 , NInd2;
  std::list<int> lChannelCol , lChannelInd1 , lChannelInd2;
  std::vector<int> vMCPDG , vMCMOMpdg , vNOF;
  std::vector<float> vMCWEI;
  std::vector<float> vMCX , vMCY , vMCZ;
  std::vector<std::string> vMCGenTag;
} TempCluster ;                            
                              
typedef struct
{ 
  float y , z , ECol , EInd1 , EInd2 , PeakTime;
  int Npoint, NCol , NInd1 , NInd2 , NOF;
  std::vector<int> vMCPDG , vMCMOMpdg;
  std::vector<float> vMCWEI;
  std::vector<float> vMCX , vMCY , vMCZ;
  std::vector<std::string> vMCGenTag;
} Cluster ; 

//stucture used for assocition between a grid coordinate and a int 
struct GridHasher {
    std::size_t operator()(const std::tuple<int, int, int>& key) const {
        // Use a combination of prime numbers to create a unique hash from 3D coordinates
        return std::get<0>(key) * 73856093 ^ std::get<1>(key) * 19349663 ^ std::get<2>(key) * 83492791;
    }
};

namespace pdvdana {

  using std::vector;
  using std::string;
  class SingleHit;
}


class pdvdana::SingleHit : public art::EDAnalyzer {
public:
  explicit SingleHit(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SingleHit(SingleHit const&) = delete;
  SingleHit(SingleHit&&) = delete;
  SingleHit& operator=(SingleHit const&) = delete;
  SingleHit& operator=(SingleHit&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  float degTOrad = 3.14159/180.; // rad / deg
  float fEltofC  = 1./1.60E-4;      // e- / fC
  float fADCtoEl = 5E-3;         // ADC x tick / e-

private:

  // Declare member data here.
  TTree *tHitTree;
  TTree *tClusterTree;

  // Hit Tree variables
  unsigned int fEventID = 0;

  geo::WireID fWire;

  int fHitNumber     = -999;
  int fPlane         = -999;
  int fChannel       = -999;
  int fHitWidth      = -999;
  int fTimeIsolation = -999;
  int fCoincidence   = -999;
  int fNearOrFarToTheBeam = -999; // 1 means hits is near the beam (ie right/top for HD/VD) -1 means far (ie left/bottom for HD/VD)

  float fEnergy         = -999.;
  float fPeakTime       = -999.;
  float fSigmaPeakTime  = -999.;
  float fRMS            = -999.;
  float fAmplitude      = -999.;
  float fSigmaAmplitude = -999.;
  float fGoodnessOfFit  = -999.;
  float fIntegral       = -999.;
  float fSigmaIntegral  = -999.;

  double CRP_T0 = -999;

  std::list<geo::WireID>   lWireInd1; //working variable
  std::list<int>           lChannelInd1;
  std::list<float>         lEnergyInd1;
  std::list<float>         lPeakTimeInd1;
  std::list<float>         lPeakAmpInd1;
  std::list<float>         lYInd1;
  std::list<float>         lZInd1;
  std::list<int>           lChIntersectInd1;
  std::list<float>         lEIntersectInd1;

  std::list<geo::WireID>   lWireInd2; //working variable
  std::list<int>           lChannelInd2;
  std::list<float>         lEnergyInd2;
  std::list<float>         lPeakTimeInd2;
  std::list<float>         lPeakAmpInd2;
  std::list<float>         lYInd2;
  std::list<float>         lZInd2;
  std::list<int>           lChIntersectInd2;
  std::list<float>         lEIntersectInd2;

  std::list<float> lYPoint;
  std::list<float> lZPoint;
  std::list<float> lEInd1Point;
  std::list<float> lEInd2Point;
  std::list<int>   lChInd1Point;
  std::list<int>   lChInd2Point;

  std::vector<int>         vMCPart_pdgCode;
  std::vector<float>       vMCPart_weight;
  std::vector<int>         vMCPart_mother;
  std::vector<int>         vMCPart_motherPdg;
  std::vector<std::string> vGenerator_tag;
  std::vector<float>       vMCPart_Endx;
  std::vector<float>       vMCPart_Endy;
  std::vector<float>       vMCPart_Endz;
  std::vector<float>       vMCPart_Startx;
  std::vector<float>       vMCPart_Starty;
  std::vector<float>       vMCPart_Startz;

  // Cluster Tree Variables
  int NCluster;
  bool fAsConverged;

  std::vector<float> vY_point;
  std::vector<float> vZ_point;

  std::vector<float> vZCluster;
  std::vector<float> vYCluster;
  std::vector<float> vEColCluster;
  std::vector<float> vEInd1Cluster;
  std::vector<float> vEInd2Cluster;
  std::vector<float> vPTCluster;

  std::vector<int>   vNoFCluster;
  std::vector<int>   vNPointCluster;
  std::vector<int>   vNColCluster; 
  std::vector<int>   vNInd1Cluster; 
  std::vector<int>   vNInd2Cluster;

  std::vector<std::string> vMCGenTagCluster;
  std::vector<int>         vMCMOMpdgCluster;
  std::vector<int>         vMCPDGCluster;
  std::vector<float>       vMCWeightCluster;
  std::vector<float>       vMCXCluster;
  std::vector<float>       vMCYCluster;
  std::vector<float>       vMCZCluster;

  //Input variables
  std::string fSpacePointLabel;
  std::string fClusterLabel;
  std::string fTrackLabel;
  std::string fHitLabel;
  std::string fG4Label;
  std::string fRDTLabel;

  int   LogLevel;
  int   fMultiplicity;
  float fTickTimeInMus;

  float fRadiusInt;
  float fRadiusExt;
  float fElectronVelocity;
  float fCoincidenceWd1_left;
  float fCoincidenceWd1_right;
  float fCoincidenceWd2_left;
  float fCoincidenceWd2_right;

  float fPitch;
  float fPitchMultiplier;

  bool  bIs3ViewsCoincidence;
  bool  bHitTree;
  bool  bIsPDVD;
  bool  bIsPDHD;
  bool  bVetoTrack;

  float fNumberInitClusters;
  float fMaxSizeCluster;
  float fMinSizeCluster;
  float fClusterSizeMulti;
  float fNumberConvStep;
  float fCovering;

  //float fCalibration;

  // geometry
  const geo::WireReadoutGeom& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  const geo::Geometry* fGeom;
  float fgeoXmin = 1e6; 
  float fgeoXmax =-1e6; 
  float fgeoYmin = 1e6; 
  float fgeoYmax =-1e6; 
  float fgeoZmin = 1e6; 
  float fgeoZmax =-1e6; 

  // working variables
 
  int fHitCounter = 0;
  int fNHits      = 0;
  int fAmbiguousHit = 0;

  std::list<int> lSingleIndex;
  std::list<int> lIsolatedIndex;

  std::vector<float> vYPointByEvent;
  std::vector<float> vZPointByEvent;

  std::vector<float> vEnergyColByEvent;
  std::vector<float> vEInd1PointByEvent;
  std::vector<float> vEInd2PointByEvent; 

  std::vector<float> vPeakTimeColByEvent;

  std::vector<int>   vChannelColByEvent;
  std::vector<int>   vChInd1PointByEvent;
  std::vector<int>   vChInd2PointByEvent;
  std::vector<int>   vNoFByEvent;

  std::vector<int>         vMCMOMpdgByEvent;
  std::vector<int>         vMCPDGByEvent;
  std::vector<float>       vMCWeightByEvent;
  std::vector<float>       vMCXByEvent;
  std::vector<float>       vMCYByEvent;
  std::vector<float>       vMCZByEvent;

  std::vector<std::string> vGeneratorTagByEvent;

  // veto track vector
  std::vector<float> vVetoTrackStartX;
  std::vector<float> vVetoTrackStartY;
  std::vector<float> vVetoTrackStartZ;
  std::vector<float> vVetoTrackEndX;
  std::vector<float> vVetoTrackEndY;
  std::vector<float> vVetoTrackEndZ;

  //function needed
  void print(std::vector<float> v);
  void print( std::vector<std::vector<int>> vv );
  bool AllSame(std::vector<int> v);

  std::vector<std::string> GetGeneratorTag(  art::Event const& e , std::string fG4Label , int LogLevel , art::ServiceHandle<cheat::BackTrackerService> bt_serv );

  bool Inside( int k , std::list<int> liste);

  float GetDist( float x0 , float y0 , float z0 , float x1 , float y1 , float z1 );

  void GetSingle(art::Event const & ev, std::string HitLabel, std::list<int> & index_list_single, int const Multiplicity);

  int NearOrFar( bool IsPDVD , bool IsPDHD , const recob::Hit & hit);

  void GetTimeIsolation(art::Event const & ev, std::string HitLabel, float const PeakTimeWdInt, float const PeakTimeWdExt, std::list<int> & index_list_single, std::list<int> & index_listIsolated);

  void GetListOfTimeCoincidenceHit( bool IsPDVD , bool IsPDHD , 
		                    art::Event const & ev, std::string HitLabel, const float CoincidenceWd1_l , const float CoincidenceWd1_r , const float CoincidenceWd2_l , float const CoincidenceWd2_r , const recob::Hit & HitCol,
                                                                                  std::list<geo::WireID> & WireInd1,
                                                                                  std::list<geo::WireID> & WireInd2,
                                                                                  std::list<int>   & ChannelInd1,
                                                                                  std::list<int>   & ChannelInd2,
                                                                                  std::list<float> & EInd1,
                                                                                  std::list<float> & EInd2,
                                                                                  std::list<float> & PTInd1,
                                                                                  std::list<float> & PTInd2,
										  std::list<float> & PAInd1,
                                                                                  std::list<float> & PAInd2);

  void GetListOfCrossingChannel(  bool IsPDVD , bool IsPDHD  , float Ymin , float Ymax , float Zmin , float Zmax ,
		              geo::WireID & WireCol , std::list<geo::WireID> & WireInd1 , std::list<geo::WireID> & WireInd2 ,
                              std::list<int>  & ChInd1 , std::list<float> & EInd1 , std::list<float> & YInd1 , std::list<float> & ZInd1 , std::list<int>  & ChIntersectInd1 , std::list<float> & EIntersectInd1 ,
                              std::list<int>  & ChInd2 , std::list<float> & EInd2 , std::list<float> & YInd2 , std::list<float> & ZInd2 , std::list<int>  & ChIntersectInd2 , std::list<float> & EIntersectInd2);

  void GetListOf3ViewsPoint( float pitch , float alpha ,
                             std::list<int> & ChIntersectInd1 , std::list<float> YInd1 , std::list<float> ZInd1 , std::list<float> EIntersectInd1,
                             std::list<int> & ChIntersectInd2 , std::list<float> YInd2 , std::list<float> ZInd2 , std::list<float> EIntersectInd2 ,
                             std::list<float> & listYSP       , std::list<float> & listZSP ,
                             std::list<float> & listEind1SP   , std::list<float> & listEind2SP ,
                             std::list<int> & listCh1SP       , std::list<int> & listCh2SP );

  std::vector<int> GetXYZIsolatedPoint( std::vector<float> vYPoint , std::vector<float> vZPoint , std::vector<float> vPeakTimeCol , std::vector<int> vNOF ,
                                                          float fElectronVelocity , float fTickToMus , float radiusInt , float radiusExt );
  
  bool IntersectOutsideOfTPC( float Ymin , float Ymax , float Zmin , float Zmax , double ChInd_start_y , double ChInd_start_z , double ChInd_end_y , double ChInd_end_z ,
		     	                          double ChCol_start_y , double ChCol_start_z , double ChCol_end_y , double ChCol_end_z ,
						  double& y , double& z );

  // CLUSTER FUNCTIONS
  point gen_yz(int size , std::vector<int> vIndex , std::vector<float> vY , std::vector<float> vZ , std::vector<int> vNOF);

  float dist2(point a, point b);

  void Plot(TCanvas * canvas , int event , int nbin , float ymin , float ymax , float zmin , float zmax , std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster,double RMS );

  float randf(float m);

  int nearest(point pt, point cent, int n_cluster, float *d2);

  int reallocate(point pt, std::vector<std::vector<float>> ClusterPosition , float threshold);

  float GetDist2D(float y0,float z0,float y1,float z1);

  float mean(float y,float z);
 
  void kpp(point pts, int len, point cent, int n_cent);

  std::vector<std::vector<float> > lloyd(point pts, int len, int n_cluster);

  std::vector<std::vector<float> > GetData(int len,point data);

  std::vector<int> CheckCompletude(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult );

  std::vector<int> CheckClusters(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult , float tmp);

  std::vector<Cluster> GetCluster( int n_point , int n_cluster , point p , 
                                   std::vector<float> vEInd1PointByEvent , std::vector<float> vEInd2PointByEvent , 
                                   std::vector<int> vChInd1PointByEvent , std::vector<int> vChInd2PointByEvent , 
                                   std::vector<float> vEnergyColByEvent , std::vector<float> vPeakTimeColByEvent , std::vector<int> vChannelColByEvent , 
                                   std::vector<int> vMCPDGByEvent , std::vector<int> vMCMOMpdgByEvent ,std::vector<float> vMCWeightByEvent , 
                                   std::vector<std::string> vGeneratorTagByEvent ,
                                   std::vector<float> vMCXByEvent , std::vector<float> vMCYByEvent , std::vector<float> vMCZByEvent , 
                                   std::vector<int> vNoFByEvent);
};


pdvdana::SingleHit::SingleHit(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSpacePointLabel(p.get<std::string>("SpacePointLabel","reco3d")),
    fClusterLabel(p.get<std::string>("ClusterLabel","pandora")),
    fTrackLabel(p.get<std::string>("TrackLabel","pandoraTrack")),
    fHitLabel(p.get<std::string>("HitLabel","hitpdune")),
    fG4Label(p.get<std::string>("G4Label","largeant")),
    fRDTLabel(p.get<std::string>("RDTLabel","tpcrawdecoder:daq")),

    LogLevel(p.get<int>("LogLevel")),
    fMultiplicity(p.get<int>("HitMultiplicity")),
   
    fRadiusInt(p.get<float>("RadiusInt")),
    fRadiusExt(p.get<float>("RadiusExt")),

    fCoincidenceWd1_left(p.get<float>("CoincidenceWindow1_left"  ,3)), //in mus,
    fCoincidenceWd1_right(p.get<float>("CoincidenceWindow1_right",3)), //in mus,
    fCoincidenceWd2_left(p.get<float>("CoincidenceWindow2_left"  ,3)), //in mus,
    fCoincidenceWd2_right(p.get<float>("CoincidenceWindow2_right",3)), //in mus,

    fPitch(p.get<float>("Pitch",0.5)), //cm
    fPitchMultiplier(p.get<float>("PitchMultiplier",1.2)), // 20% error    
  
    bIs3ViewsCoincidence(p.get<bool>("Is3ViewsCoincidence")),
    bHitTree(p.get<bool>("HitTree")),
    bIsPDVD(p.get<bool>("IsPDVD",false)),
    bIsPDHD(p.get<bool>("IsPDHD",false)),
    bVetoTrack(p.get<bool>("VetoTrack",false)),

    fNumberInitClusters(p.get<int>("NumberInitClusters",25)),
    fMaxSizeCluster(p.get<float>("MaxSizeCluster",0)),
    fMinSizeCluster(p.get<float>("MinSizeCluster",0)),
    fClusterSizeMulti(p.get<float>("ClusterSizeMulti",1.2)),
    fNumberConvStep(p.get<int>("NumberConvStep",300)),
    fCovering(p.get<float>("Covering",0.99))

    // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
  
  fTickTimeInMus  = sampling_rate(clockData)/1000; //caution sampling_rate(clockData) is in ns
  fCoincidenceWd1_right = fCoincidenceWd1_right/fTickTimeInMus; // in tt
  fCoincidenceWd1_left  = fCoincidenceWd1_left/fTickTimeInMus; // in tt
  fCoincidenceWd2_left  = fCoincidenceWd2_left/fTickTimeInMus; // in tt
  fCoincidenceWd2_right = fCoincidenceWd2_right/fTickTimeInMus; // in tt

  fElectronVelocity   = detProp.DriftVelocity();
 
  if(LogLevel>0){
    std::cout << " -- timing parameters -- "                             << std::endl;
    std::cout << "   sampling_rate         " << sampling_rate(clockData) << std::endl;
    std::cout << "   fTickTimeInMus:       " << fTickTimeInMus           << std::endl;
    std::cout << "   CoincidenceWd1:       " << fCoincidenceWd1_right + fCoincidenceWd1_left << std::endl;
    std::cout << "   CoincidenceWd2:       " << fCoincidenceWd2_right + fCoincidenceWd2_left << std::endl;
  } 
  //Get detector Boundaries
  unsigned fNtpcs = fGeom->NTPC();

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++)
  {
    geo::TPCID tpcid{0, t_tpc_id};
    if(fgeoXmin > fGeom->TPC(tpcid).BoundingBox().MinX()) fgeoXmin = fGeom->TPC(tpcid).BoundingBox().MinX();
    if(fgeoXmax < fGeom->TPC(tpcid).BoundingBox().MaxX()) fgeoXmax = fGeom->TPC(tpcid).BoundingBox().MaxX();
    if(fgeoYmin > fGeom->TPC(tpcid).BoundingBox().MinY()) fgeoYmin = fGeom->TPC(tpcid).BoundingBox().MinY();
    if(fgeoYmax < fGeom->TPC(tpcid).BoundingBox().MaxY()) fgeoYmax = fGeom->TPC(tpcid).BoundingBox().MaxY();
    if(fgeoZmin > fGeom->TPC(tpcid).BoundingBox().MinZ()) fgeoZmin = fGeom->TPC(tpcid).BoundingBox().MinZ();
    if(fgeoZmax < fGeom->TPC(tpcid).BoundingBox().MaxZ()) fgeoZmax = fGeom->TPC(tpcid).BoundingBox().MaxZ();
  }
  if(LogLevel>0) 
  {
    std::cout << " -- detector boundaries -- " << std::endl;
    std::cout << "  " << fgeoXmin << " < X < " << fgeoXmax << std::endl;
    std::cout << "  " << fgeoYmin << " < Y < " << fgeoYmax << std::endl;
    std::cout << "  " << fgeoZmin << " < Z < " << fgeoZmax << std::endl;
  }

  if (fMaxSizeCluster == 0) fMaxSizeCluster = fRadiusInt*1.75;
  if (fMinSizeCluster == 0) fMinSizeCluster = fRadiusInt;
}

void pdvdana::SingleHit::analyze(art::Event const& e)
{
  //Set event ID
  fEventID = e.id().event();

  // Initializing colection Hit counter
  fHitCounter = 0;

  //clock 
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  //Get DAQ Time Stamp
  if ( e.isRealData())
  {
    auto const& rdts = *e.getValidHandle<vector<raw::RDTimeStamp>>(fRDTLabel);
    CRP_T0 = rdts[0].GetTimeStamp();
  }
  else CRP_T0 = 0;
  
  //retrieve map trackID MC particle to genrator tag
  std::vector<std::string> vTrackIDToGeneratorTag;
  if (!e.isRealData()) vTrackIDToGeneratorTag = GetGeneratorTag( e , fG4Label , LogLevel , bt_serv );
 
  //retrieve hit list
  //art::InputTag hittag(fHitLabel);
  auto const HitList = e.getValidHandle<vector<recob::Hit>>(fHitLabel);
  fNHits = HitList->size();
 
  if (fNHits == 0)
  {
    if( LogLevel > 2) std::cout << " NO HIT IN EVENT " << fEventID << std::endl;

    fChannel        = -999;
    fPlane          = -999;
    fHitWidth       = -999;
    fNearOrFarToTheBeam = -999;
	  
    fEnergy         = -999;
    fPeakTime       = -999;
    fSigmaPeakTime  = -999;
    fRMS            = -999;
    fAmplitude      = -999;
    fSigmaAmplitude = -999;
    fGoodnessOfFit  = -999;
    fIntegral       = -999;
    fSigmaIntegral  = -999;
   
    fHitNumber   = -999;
    fCoincidence = -999;
    lWireInd1.clear();
    lWireInd2.clear();
    lChannelInd1.clear();
    lChannelInd1.push_back(-999);
    lChannelInd2.clear();
    lChannelInd2.push_back(-999);
    lEnergyInd1.clear();
    lEnergyInd1.push_back(-999);
    lEnergyInd2.clear();
    lEnergyInd2.push_back(-999);
    lPeakTimeInd1.clear();
    lPeakTimeInd1.push_back(-999);
    lPeakTimeInd2.clear();
    lPeakTimeInd2.push_back(-999);
    lPeakAmpInd1.clear();
    lPeakAmpInd1.push_back(-999);
    lPeakAmpInd2.clear();
    lPeakAmpInd2.push_back(-999);
    lYInd1.clear();
    lYInd1.push_back(-999);
    lZInd1.clear();
    lZInd1.push_back(-999);
    lYInd2.clear();
    lYInd2.push_back(-999);
    lZInd2.clear();
    lZInd2.push_back(-999);
    lChIntersectInd1.clear();
    lChIntersectInd1.push_back(-999);
    lChIntersectInd2.clear();
    lChIntersectInd2.push_back(-999);
    lYPoint.clear();
    lYPoint.push_back(-999);
    lZPoint.clear();
    lZPoint.push_back(-999);
    lEInd1Point.clear();
    lEInd1Point.push_back(-999);
    lEInd2Point.clear();
    lEInd2Point.push_back(-999);
    lChInd1Point.clear();
    lChInd1Point.push_back(-999);
    lChInd2Point.clear();
    lChInd2Point.push_back(-999);
    lEIntersectInd1.clear();
    lEIntersectInd1.push_back(-999);
    lEIntersectInd2.clear();
    lEIntersectInd2.push_back(-999);
    vMCPart_pdgCode.clear();
    vMCPart_pdgCode.push_back(-999);
    vMCPart_mother.clear();
    vMCPart_mother.push_back(-999);
    vMCPart_weight.clear();
    vMCPart_weight.push_back(-999);
    vGenerator_tag.clear();
    vGenerator_tag.push_back("nothing");
    vMCPart_Endx.clear();
    vMCPart_Endx.push_back(-9999);
    vMCPart_Endy.clear();
    vMCPart_Endy.push_back(-9999);
    vMCPart_Endz.clear();
    vMCPart_Endz.push_back(-9999);
    vMCPart_Startx.clear();
    vMCPart_Startx.push_back(-9999);
    vMCPart_Starty.clear();
    vMCPart_Starty.push_back(-9999);
    vMCPart_Startz.clear();
    vMCPart_Startz.push_back(-9999);

    if (bHitTree) tHitTree->Fill();

    NCluster = 0;
    vZCluster.push_back(-999);
    vYCluster.push_back(-999);
    vEColCluster.push_back(-999);
    vEInd1Cluster.push_back(-999);
    vEInd2Cluster.push_back(-999);
    vPTCluster.push_back(-999);

    vNoFCluster.push_back(-999);
    vNPointCluster.push_back(-999);
    vNColCluster.push_back(-999);
    vNInd1Cluster.push_back(-999);
    vNInd2Cluster.push_back(-999);

    vMCPDGCluster.push_back( -999 );
    vMCMOMpdgCluster.push_back( -999 );
    vMCWeightCluster.push_back( -999 );
    vMCGenTagCluster.push_back( "nothing" );
    vMCXCluster.push_back(-9999);
    vMCYCluster.push_back(-9999);
    vMCZCluster.push_back(-9999);

    vVetoTrackStartX.push_back(-9999.0);
    vVetoTrackStartY.push_back(-9999.0);
    vVetoTrackStartZ.push_back(-9999.0);
    vVetoTrackEndX.push_back(-9999.0);
    vVetoTrackEndY.push_back(-9999.0);
    vVetoTrackEndZ.push_back(-9999.0);

    tClusterTree->Fill();

    return;
  } 

  if( !lSingleIndex.empty()   ) lSingleIndex.clear();
  if( !lIsolatedIndex.empty() ) lIsolatedIndex.clear();

  // Single Isolated concideration
  GetSingle(   e, fHitLabel, lSingleIndex, fMultiplicity );

  float fPeakTimeWdInt = ( fRadiusInt / fElectronVelocity )/fTickTimeInMus;
  float fPeakTimeWdExt = ( fRadiusExt / fElectronVelocity )/fTickTimeInMus;

  if ( LogLevel > 2)  std::cout << "HIT BASED ANALYSIS --> Isolation, Spatialisation etc..." << std::endl;

  GetTimeIsolation( e, fHitLabel, fPeakTimeWdInt, fPeakTimeWdExt, lSingleIndex, lIsolatedIndex);

  if ( !vYPointByEvent.empty() )       vYPointByEvent.clear();
  if ( !vZPointByEvent.empty() )       vZPointByEvent.clear();
  if ( !vEInd1PointByEvent.empty() )   vEInd1PointByEvent.clear();
  if ( !vEInd2PointByEvent.empty() )   vEInd2PointByEvent.clear();
  if ( !vChInd1PointByEvent.empty() )  vChInd1PointByEvent.clear();
  if ( !vChInd2PointByEvent.empty() )  vChInd2PointByEvent.clear();
  if ( !vNoFByEvent.empty() )          vNoFByEvent.clear();
  if ( !vEnergyColByEvent.empty() )    vEnergyColByEvent.clear();
  if ( !vPeakTimeColByEvent.empty() )  vPeakTimeColByEvent.clear();
  if ( !vChannelColByEvent.empty() )   vChannelColByEvent.clear();
  if ( !vMCPDGByEvent.empty() )        vMCPDGByEvent.clear();
  if ( !vMCMOMpdgByEvent.empty() )     vMCMOMpdgByEvent.clear();
  if ( !vMCWeightByEvent.empty() )     vMCWeightByEvent.clear();
  if ( !vMCXByEvent.empty() )          vMCXByEvent.clear();
  if ( !vMCYByEvent.empty() )          vMCYByEvent.clear();
  if ( !vMCZByEvent.empty() )          vMCZByEvent.clear();
  if ( !vGeneratorTagByEvent.empty() ) vGeneratorTagByEvent.clear();

  for(int index =0 ; index<fNHits; index++)
  {

    const recob::Hit& hit = HitList->at(index);
    if (!e.isRealData()) 
    {

      vMCPart_pdgCode.clear();
      vMCPart_mother.clear();
      vMCPart_motherPdg.clear();
      vMCPart_weight.clear();
      vGenerator_tag.clear();
      vMCPart_Endx.clear();
      vMCPart_Endy.clear();
      vMCPart_Endz.clear();
      vMCPart_Startx.clear();
      vMCPart_Starty.clear();
      vMCPart_Startz.clear();

      float sumMCEnergy = 0 ;

      using weightedMCPair = std::pair<const simb::MCParticle*, float>;
      std::vector<weightedMCPair> tempMCPair;

      for (const auto & ide : bt_serv->HitToEveTrackIDEs(clockData, hit))
      {
        const simb::MCParticle* curr_part = pi_serv->TrackIdToParticle_P(ide.trackID);
        tempMCPair.push_back(std::make_pair(curr_part,ide.energy));
        sumMCEnergy += ide.energy;
      }

      
      std::sort(tempMCPair.begin(), tempMCPair.end(), [](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});
      
      for (weightedMCPair& p : tempMCPair ) 
      {   
        vMCPart_pdgCode.push_back( (p.first)->PdgCode()   );
        vMCPart_mother.push_back(  (p.first)->Mother()    );
        vMCPart_weight.push_back(  (p.second)/sumMCEnergy );
	vGenerator_tag.push_back( vTrackIDToGeneratorTag[(p.first)->TrackId()] );

	vMCPart_Endx.push_back( (float) (p.first)->EndX() );
        vMCPart_Endy.push_back( (float) (p.first)->EndY() );
        vMCPart_Endz.push_back( (float) (p.first)->EndZ() );

	vMCPart_Startx.push_back( (float) (p.first)->Vx(0) );
	vMCPart_Starty.push_back( (float) (p.first)->Vy(0) );
	vMCPart_Startz.push_back( (float) (p.first)->Vz(0) );

        if ( (p.first)->Mother() == 0 )
        {
          vMCPart_motherPdg.push_back( 0 ); //primary particle
          continue;
        }
        const simb::MCParticle* curr_part_mom = pi_serv->TrackIdToParticle_P((p.first)->Mother());
        vMCPart_motherPdg.push_back( curr_part_mom->PdgCode() );
      }
      
      
      //if( LogLevel > 2) print(vMCPart_Endz); 
    }// end if event != real data
    else
    {
      vMCPart_pdgCode.clear();
      vMCPart_pdgCode.push_back(-999);
      vMCPart_mother.clear();
      vMCPart_mother.push_back(-999);
      vMCPart_motherPdg.clear();
      vMCPart_motherPdg.push_back(-999);
      vMCPart_weight.clear();
      vMCPart_weight.push_back(-999);
      vGenerator_tag.clear();
      vGenerator_tag.push_back("data");
      vMCPart_Endx.clear();
      vMCPart_Endx.push_back(-9999);
      vMCPart_Endy.clear();
      vMCPart_Endy.push_back(-9999);
      vMCPart_Endz.clear();
      vMCPart_Endz.push_back(-9999);
      vMCPart_Startx.clear();
      vMCPart_Startx.push_back(-9999);
      vMCPart_Starty.clear();
      vMCPart_Starty.push_back(-9999);
      vMCPart_Startz.clear();
      vMCPart_Startz.push_back(-9999);
    }

    fWire           = hit.WireID();
    fChannel        = hit.Channel();
    fPlane          = hit.WireID().Plane;
    fHitWidth       = hit.EndTick() - hit.StartTick();
	 
    fNearOrFarToTheBeam = NearOrFar(bIsPDVD , bIsPDHD , hit);

    fEnergy         = hit.ROISummedADC();///fADCtoEl;
    fPeakTime       = hit.PeakTime();//*ftick_in_mus;
    fSigmaPeakTime  = hit.SigmaPeakTime();//*ftick_in_mus;
    fRMS            = hit.RMS();
    fAmplitude      = hit.PeakAmplitude();
    fSigmaAmplitude = hit.SigmaPeakAmplitude();
    fGoodnessOfFit  = hit.GoodnessOfFit();
    fIntegral       = hit.Integral();//fADCtoEl;
    fSigmaIntegral  = hit.SigmaIntegral();

    if( !lWireInd1.empty()         ) lWireInd1.clear();
    if( !lWireInd2.empty()         ) lWireInd2.clear();
    if( !lChannelInd1.empty()      ) lChannelInd1.clear();
    if( !lChannelInd2.empty()      ) lChannelInd2.clear();
    if( !lEnergyInd1.empty()       ) lEnergyInd1.clear();
    if( !lEnergyInd2.empty()       ) lEnergyInd2.clear();
    if( !lPeakTimeInd1.empty()     ) lPeakTimeInd1.clear();
    if( !lPeakTimeInd2.empty()     ) lPeakTimeInd2.clear();
    if( !lPeakAmpInd1.empty()      ) lPeakAmpInd1.clear();
    if( !lPeakAmpInd2.empty()      ) lPeakAmpInd2.clear();
    if( !lYInd1.empty()            ) lYInd1.clear();
    if( !lZInd1.empty()            ) lZInd1.clear();
    if( !lYInd2.empty()            ) lYInd2.clear();
    if( !lZInd2.empty()            ) lZInd2.clear();
    if( !lChIntersectInd1.empty()  ) lChIntersectInd1.clear();
    if( !lChIntersectInd2.empty()  ) lChIntersectInd2.clear();
    if( !lYPoint.empty()           ) lYPoint.clear();
    if( !lZPoint.empty()           ) lZPoint.clear();
    if( !lEInd1Point.empty()       ) lEInd1Point.clear();
    if( !lEInd2Point.empty()       ) lEInd2Point.clear();
    if( !lChInd1Point.empty()      ) lChInd1Point.clear();
    if( !lChInd2Point.empty()      ) lChInd2Point.clear();

    if (fPlane == 2)
    {
      fHitNumber = fHitCounter;
      ++fHitCounter;

      // Coincidence research

      GetListOfTimeCoincidenceHit( bIsPDVD, bIsPDHD , e, fHitLabel, fCoincidenceWd1_left, fCoincidenceWd1_right ,fCoincidenceWd2_left, fCoincidenceWd2_right , hit, lWireInd1, lWireInd2, lChannelInd1, lChannelInd2, lEnergyInd1, lEnergyInd2, lPeakTimeInd1, lPeakTimeInd2, lPeakAmpInd1, lPeakAmpInd2);
      //GetListOfTimeCoincidenceHit( e, fHitLabel, 10 , 0 , hit, lWireInd1, lWireInd2, lChannelInd1, lChannelInd2, lEnergyInd1, lEnergyInd2, lPeakTimeInd1, lPeakTimeInd2);

      fCoincidence = 0;
      if ( !lWireInd1.empty() || !lWireInd2.empty() ) fCoincidence += 1;
      if ( !lWireInd1.empty() && !lWireInd2.empty() ) fCoincidence += 1;

      if ( fCoincidence > 0 )
      {
        GetListOfCrossingChannel(  bIsPDVD, bIsPDHD , fgeoYmin , fgeoYmax , fgeoZmin , fgeoZmax , fWire , lWireInd1 , lWireInd2 , 
			          lChannelInd1 , lEnergyInd1 , lYInd1 , lZInd1 , lChIntersectInd1 , lEIntersectInd1 , 
				  lChannelInd2 , lEnergyInd2 , lYInd2 , lZInd2 , lChIntersectInd2 , lEIntersectInd2); 
        if ( bIs3ViewsCoincidence ) GetListOf3ViewsPoint( fPitch , fPitchMultiplier , lChIntersectInd1 , lYInd1 , lZInd1 , lEIntersectInd1 , lChIntersectInd2 , lYInd2 , lZInd2 , lEIntersectInd2 , lYPoint , lZPoint , lEInd1Point , lEInd2Point , lChInd1Point , lChInd2Point);
	else
	{
	  //induction 1
	  lYPoint.insert( lYPoint.end() , lYInd1.begin() , lYInd1.end() );
	  lZPoint.insert( lZPoint.end() , lZInd1.begin() , lZInd1.end() );
	  lEInd1Point.insert( lEInd1Point.end() , lEIntersectInd1.begin() , lEIntersectInd1.end() );
	  lChInd1Point.insert( lChInd1Point.end() , lChIntersectInd1.begin() , lChIntersectInd1.end() );
	  int N1 = lYInd1.size();
	  lEInd2Point.insert( lEInd2Point.end() , N1 , 0 );
	  lChInd2Point.insert( lChInd2Point.end() , N1 , -1 );

	  //induction 2
	  lYPoint.insert( lYPoint.end() , lYInd2.begin() , lYInd2.end() );
	  lZPoint.insert( lZPoint.end() , lZInd2.begin() , lZInd2.end() );
	  lEInd2Point.insert( lEInd2Point.end() , lEIntersectInd2.begin() , lEIntersectInd2.end() );
	  lChInd2Point.insert( lChInd2Point.end() , lChIntersectInd2.begin() , lChIntersectInd2.end() );
	  int N2 = lYInd2.size();
	  lEInd1Point.insert( lEInd1Point.end() , N2 , 0 );
	  lChInd1Point.insert( lChInd1Point.end() , N2 , -1 );
	}
	  
      }
    }

    if (fPlane != 2)
    {
      fHitNumber   = -1;
      fCoincidence = -999;

      lWireInd1.clear();
      lWireInd2.clear();
      lChannelInd1.clear();
      lChannelInd1.push_back(-999);
      lChannelInd2.clear();
      lChannelInd2.push_back(-999);
      lEnergyInd1.clear();
      lEnergyInd1.push_back(-999);
      lEnergyInd2.clear();
      lEnergyInd2.push_back(-999);
      lPeakTimeInd1.clear();
      lPeakTimeInd1.push_back(-999);
      lPeakTimeInd2.clear();
      lPeakTimeInd2.push_back(-999);
      lPeakAmpInd1.clear();
      lPeakAmpInd1.push_back(-999);
      lPeakAmpInd2.clear();
      lPeakAmpInd2.push_back(-999);
      lYInd1.clear();
      lYInd1.push_back(-999);
      lZInd1.clear();
      lZInd1.push_back(-999);
      lYInd2.clear();
      lYInd2.push_back(-999);
      lZInd2.clear();
      lZInd2.push_back(-999);
      lChIntersectInd1.clear();
      lChIntersectInd1.push_back(-999);
      lChIntersectInd2.clear();
      lChIntersectInd2.push_back(-999);
      lYPoint.clear();
      lYPoint.push_back(-999);
      lZPoint.clear();
      lZPoint.push_back(-999);
      lEInd1Point.clear();
      lEInd1Point.push_back(-999);
      lEInd2Point.clear();
      lEInd2Point.push_back(-999);
      lChInd1Point.clear();
      lChInd1Point.push_back(-999);
      lChInd2Point.clear();
      lChInd2Point.push_back(-999);
      lEIntersectInd1.clear();
      lEIntersectInd1.push_back(-999);
      lEIntersectInd2.clear();
      lEIntersectInd2.push_back(-999);
    }
  
    fTimeIsolation = 0;
    if ( Inside(index, lSingleIndex)   ) fTimeIsolation+=1;
    if ( Inside(index, lIsolatedIndex) ) fTimeIsolation+=1;
  
    //Filling Hit tree
    if (bHitTree) tHitTree->Fill();

    // Retreiving hit info by event

    std::list<int>::iterator ch1  = lChInd1Point.begin();
    std::list<int>::iterator ch2  = lChInd2Point.begin();
    std::list<float>::iterator e1 = lEInd1Point.begin();
    std::list<float>::iterator e2 = lEInd2Point.begin();
    std::list<float>::iterator z  = lZPoint.begin();

    
    for ( auto const y : lYPoint)
    {
      if(( y == -999) || (*z == -999) || (*e1 == -999) || (*e2 == -999) || (*ch1 == -999) || (*ch2 == -999))
      {
        z++;
        ch1++;
        ch2++;
        e1++;
        e2++;
	continue;
      }

      vYPointByEvent.push_back( y );
      vZPointByEvent.push_back( *z );
      vEInd1PointByEvent.push_back( *e1 );
      vEInd2PointByEvent.push_back( *e2 );
      vChInd1PointByEvent.push_back( *ch1 );
      vChInd2PointByEvent.push_back( *ch2 );
      vNoFByEvent.push_back( fNearOrFarToTheBeam );
      vEnergyColByEvent.push_back( fEnergy   );
      vPeakTimeColByEvent.push_back(  fPeakTime );
      vChannelColByEvent.push_back( fChannel );

      if ( (int) vMCPart_pdgCode.size() > 1 ) fAmbiguousHit += 1 ;

      // take first origin MC truth for now
      vMCPDGByEvent.push_back( vMCPart_pdgCode[0] );
      vMCMOMpdgByEvent.push_back( vMCPart_motherPdg[0] );
      vMCWeightByEvent.push_back( vMCPart_weight[0] );
      vMCXByEvent.push_back( vMCPart_Endx[0] );
      vMCYByEvent.push_back( vMCPart_Endy[0] );
      vMCZByEvent.push_back( vMCPart_Endz[0] );

      vGeneratorTagByEvent.push_back( vGenerator_tag[0] );

      z++;
      ch1++;
      ch2++;
      e1++;
      e2++;

    }

    // setting hit variables values to initiale values
    fChannel        = -999;
    fPlane          = -999;
    fHitWidth       = -999;
    fHitNumber      = -999;
    fCoincidence    = -999;
    fTimeIsolation  = -999;
    fAmbiguousHit   = 0;
    fNearOrFarToTheBeam = -999;
	  
    fEnergy         = -999.;
    fPeakTime       = -999.;
    fSigmaPeakTime  = -999.;
    fRMS            = -999.;
    fAmplitude      = -999.;
    fSigmaAmplitude = -999.;
    fGoodnessOfFit  = -999.;
    fIntegral       = -999.;
    fSigmaIntegral  = -999.;

    lSingleIndex.clear();
    lIsolatedIndex.clear();

    lWireInd1.clear();
    lWireInd2.clear();
    lChannelInd1.clear();
    lChannelInd2.clear();
    lEnergyInd1.clear();
    lEnergyInd2.clear();
    lPeakTimeInd1.clear();
    lPeakTimeInd2.clear();
    lPeakAmpInd1.clear();
    lPeakAmpInd2.clear();
    lYInd1.clear();
    lZInd1.clear();
    lYInd2.clear();
    lZInd2.clear();
    lChIntersectInd1.clear();
    lChIntersectInd2.clear();
    lYPoint.clear();
    lZPoint.clear();
    lEInd1Point.clear();
    lEInd2Point.clear();
    lChInd1Point.clear();
    lChInd2Point.clear();
    lEIntersectInd1.clear();
    lEIntersectInd2.clear();
    vMCPart_pdgCode.clear();
    vMCPart_mother.clear();
    vMCPart_motherPdg.clear();
    vMCPart_weight.clear();
    vMCPart_Endx.clear();
    vMCPart_Endy.clear();
    vMCPart_Endz.clear();
    vMCPart_Startx.clear();
    vMCPart_Starty.clear();
    vMCPart_Startz.clear();

  }// end hit loop

  if ( LogLevel > 4) std::cout << "THERE ARE " << fAmbiguousHit << " HIT(s) WITH AMBIGUOUS ORIGIN " << std::endl;

  std::vector<int> vIso = GetXYZIsolatedPoint( vYPointByEvent , vZPointByEvent , vPeakTimeColByEvent , vNoFByEvent , fElectronVelocity , fTickTimeInMus , fRadiusInt , fRadiusExt );
  int PTSIsolated = (int) vIso.size();

  if (PTSIsolated == 0)
  {
    if ( LogLevel > 0) std::cout << "EXEPTION ERROR : THERE IS NO ISOLATED POINT IN EVENT " << fEventID << std::endl;

    NCluster = 0;
    vZCluster.push_back(-999);
    vYCluster.push_back(-999);
    vEColCluster.push_back(-999);
    vEInd1Cluster.push_back(-999);
    vEInd2Cluster.push_back(-999);
    vPTCluster.push_back(-999);

    vNoFCluster.push_back(-999);
    vNPointCluster.push_back(-999);
    vNColCluster.push_back(-999);
    vNInd1Cluster.push_back(-999);
    vNInd2Cluster.push_back(-999);

    vMCPDGCluster.push_back( -999 );
    vMCMOMpdgCluster.push_back( -999 );
    vMCWeightCluster.push_back( -999 );
    vMCGenTagCluster.push_back( "nothing" ); 
    vMCXCluster.push_back(-9999);
    vMCYCluster.push_back(-9999);
    vMCZCluster.push_back(-9999);

    vVetoTrackStartX.push_back(-9999.0);
    vVetoTrackStartZ.push_back(-9999.0);
    vVetoTrackStartY.push_back(-9999.0);
    vVetoTrackEndX.push_back(-9999.0);
    vVetoTrackEndZ.push_back(-9999.0);
    vVetoTrackEndY.push_back(-9999.0);

    tClusterTree->Fill();
    return;
  }

  if( LogLevel > 5)
  {
  std::cout << " THERE ARE " << vYPointByEvent.size() << " POINTS IN EVENT " << fEventID << std::endl;
  std::cout << " THERE ARE " << PTSIsolated << " ISOLATED POINTS IN EVENT " << fEventID << std::endl;
  }

  point v;
  v = gen_yz( PTSIsolated , vIso , vYPointByEvent , vZPointByEvent , vNoFByEvent );


  std::vector<std::vector<float> > dataPos = GetData(PTSIsolated,v);
  std::vector<std::vector<float> > clustersPos;

  std::vector<int> vchecks(2,0);

  int K = fNumberInitClusters;

  int check = 0;
  float threshold = fMinSizeCluster;

  fAsConverged = false;
  for( int i = 0 ; i < fNumberConvStep ; i++ )
  {
    if ( LogLevel > 0) printf("%d %d %d ",i,K,check);

    clustersPos = lloyd(v, PTSIsolated, K);
    vchecks = CheckClusters( dataPos , clustersPos , threshold , fClusterSizeMulti, fCovering);
    check = vchecks[0];

    if(check == 1)  
    {
      fAsConverged = true;	     
      break;
    }
    if(check == 2)
    {
      if (threshold < fMaxSizeCluster)
      {
        threshold *= fClusterSizeMulti;
        if ( LogLevel > 3) printf("Threshold increased \n");
        K = vchecks[1];
      }
      else
      {
        K = vchecks[1] + 1;
        if ( LogLevel > 2) printf("Threshold Max reached \n");
      }
    }
    else K = vchecks[1] + 5;
  }

  if ( LogLevel > 2) printf("Data size : %lu x %lu \n",dataPos.size(),dataPos[0].size());
  if ( LogLevel > 2) printf("Cluster size : %lu x %lu \n",clustersPos.size(),clustersPos[0].size());

  if ( LogLevel > 0) printf("Data clustering ended successfully \n");

  int j = 0;
  point p;
  for (j = 0, p = v; j < PTSIsolated ; j++, p++)
  {
    p->group = reallocate( p , clustersPos , threshold);
  }
  if (LogLevel > 3) std::cout<< " reallocation done " << std::endl;


  if ( LogLevel > 2)
  {
    art::ServiceHandle<art::TFileService> tfs;
    TString canvasName = Form("clusterisation_%i", fEventID);
    TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);
    Plot(canvas , fEventID ,2000 , fgeoYmin-50 , fgeoYmax+50 , -2*fgeoZmax-50 , fgeoZmax+50, dataPos ,clustersPos,threshold );
    canvas->Write(canvasName);
  } 

  std::vector<Cluster> vCluster = GetCluster( PTSIsolated , clustersPos[0].size() , v , vEInd1PointByEvent , vEInd2PointByEvent , vChInd1PointByEvent , vChInd2PointByEvent , vEnergyColByEvent , vPeakTimeColByEvent , vChannelColByEvent , vMCPDGByEvent , vMCMOMpdgByEvent , vMCWeightByEvent , vGeneratorTagByEvent , vMCXByEvent , vMCYByEvent , vMCZByEvent , vNoFByEvent);
  NCluster = vCluster.size();

  if ( !vYCluster.empty() )        vYCluster.clear();
  if ( !vZCluster.empty() )        vZCluster.clear();
  if ( !vEColCluster.empty() )     vEColCluster.clear();
  if ( !vEInd1Cluster.empty() )    vEInd1Cluster.clear();
  if ( !vEInd2Cluster.empty() )    vEInd2Cluster.clear();
  if ( !vPTCluster.empty() )       vPTCluster.clear();
  if ( !vNoFCluster.empty() )      vNoFCluster.clear();
  if ( !vNPointCluster.empty() )   vNPointCluster.clear();
  if ( !vNColCluster.empty() )     vNColCluster.clear();
  if ( !vNInd1Cluster.empty() )    vNInd1Cluster.clear();
  if ( !vNInd2Cluster.empty() )    vNInd2Cluster.clear();
  if ( !vMCMOMpdgCluster.empty() ) vMCMOMpdgCluster.clear();
  if ( !vMCPDGCluster.empty() )    vMCPDGCluster.clear();
  if ( !vMCWeightCluster.empty() ) vMCWeightCluster.clear();
  if ( !vMCXCluster.empty() )      vMCXCluster.clear();
  if ( !vMCYCluster.empty() )      vMCYCluster.clear();
  if ( !vMCZCluster.empty() )      vMCZCluster.clear();
  if ( !vMCGenTagCluster.empty() ) vMCGenTagCluster.clear();
  
  int NEmpty_cluster = 0;
  for( int j = 0 ; j < NCluster ; j++ )
  {
    // Filling tree variables
    if ( vCluster[j].Npoint == 0) 
    {
      NEmpty_cluster += 1; 
      if( LogLevel > 4) std::cout << "no point in cluster " << j << std::endl;
      continue;
    }

    vYCluster.push_back(vCluster[j].y);
    vZCluster.push_back(vCluster[j].z);
    vEColCluster.push_back(vCluster[j].ECol);
    vEInd1Cluster.push_back(vCluster[j].EInd1);
    vEInd2Cluster.push_back(vCluster[j].EInd2);
    vPTCluster.push_back(vCluster[j].PeakTime);
	  
    vNoFCluster.push_back(vCluster[j].NOF);
    vNPointCluster.push_back(vCluster[j].Npoint);
    vNColCluster.push_back(vCluster[j].NCol);
    vNInd1Cluster.push_back(vCluster[j].NInd1);
    vNInd2Cluster.push_back(vCluster[j].NInd2);
    
    vMCMOMpdgCluster.insert( vMCMOMpdgCluster.end() , (vCluster[j].vMCMOMpdg).begin() , (vCluster[j].vMCMOMpdg).end() );
    vMCPDGCluster.insert(    vMCPDGCluster.end()    , (vCluster[j].vMCPDG   ).begin() , (vCluster[j].vMCPDG   ).end() );
    vMCWeightCluster.insert( vMCWeightCluster.end() , (vCluster[j].vMCWEI   ).begin() , (vCluster[j].vMCWEI   ).end() );

    vMCXCluster.insert( vMCXCluster.end() , (vCluster[j].vMCX).begin() , (vCluster[j].vMCX).end() );
    vMCYCluster.insert( vMCYCluster.end() , (vCluster[j].vMCY).begin() , (vCluster[j].vMCY).end() );
    vMCZCluster.insert( vMCZCluster.end() , (vCluster[j].vMCZ).begin() , (vCluster[j].vMCZ).end() );

    vMCGenTagCluster.insert( vMCGenTagCluster.end()    , (vCluster[j].vMCGenTag).begin() , (vCluster[j].vMCGenTag).end() );
  }

  NCluster = NCluster - NEmpty_cluster;
  // veto part
  if (bVetoTrack)
  {
    auto const tracklist = e.getValidHandle<vector<recob::Track>>(fTrackLabel);

    vVetoTrackStartX.clear();
    vVetoTrackStartY.clear();
    vVetoTrackStartZ.clear();
    vVetoTrackEndX.clear();
    vVetoTrackEndY.clear();
    vVetoTrackEndZ.clear();

    for (unsigned itrk = 0; itrk < tracklist->size(); ++itrk) 
    {
      const recob::Track& track = tracklist->at(itrk);
      vVetoTrackStartX.push_back( track.Start().X() );
      vVetoTrackStartY.push_back( track.Start().Y() );
      vVetoTrackStartZ.push_back( track.Start().Z() );
      vVetoTrackEndX.push_back( track.End().X() );
      vVetoTrackEndY.push_back( track.End().Y() );
      vVetoTrackEndZ.push_back( track.End().Z() );
    }
  }
  tClusterTree->Fill();
   

  // setting variables values to initiale values
  fChannel        = -999;
  fPlane          = -999;
  fHitWidth       = -999;
  fHitNumber      = -999;
  fCoincidence    = -999;
  fTimeIsolation  = -999;
  fAmbiguousHit   = 0;
  fNearOrFarToTheBeam = -999;
	
  fEnergy         = -999.;
  fPeakTime       = -999.;
  fSigmaPeakTime  = -999.;
  fRMS            = -999.;
  fAmplitude      = -999.;
  fSigmaAmplitude = -999.;
  fGoodnessOfFit  = -999.;
  fIntegral       = -999.;
  fSigmaIntegral  = -999.;
 
  lSingleIndex.clear();
  lIsolatedIndex.clear();

  lWireInd1.clear();
  lWireInd2.clear();
  lChannelInd1.clear();
  lChannelInd2.clear();
  lEnergyInd1.clear();
  lEnergyInd2.clear();
  lPeakTimeInd1.clear();
  lPeakTimeInd2.clear();
  lPeakAmpInd1.clear();
  lPeakAmpInd2.clear();
  lYInd1.clear();
  lZInd1.clear();
  lYInd2.clear();
  lZInd2.clear();
  lChIntersectInd1.clear();
  lChIntersectInd2.clear();
  lYPoint.clear();
  lZPoint.clear();
  lEInd1Point.clear();
  lEInd2Point.clear();
  lChInd1Point.clear();
  lChInd2Point.clear();
  vMCPart_pdgCode.clear();
  vMCPart_mother.clear();
  vMCPart_motherPdg.clear();
  vMCPart_weight.clear();
  vMCPart_Endx.clear();
  vMCPart_Endy.clear();
  vMCPart_Endz.clear();
  vMCPart_Startx.clear();
  vMCPart_Starty.clear();
  vMCPart_Startz.clear();
  lEIntersectInd1.clear();
  lEIntersectInd2.clear();

  NCluster = -999;

  vYCluster.clear();
  vZCluster.clear();
  vEColCluster.clear();
  vEInd1Cluster.clear();
  vEInd2Cluster.clear();
  vPTCluster.clear();

  vNoFCluster.clear();

  vNPointCluster.clear();
  vNColCluster.clear();
  vNInd1Cluster.clear();
  vNInd2Cluster.clear();

  vMCMOMpdgCluster.clear();
  vMCMOMpdgCluster.clear();
  vMCWeightCluster.clear();
  vMCGenTagCluster.clear();
  vMCXCluster.clear();
  vMCYCluster.clear();
  vMCZCluster.clear();

  vYPointByEvent.clear();
  vZPointByEvent.clear();

  vNoFByEvent.clear();

  vEnergyColByEvent.clear();
  vEInd1PointByEvent.clear();
  vEInd2PointByEvent.clear();

  vPeakTimeColByEvent.clear();

  vChannelColByEvent.clear();
  vChInd1PointByEvent.clear();
  vChInd2PointByEvent.clear();

  vMCMOMpdgByEvent.clear();
  vMCPDGByEvent.clear();
  vMCWeightByEvent.clear();
  vGeneratorTagByEvent.clear();
  vMCXByEvent.clear();
  vMCYByEvent.clear();
  vMCZByEvent.clear();

  vVetoTrackStartX.clear();
  vVetoTrackStartY.clear();
  vVetoTrackStartZ.clear();
  vVetoTrackEndX.clear();
  vVetoTrackEndY.clear();
  vVetoTrackEndZ.clear();

  // Implementation of required member function here.
}

void pdvdana::SingleHit::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  
  // HIT TREE

  if (bHitTree) 
  {
    tHitTree = tfs->make<TTree>("hitTree", "Output Tree");

    //add branches
    tHitTree->Branch("eventID"       , &fEventID       );
    tHitTree->Branch("CRP_T0"        , &CRP_T0         );
    tHitTree->Branch("wireID"        , &fWire          );
    tHitTree->Branch("hitNumber"     , &fHitNumber     );
    tHitTree->Branch("plane"         , &fPlane         );
    tHitTree->Branch("channel"       , &fChannel       );
    tHitTree->Branch("timeIsolation" , &fTimeIsolation );
    tHitTree->Branch("coincidence"   , &fCoincidence   );
    tHitTree->Branch("NearOrFarToTheBeam" , &fNearOrFarToTheBeam );
	
    tHitTree->Branch("energy"           , &fEnergy        );
    tHitTree->Branch("peakTime"         , &fPeakTime      );
    tHitTree->Branch("hitWidth"         , &fHitWidth      );
    tHitTree->Branch("sigmaPeakTime"    , &fSigmaPeakTime );
    tHitTree->Branch("amplitudePeaktime", &fAmplitude     );
    tHitTree->Branch("sigmaAmplitude"   , &fSigmaAmplitude);
    tHitTree->Branch("rms"              , &fRMS           );
    tHitTree->Branch("goodnessOfFit"    , &fGoodnessOfFit );
    tHitTree->Branch("integral"         , &fIntegral      );
    tHitTree->Branch("sigmaIntegral"    , &fSigmaIntegral );

    tHitTree->Branch("listChannelInd1"        , &lChannelInd1     );
    tHitTree->Branch("listEnergyInd1"         , &lEnergyInd1      );
    tHitTree->Branch("listPeakTimeInd1"       , &lPeakTimeInd1    );
    tHitTree->Branch("listPeakAmpInd1"        , &lPeakAmpInd1     );
    tHitTree->Branch("listYIntersecPointInd1" , &lYInd1           );
    tHitTree->Branch("listZIntersecPointInd1" , &lZInd1           );
    tHitTree->Branch("listChIntersecInd1"     , &lChIntersectInd1 );

    tHitTree->Branch("listChannelInd2"        , &lChannelInd2     );
    tHitTree->Branch("listEnergyInd2"         , &lEnergyInd2      );
    tHitTree->Branch("listPeakTimeInd2"       , &lPeakTimeInd2    );
    tHitTree->Branch("listPeakAmpInd2"        , &lPeakAmpInd2     );
    tHitTree->Branch("listYIntersecPointInd2" , &lYInd2           );
    tHitTree->Branch("listZIntersecPointInd2" , &lZInd2           );
    tHitTree->Branch("listChIntersecInd2"     , &lChIntersectInd2 );

    tHitTree->Branch("listOfYPoint"           , &lYPoint      );
    tHitTree->Branch("listOfZPoint"           , &lZPoint      );
    tHitTree->Branch("listEInd1Point"         , &lEInd1Point  );
    tHitTree->Branch("listEInd2Point"         , &lEInd2Point  );
    tHitTree->Branch("listChInd1Point"        , &lChInd1Point );
    tHitTree->Branch("listChInd2Point"        , &lChInd2Point );

    tHitTree->Branch("vectorOFMCParticlePDG"     , &vMCPart_pdgCode   );
    tHitTree->Branch("vectorOFMCParticleMom"     , &vMCPart_mother    );
    tHitTree->Branch("vectorOFMCParticleMomPDG"  , &vMCPart_motherPdg );
    tHitTree->Branch("vectorOFMCParticleWeight"  , &vMCPart_weight    );
    tHitTree->Branch("vectorOFMCParticleEX"      , &vMCPart_Endx      );
    tHitTree->Branch("vectorOFMCParticleEY"      , &vMCPart_Endy      );
    tHitTree->Branch("vectorOFMCParticleEZ"      , &vMCPart_Endz      );
    tHitTree->Branch("vectorOFMCParticleSX"      , &vMCPart_Startx    );
    tHitTree->Branch("vectorOFMCParticleSY"      , &vMCPart_Starty    );
    tHitTree->Branch("vectorOFMCParticleSZ"      , &vMCPart_Startz    );
    tHitTree->Branch("vectorOFGeneratorTag"      , &vGenerator_tag    );

  }

  // CLUSTER TREE
  tClusterTree = tfs->make<TTree>("ClusterTree","ClusterTree");

  tClusterTree->Branch("AsConverged"        , &fAsConverged     );
  tClusterTree->Branch("eventID"            , &fEventID         );
  tClusterTree->Branch("CRP_T0"             , &CRP_T0           );
  tClusterTree->Branch("NumberOfCluster"    , &NCluster         );
  tClusterTree->Branch("Z"                  , &vZCluster        );
  tClusterTree->Branch("Y"                  , &vYCluster        );
  tClusterTree->Branch("EnergyCollection"   , &vEColCluster     );
  tClusterTree->Branch("NearOrFarToTheBeam" , &vNoFCluster      );
  tClusterTree->Branch("EnergyPlane0"       , &vEInd1Cluster    );
  tClusterTree->Branch("EnergyPlane1"       , &vEInd2Cluster    );
  tClusterTree->Branch("NumberOfPoint"      , &vNPointCluster   );
  tClusterTree->Branch("NumberOfCollection" , &vNColCluster     );
  tClusterTree->Branch("NumberOfPlane0"     , &vNInd1Cluster    );
  tClusterTree->Branch("NumberOfPlane1"     , &vNInd2Cluster    );
  tClusterTree->Branch("PeakTime"           , &vPTCluster       );
  tClusterTree->Branch("MCParticlePDG"      , &vMCPDGCluster    );
  tClusterTree->Branch("MCParticleMOMPdg"   , &vMCMOMpdgCluster );
  tClusterTree->Branch("MCParticleWeight"   , &vMCWeightCluster );
  tClusterTree->Branch("MCParticleX"        , &vMCXCluster      );
  tClusterTree->Branch("MCParticleY"        , &vMCYCluster      );
  tClusterTree->Branch("MCParticleZ"        , &vMCZCluster      );
  tClusterTree->Branch("HitGenerationTag"   , &vMCGenTagCluster );
  if (LogLevel > 2)
  {
    tClusterTree->Branch("Y_isolated_point" , &vY_point          );
    tClusterTree->Branch("Z_isolated_point" , &vZ_point          );
  }
  if (bVetoTrack)
  {
    tClusterTree->Branch("VetoTrackStartX"    , &vVetoTrackStartX );
    tClusterTree->Branch("VetoTrackStartY"    , &vVetoTrackStartY );
    tClusterTree->Branch("VetoTrackStartZ"    , &vVetoTrackStartZ );
    tClusterTree->Branch("VetoTrackEndX"      , &vVetoTrackEndX   );
    tClusterTree->Branch("VetoTrackEndY"      , &vVetoTrackEndY   );
    tClusterTree->Branch("VetoTrackEndZ"      , &vVetoTrackEndZ   );
  }
}


void pdvdana::SingleHit::endJob()
{
  // Implementation of optional member function here.
}

void pdvdana::SingleHit::print(std::vector<std::vector<int>> vv )
{
  int size1 =  vv.size();
  for( int k = 0 ; k < size1 ; k++ )
  {
    int size2 = vv[k].size();
    std::cout << "[ " ;
    for( int i = 0 ; i < size2 ; i++ )
    {
      std::cout << vv[k][i] << " , ";
    }
    std::cout << " ]" << std::endl;
  }
}
void pdvdana::SingleHit::print(std::vector<float> v )
{
  int size = v.size();
  std::cout << "[ ";
  for( int k = 0 ; k < size ; k++)
  {
    std::cout << v[k] << " , ";
  }
  std::cout << " ]" << std::endl;
}

float pdvdana::SingleHit::GetDist( float x0 , float y0 , float z0 , float x1 , float y1 , float z1 )
{
  float x = x0-x1;
  float z = z0-z1;
  float y = y0-y1;
  return sqrt(x*x+z*z+y*y);
}

bool pdvdana::SingleHit::Inside( int k , std::list<int> list){
  return (std::find(list.begin(), list.end(), k) != list.end());
}

bool pdvdana::SingleHit::AllSame( std::vector<int> v)
{
  if (v.size() == 0) return false;
  bool allsame = true;
  for( int i = 0; i < (int) v.size() ; i++)
  {
    if( v[i] == v[0] ) continue;
    allsame = false;
    break;
  }
  return allsame;
}
int pdvdana::SingleHit::NearOrFar( bool IsPDVD , bool IsPDHD , const recob::Hit & hit)
{
  if (bIsPDHD)
  {
    if ( (hit.WireID().TPC == 2)|| (hit.WireID().TPC == 6) || (hit.WireID().TPC == 3) || (hit.WireID().TPC == 7) ) return  -1; // 3 and 7 are dumie TPC
    if ( (hit.WireID().TPC == 1)|| (hit.WireID().TPC == 5) || (hit.WireID().TPC == 0) || (hit.WireID().TPC == 4) ) return  1;  // 0 and 4 are dumie TPC
  }
  if (bIsPDVD)
  {
    if (hit.WireID().TPC <= 7 ) return -1;
    if (hit.WireID().TPC  > 7 ) return 1;
  }
  return -999;
}

void pdvdana::SingleHit::GetSingle(art::Event const & ev, std::string HitLabel, std::list<int> & index_list_single , int const Multiplicity)
{
  auto const hitlist = ev.getValidHandle<vector<recob::Hit>>(HitLabel);
  recob::Hit hit = hitlist->at(0);

  for (int i=0, sz=hitlist->size(); i!=sz; ++i)
  {
    hit = hitlist->at(i);
    if(hit.Multiplicity() > Multiplicity) continue;
    index_list_single.push_back(i);    
  } 
}

void pdvdana::SingleHit::GetTimeIsolation(art::Event const & ev, std::string HitLabel, float const PeakTimeWdInt, float const PeakTimeWdExt,  std::list<int> & index_list_single, std::list<int> & index_listIsolated)
{
  auto const hitlist = ev.getValidHandle<vector<recob::Hit>>(HitLabel);
  index_listIsolated = index_list_single;

  recob::Hit hit = hitlist->at(0);
  float PeakTime = -999;

  float PeakTimeSingle     = -999;
  float PeakTimeMinInt     = -999;
  float PeakTimeMaxInt     = -999;
  float PeakTimeMinExt     = -999;
  float PeakTimeMaxExt     = -999;
   
  if(not( index_listIsolated.empty()))
  {
    for (int i=0, sz=hitlist->size(); i!=sz; ++i)
    {
      hit      = hitlist->at(i);
      PeakTime = hit.PeakTime();

      PeakTimeSingle     = -999;
      PeakTimeMinInt     = -999;
      PeakTimeMaxInt     = -999;
      PeakTimeMinExt     = -999;
      PeakTimeMaxExt     = -999;

      std::list<int>::iterator elem = index_listIsolated.begin();

      while( elem != index_listIsolated.end())
      {
        PeakTimeSingle  = (hitlist->at(*elem)).PeakTime();
      
        PeakTimeMinInt     = PeakTimeSingle - PeakTimeWdInt;
        PeakTimeMaxInt     = PeakTimeSingle + PeakTimeWdInt;

        PeakTimeMinExt     = PeakTimeMinInt - PeakTimeWdExt;
        PeakTimeMaxExt     = PeakTimeMaxInt + PeakTimeWdExt;
        
        if ( ((PeakTime >= PeakTimeMinExt)&&(PeakTime < PeakTimeMinInt)) || ((PeakTime > PeakTimeMaxInt)&&(PeakTime <= PeakTimeMaxExt)) )
        {
          if (i != *elem) //normally always true now
          {
            elem = index_listIsolated.erase(elem);
          }
        }
        ++elem;
        
        PeakTimeSingle  = -999;
        PeakTimeMinInt     = -999;
        PeakTimeMaxInt     = -999;
        PeakTimeMinExt     = -999;
        PeakTimeMaxExt     = -999;
      }
      
      PeakTime        = -999;
    }
  }
}

void pdvdana::SingleHit::GetListOfTimeCoincidenceHit( bool IsPDVD , bool IsPDHD , art::Event const & ev, std::string HitLabel, float const CoincidenceWd1_l , float const CoincidenceWd1_r, float const CoincidenceWd2_l , float const CoincidenceWd2_r, const recob::Hit & HitCol, 
                                                                                  std::list<geo::WireID> & WireInd1,
                                                                                  std::list<geo::WireID> & WireInd2,
                                                                                  std::list<int>   & ChannelInd1,
                                                                                  std::list<int>   & ChannelInd2,
                                                                                  std::list<float> & EInd1,
                                                                                  std::list<float> & EInd2,
                                                                                  std::list<float> & PTInd1,
                                                                                  std::list<float> & PTInd2,
										  std::list<float> & PAInd1,
                                                                                  std::list<float> & PAInd2)
{
  auto const hitlist = ev.getValidHandle<vector<recob::Hit>>(HitLabel);

  recob::Hit hit = hitlist->at(0);

  float PeakTimeCol    = HitCol.PeakTime();
  //float RMSPeakTimeCol = HitCol.RMS();

  float EndTime1   = PeakTimeCol + CoincidenceWd1_r;
  float EndTime2   = PeakTimeCol + CoincidenceWd2_r;
  float StartTime1 = PeakTimeCol - CoincidenceWd1_l;
  float StartTime2 = PeakTimeCol - CoincidenceWd2_l;

  float PeakTime = -999;
  int   Plane    = -999;
  int NoFCol = NearOrFar(IsPDVD,IsPDHD,HitCol);
  int NoF = -4;

  for (int i=0, sz=hitlist->size(); i!=sz; ++i)
  { 
    hit   = hitlist->at(i);
    Plane = hit.WireID().Plane;
    if (Plane == 2) continue;

    NoF = NearOrFar(IsPDVD,IsPDHD,hit);
    if (NoF != NoFCol) continue;

    PeakTime = hit.PeakTime();
    if (Plane == 0)
    {
      if ((PeakTime < StartTime1)||(PeakTime > EndTime1)) continue;

      WireInd1.push_back(hit.WireID());
      ChannelInd1.push_back(hit.Channel());
      EInd1.push_back(hit.ROISummedADC());
      PTInd1.push_back(PeakTime);
      PAInd1.push_back(hit.PeakAmplitude());
      continue;
    }
    if (Plane == 1)
    {
      if ((PeakTime < StartTime2)||(PeakTime > EndTime2)) continue;

      WireInd2.push_back(hit.WireID());
      ChannelInd2.push_back(hit.Channel());
      EInd2.push_back(hit.ROISummedADC());
      PTInd2.push_back(PeakTime);
      PAInd2.push_back(hit.PeakAmplitude());
    }
  }
}

bool pdvdana::SingleHit::IntersectOutsideOfTPC( float Ymin , float Ymax , float Zmin , float Zmax ,
		                                double ChInd_start_y , double ChInd_start_z ,
					       	double ChInd_end_y ,double ChInd_end_z ,
						double ChCol_start_y , double ChCol_start_z ,
					    	double ChCol_end_y , double ChCol_end_z ,
		     				double& y , double& z )
{
  // Equation from http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection

  double const denom = (ChInd_start_y - ChInd_end_y) * (ChCol_start_z - ChCol_end_z) - (ChInd_start_z - ChInd_end_z) * (ChCol_start_y - ChCol_end_y);

  if (denom == 0) return false;

  double const A = (ChInd_start_y * ChInd_end_z - ChInd_start_z * ChInd_end_y) / denom;
  double const B = (ChCol_start_y * ChCol_end_z - ChCol_start_z * ChCol_end_y) / denom;

  y = (ChCol_start_y - ChCol_end_y) * A - (ChInd_start_y - ChInd_end_y) * B;
  z = (ChCol_start_z - ChCol_end_z) * A - (ChInd_start_z - ChInd_end_z) * B;

  bool drap = ( y > Ymin ) && ( y < Ymax ) && ( z > Zmin ) && ( z < Zmax ) ;
  
  return drap;
}

void pdvdana::SingleHit::GetListOfCrossingChannel(   bool IsPDVD , bool IsPDHD , float Ymin , float Ymax , float Zmin , float Zmax ,
		                                    geo::WireID & WireCol , std::list<geo::WireID> & WireInd1 , std::list<geo::WireID> & WireInd2 , 
						    std::list<int>  & ChInd1 , std::list<float> & EInd1 , std::list<float> & YInd1 , std::list<float> & ZInd1 , 
						    std::list<int>  & ChIntersectInd1 , std::list<float> & EIntersectInd1 ,
                                                    std::list<int>  & ChInd2 , std::list<float> & EInd2 , std::list<float> & YInd2 , std::list<float> & ZInd2 , 
						    std::list<int>  & ChIntersectInd2 , std::list<float> & EIntersectInd2)
{
  geo::Point_t point = geo::Point_t(-999,-999,-999);
  bool drap;

  double y = -999. , z = -999.;
  //auto const wcol = fGeom->WireEndPoints(WireCol);
  auto const [wcolstart, wcolend] = fWireReadout.WireEndPoints(WireCol);

  ChIntersectInd1.clear();
  EIntersectInd1.clear();
  ChIntersectInd2.clear();
  EIntersectInd2.clear();
  
  std::list<int>::iterator ch1  = ChInd1.begin();
  std::list<float>::iterator e1   = EInd1.begin();

  for (auto const elementInd1 : WireInd1)
  {
    if (WireCol.TPC != elementInd1.TPC)
    {
    
      if (bIsPDVD)
      {  
        //auto const wind1 = fGeom->WireEndPoints(elementInd1);
        auto const [wind1start, wind1end] = fWireReadout.WireEndPoints(elementInd1);
        bool flag = IntersectOutsideOfTPC( Ymin , Ymax , Zmin ,Zmax , wind1start.Y() , wind1start.Z() , wind1end.Y() , wind1end.Z() , wcolstart.Y() , wcolstart.Z() , wcolend.Y() , wcolend.Z() , y , z);
        if (flag)
        {
	  YInd1.push_back(y);
          ZInd1.push_back(z);
          EIntersectInd1.push_back(*e1);
	  ChIntersectInd1.push_back(*ch1);
          ++e1; 
	  ++ch1;
	  continue;
        }
      }
      ++e1;
      ++ch1;
      continue ;
    }

    //drap = fGeom->WireIDsIntersect( WireCol , elementInd1 , point);
    drap = fWireReadout.WireIDsIntersect( WireCol , elementInd1 , point);
    if ( drap )
    {
      YInd1.push_back(point.Y());
      ZInd1.push_back(point.Z());
      ChIntersectInd1.push_back(*ch1);
      EIntersectInd1.push_back(*e1);
    }
    ++e1;
    ++ch1;
  }
  std::list<int>::iterator ch2  = ChInd2.begin();
  std::list<float>::iterator e2   = EInd2.begin();

  for (auto const elementInd2 : WireInd2)
  { 
    if (WireCol.TPC != elementInd2.TPC ) 
    { 

      if (bIsPDVD)
      {
        //auto const wind2 = fGeom->WireEndPoints(elementInd2);
        auto const [wind2start, wind2end] = fWireReadout.WireEndPoints(elementInd2);
        bool flag = IntersectOutsideOfTPC( Ymin , Ymax , Zmin ,Zmax , wind2start.Y() , wind2start.Z() , wind2end.Y() , wind2end.Z() , wcolstart.Y() , wcolstart.Z() , wcolend.Y() , wcolend.Z() , y , z);
        if (flag)
        {
          YInd2.push_back(y);
          ZInd2.push_back(z);
          ChIntersectInd2.push_back(*ch2);
	  EIntersectInd2.push_back(*e2);
	  ++e2;
          ++ch2;
          continue;
        }
      }
      ++e2;
      ++ch2; 
      continue ;
    }
    //drap = fGeom->WireIDsIntersect( WireCol , elementInd2 , point);
    drap = fWireReadout.WireIDsIntersect( WireCol , elementInd2 , point);
    if ( drap ) 
    { 
      YInd2.push_back(point.Y());
      ZInd2.push_back(point.Z());
      ChIntersectInd2.push_back(*ch2);
      EIntersectInd2.push_back(*e2);
    }
    ++e2;
    ++ch2;
  }
}

void pdvdana::SingleHit::GetListOf3ViewsPoint( float pitch , float alpha , 
                                        std::list<int> & ChIntersectInd1 , std::list<float> YInd1 , std::list<float> ZInd1 , std::list<float> EIntersectInd1, 
                                        std::list<int> & ChIntersectInd2 , std::list<float> YInd2 , std::list<float> ZInd2 , std::list<float> EIntersectInd2 , 
                                        std::list<float> & listYSP       , std::list<float> & listZSP , 
                                        std::list<float> & listEind1SP   , std::list<float> & listEind2SP , 
                                        std::list<int> & listCh1SP       , std::list<int> & listCh2SP)
{

  std::list<int>::iterator  ch1t = ChIntersectInd1.begin();
  std::list<float>::iterator z1t = ZInd1.begin();
  std::list<float>::iterator e1t = EIntersectInd1.begin();

  float dy, dz, dr;

  for( auto const yind1 : YInd1)
  {
    std::list<int>::iterator  ch2t = ChIntersectInd2.begin();
    std::list<float>::iterator z2t = ZInd2.begin();
    std::list<float>::iterator e2t = EIntersectInd2.begin();

    for ( auto const yind2 : YInd2)
    {
      dy = yind1 - yind2;
      dz = *z1t - *z2t  ;
      dr = TMath::Sqrt( dy*dy + dz*dz );

      if ( dr <= pitch*alpha )
      {
        float y = ( (*e1t)*(yind1) + (*e2t)*(yind2) )/( *e1t + *e2t );
        float z = ( (*e1t)*(*z1t) + (*e2t)*(*z2t) )/( *e1t + *e2t );

        listYSP.push_back( y );     
        listZSP.push_back( z );
        listEind1SP.push_back( *e1t );
        listEind2SP.push_back( *e2t );
        listCh1SP.push_back( *ch1t );
        listCh2SP.push_back( *ch2t );
      }
      ++e2t;
      ++z2t;
      ++ch2t;
    }
    ++e1t;
    ++z1t;
    ++ch1t;
  }
}


std::vector<int> pdvdana::SingleHit::GetXYZIsolatedPoint( std::vector<float> vYPoint , std::vector<float> vZPoint , std::vector<float> vPeakTimeCol , std::vector<int> vNOF ,
					                  float fElectronVelocity , float fTickToMus , float radiusInt , float radiusExt )
{

  if (vYPoint.size() != vZPoint.size())  throw std::invalid_argument( "BIG PROBLEM" );

  std::vector<int> vIso;
  vIso.clear();

  if (radiusInt >= radiusExt) 
  {
    if (LogLevel > 4) std::cout << " Size of DONUT NON PHYSIQUE " << std::endl;   
    return vIso;
  }
  
  int npoint = vYPoint.size();
  std::vector<int> vIsIsolated( npoint , -1 );
  vIso.resize( npoint, -1);

  float zIs = 0;
  float yIs = 0;
  float xIs = 0;
  float xDiff = 0;
  float yDiff = 0;
  float zDiff = 0;

  int cellX = 0;
  int cellY = 0;
  int cellZ = 0;

  bool IsIsolated = true;

  float distSq = 0;
  int iso_count = 0;

  float electronDriftScale = fElectronVelocity * fTickToMus;
  float radiusIntSq = radiusInt*radiusInt;
  float radiusExtSq = radiusExt*radiusExt;

  int nof = 0;

  float gridCellSize = 2*radiusExt; // The size of each grid cell is based on the external radius
  std::unordered_map<std::tuple<int, int, int>, std::vector<int>, GridHasher> grid;

  // Populate the grid with points
  for (int k = 0; k < npoint; k++) 
  {	  
    xIs = vPeakTimeCol[k] * electronDriftScale;
    yIs = vYPoint[k];
    zIs = vZPoint[k];

    // Skip invalid points
    if ((yIs == -999) || (zIs == -999)) 
    {
      vIsIsolated[k] = 0;
      continue;
    }

    // Determine which cell the point belongs to
    cellX = (int)(xIs / gridCellSize);
    cellY = (int)(yIs / gridCellSize);
    cellZ = (int)(zIs / gridCellSize);

    // Store the point in the appropriate grid cell
    grid[std::make_tuple(cellX, cellY, cellZ)].push_back(k);
  }

  for( int k = 0 ; k<npoint ; k++)
  {
    if (vIsIsolated[k] == 0)
    {
      continue;
    }

    zIs = vZPoint[k];
    yIs = vYPoint[k];
    xIs = vPeakTimeCol[k]*electronDriftScale;

    IsIsolated = true;

    nof = vNOF[k];
    if (( yIs == -999) || (zIs == -999))
    {
      vIsIsolated[k] = 0;
      continue;
    }

    // Determine the grid cell of the point
    cellX = (int)(xIs / gridCellSize);
    cellY = (int)(yIs / gridCellSize);
    cellZ = (int)(zIs / gridCellSize);

    for (int dx = -1; dx <= 1; dx++) 
    {
      for (int dy = -1; dy <= 1; dy++) 
      {
        for (int dz = -1; dz <= 1; dz++) 
	{
          auto neighborCell = std::make_tuple(cellX + dx, cellY + dy, cellZ + dz);

          // If the neighboring cell contains points
          if (grid.find(neighborCell) != grid.end()) 
	  {
            for (int i : grid[neighborCell]) // i indice in the vector of indices of point in neighborCell
	    {
              if (i == k) continue; // Skip comparing the point with itself
              if (nof != vNOF[i]) continue; 

              xDiff = xIs - vPeakTimeCol[i] * electronDriftScale;
              yDiff = yIs - vYPoint[i];
              zDiff = zIs - vZPoint[i];

              distSq = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff;

              if (distSq > radiusIntSq && distSq < radiusExtSq) 
	      {
                vIsIsolated[k] = 0; // Mark both points as non-isolated
                vIsIsolated[i] = 0;
                IsIsolated = false;
                break;
              }
            }
          }
          if (!IsIsolated) break;
        }
        if (!IsIsolated) break;
      }
      if (!IsIsolated) break;
    }

    if (IsIsolated)
    {
      vIsIsolated[k] = 1;
      vIso[iso_count] = k;
      iso_count++;
    }
  }//end point loop

  vIso.resize(iso_count);
  return vIso;
}

///////////////////////////////////////////////////////////////////////////////////
// CLUSTER FUNCTIONS

float pdvdana::SingleHit::dist2(point a, point b)
{
    float z = a->z - b->z, y = a->y - b->y;
    return z*z + y*y;
}

float pdvdana::SingleHit::randf(float m)
{
    return m * rand() / (RAND_MAX - 1.);
}

point pdvdana::SingleHit::gen_yz(int size , std::vector<int> vIndex , std::vector<float> vY , std::vector<float> vZ , std::vector<int> vNOF)
{
  int i = 0;
  point p, pt = (point) malloc(sizeof(point_t) * size);

  for (p = pt + size; p-- > pt;)
  {
    p->y = vY[vIndex[i]];
    p->index = vIndex[i];

    if (vNOF[vIndex[i]] == -1)
    {
      p->z = -1*vZ[vIndex[i]] - fgeoZmax;
    }
    else p->z = vZ[vIndex[i]];
    i++;
  }

  return pt;
}

void pdvdana::SingleHit::Plot(TCanvas* canvas ,int event , int nbin , float ymin , float ymax , float zmin , float zmax , std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster,double RMS ){

    TGraph* grData = new TGraph();
    TGraph* grCluster = new TGraph();

    TH2F* Background = new TH2F(Form("Background_%d",event),"",nbin,ymin,ymax,nbin,zmin,zmax);
    Background->SetXTitle("y[cm]");
    Background->SetYTitle("z[cm]");
    Background->SetStats(0);
    std::vector<TEllipse*> Ell;

    for(int i = 0;i< (int) data[0].size();i++) grData->SetPoint(i,data[1][i],data[0][i]);

    for(int i = 0;i< (int) cluster[0].size();i++){
        //printf("Final cluster %i position : %f %f \n",i+1,cluster[1][i],cluster[0][i]);

        grCluster->SetPoint(i,cluster[1][i],cluster[0][i]);

        TEllipse* ell = new TEllipse();
        ell->SetX1(cluster[1][i]);
        ell->SetY1(cluster[0][i]);
        ell->SetR1(RMS);
        ell->SetR2(RMS);
        ell->SetLineColor(i);
        ell->SetFillStyle(3001);
        ell->SetFillColor(i%9+1);

        Ell.push_back(ell);
    }

    grData->SetMarkerStyle(50);
    grData->SetMarkerSize(0.9);
    grData->SetMarkerColor(kGray+2);

    grCluster->SetMarkerStyle(3);
    grCluster->SetMarkerColor(kBlack);

    canvas->cd();
    Background->Draw();
    grData->Draw("P SAME");
    grCluster->Draw("P SAME");

    for(int i = 0;i< (int) Ell.size();i++) 
    {
      canvas->cd();	 
      Ell[i]->Draw("same");
    }
}

int pdvdana::SingleHit::nearest(point pt, point cent, int n_cluster, float *d2)
{
    int i = 0;
    int  min_i = 0;
    point c;
    float d = 0.0;
    float  min_d = 0.0;

#       define for_n for (c = cent, i = 0; i < n_cluster; i++, c++)
    for_n {
        min_d = HUGE_VAL;
        min_i = pt->group;
        for_n {
            if (min_d > (d = dist2(c, pt))) {
                min_d = d; min_i = i;
            }
        }
    }
    if (d2) *d2 = min_d;
    return min_i;
}

int pdvdana::SingleHit::reallocate(point pt, std::vector<std::vector<float>> ClusterPosition , float threshold)
{
    int  min_i = pt->group;
    float  min_d = HUGE_VAL;

    for( int k = 0 ; k < (int) ClusterPosition[0].size() ; k++) 
    {

      float dist = sqrt(GetDist2D(pt->z,pt->y,ClusterPosition[0][k],ClusterPosition[1][k]));
      if (min_d > dist ) 
      {
         min_d = dist; 
	 min_i = k;
      }
     
    }
    if ( min_d > threshold ) return -1;
    return min_i;
}


float pdvdana::SingleHit::GetDist2D(float y0,float z0,float y1,float z1){
    float z = z0-z1;
    float y = y0-y1;
    return z*z+y*y;
}

float pdvdana::SingleHit::mean(float y,float z){
    return (z+y)/2.;
}

void pdvdana::SingleHit::kpp(point pts, int len, point cent, int n_cent)
{
#       define for_len for (j = 0, p = pts; j < len; j++, p++)
    int j;
    int n_cluster;
    float sum, *d = (float*)malloc(sizeof(float) * len);

    point p;
    cent[0] = pts[ rand() % len ];
    for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
        sum = 0;
        for_len {
            nearest(p, cent, n_cluster, d + j);
            sum += d[j];
        }
        sum = randf(sum);
        for_len {
            if ((sum -= d[j]) > 0) continue;
            cent[n_cluster] = pts[j];
            break;
        }
    }
    for_len p->group = nearest(p, cent, n_cluster, 0);
    free(d);
}

std::vector<std::vector<float>> pdvdana::SingleHit::lloyd(point pts, int len, int n_cluster)
{
    int i, j, min_i;
    int changed;

    point cent = (point)malloc(sizeof(point_t) * n_cluster), p, c;

    /* assign init grouping randomly */
    //for_len p->group = j % n_cluster;

    /* or call k++ init */
    kpp(pts, len, cent, n_cluster);

    do {
        /* group element for centroids are used as counters */
        for_n { c->group = 0; c->z = c->y = 0; }
        for_len {
            c = cent + p->group;
            c->group++;
            c->z += p->z; c->y += p->y;
        }
        for_n { c->z /= c->group; c->y /= c->group; }

        changed = 0;
        /* fInd closest centroid of each point */
        for_len {
            min_i = nearest(p, cent, n_cluster, 0);
            if (min_i != p->group) {
                changed++;
                p->group = min_i;
            }
        }
    } while (changed > (len >> 10)); /* stop when 99.9% of points are good */

    for_n { c->group = i; }

    std::vector<std::vector<float> > clusterPos;
    std::vector<float> clusterPosY,clusterPosZ;

    point result;

    for(i = 0, result = cent; i < n_cluster; i++, result++) {
        clusterPosZ.push_back(result->z);
        clusterPosY.push_back(result->y);
    }
    clusterPos.push_back(clusterPosZ);
    clusterPos.push_back(clusterPosY);

    return clusterPos;
}


std::vector<std::vector<float>> pdvdana::SingleHit::GetData(int len,point data){

    std::vector<std::vector<float> > dataPos;
    std::vector<float> dataPosZ,dataPosY;


  for(int i = 0; i < len; i++, data++) {
        dataPosZ.push_back(data->z);
        dataPosY.push_back(data->y);
    }
    dataPos.push_back(dataPosZ);
    dataPos.push_back(dataPosY);

    return dataPos;
}

std::vector<int> pdvdana::SingleHit::CheckCompletude(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult )
{
    int Npts = data[0].size();
    int Ncls = cluster[0].size();

    int Nin = 0, Nin2 = 0, Nout = 0;
    float dist;
    std::vector<int> IDin(Npts,0);

    for(int i = 0 ; i < Ncls ; i++ )
    {
        Nout = 0;
        for(int j = 0;j<Npts;j++)
        {
            if(IDin[j] == 0)
            {
                dist = sqrt(GetDist2D(data[0][j],data[1][j],cluster[0][i],cluster[1][i]));
                // printf("Distance to cluster : %f %f \n",i,Nin,Nout);
                if(dist <= RMS){ Nin++; IDin[j] = 1;}
                else if(dist <= mult*RMS )
                {
                  Nin2++;
                  IDin[j] = 1;
                }
                else{ Nout++; }
            }
        }
    }

    //std::cout << "!!!!!!!!!!!!!!! COMPLETUDE : " << Nin << " " << Nin2 << " " << Nout << std::endl;
    std::vector<int> N(3,0);
    N[0] = Nin;
    N[1] = Nin2;
    N[2] = Nout;

    return N;
}

std::vector<int> pdvdana::SingleHit::CheckClusters(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult , float tmp)
{

    std::vector<int> v(2,0);
    int Npts = data[0].size();
    int Ncls = cluster[0].size();
    //std::cout << " nombre of point for check cluster : " << Npts << std::endl;

    int Nin = 0, Nin2 = 0, Nout = 0;

    std::vector<int> NComp(3,0);
    NComp  =  CheckCompletude( data, cluster , RMS , mult);
    Nin  = NComp[0];
    Nin2 = NComp[1];
    Nout = NComp[2];

    //std::cout<< " Nin " << Nin << " Nin2 " << Nin2 << " Nout " << Nout << " Npst " << Npts << " Ncls " << Ncls << std::endl;

    if ( LogLevel > 0) printf("Counting : %.03f %.03f %.03f sum : %.02f \n",float(Nin)/float(Npts),float(Nin2)/float(Npts),float(Nout)/float(Npts),float(Nin+Nin2+Nout)/float(Npts));


    if((float(Nin+Nin2)/float(Npts) > tmp -0.04) && float(Nin)/float(Npts) < tmp)
    {
      v[0] = 2;
      v[1] = cluster[0].size();
      return v;
    }
    else if(float(Nin)/float(Npts) < tmp)
    {
      v[0] = 0;
      v[1] = cluster[0].size();
      return v;
    }
    //if(float(Nin)/float(Npts) < 0.99) return 0;


    float dist;
    std::vector<std::vector<int>> IDoverlap( Ncls , std::vector<int> (Ncls ,0));
    int overlap = 0;

    for(int i = 0 ; i < Ncls ; i++)
    {
        for(int j = i ; j < Ncls ; j++)
        {
            if(j > i)
            {
                dist = sqrt(GetDist2D(cluster[0][j],cluster[1][j],cluster[0][i],cluster[1][i]));
                if(dist < 2.*RMS)
                {
                    IDoverlap[i][j] = 1;
                    overlap++;
                }
            }
        }
    }
    int overlap_counter = 0;

    vector<float> newclusterZ, newclusterY;
    float meanZ, meanY;
    while( (overlap_counter<10)&&(overlap>0) )
    {
      newclusterZ.clear(); 
      newclusterY.clear();
      meanZ = 0;
      meanY = 0;

      for(int i = 0;i<Ncls;i++)
      {
        meanZ = 0.;
        meanY = 0.;
        overlap = 0;

        for(int j = i;j<Ncls;j++)
        {
          if(IDoverlap[i][j] == 1 && cluster[0][j] != -999)
          {
            meanZ += mean(cluster[0][i],cluster[0][j]);
            meanY += mean(cluster[1][i],cluster[1][j]);
            cluster[0][j] = -999;
            cluster[1][j] = -999;
            overlap++;
          }
        }
        if(overlap == 0 && cluster[0][i] != -999)
        {
          newclusterZ.push_back(cluster[0][i]);
          newclusterY.push_back(cluster[1][i]);
        }
        else if(cluster[0][i] != -999)
        {
          newclusterZ.push_back(meanZ/float(overlap));
          newclusterY.push_back(meanY/float(overlap));
        }
      }

      cluster.clear();
      cluster.push_back(newclusterZ);
      cluster.push_back(newclusterY);

      if ( LogLevel > 3) printf("%lu clusters has been removed at iteration %d \n",Ncls-cluster[0].size(),overlap_counter);

      dist = 0;
      Ncls = cluster[0].size();
      IDoverlap.clear();
      overlap = 0;

      for(int i = 0 ; i < Ncls ; i++)
      {
	std::vector<int> v(Ncls , 0); 
        for(int j = i ; j < Ncls ; j++)
        {
          if(j > i)
          {
            dist = sqrt(GetDist2D(cluster[0][j],cluster[1][j],cluster[0][i],cluster[1][i]));
            if(dist < 2.*RMS)
            {
              v[j] = 1;
              overlap++;
            }
          }
        }
	IDoverlap.push_back(v);
      }

      overlap_counter++;
    }

    NComp = CheckCompletude( data, cluster , RMS , mult );
    Nin  = NComp[0];
    Nin2 = NComp[1];
    Nout = NComp[2];

    if ( LogLevel > 3) printf("Counting : %.03f %.03f %.03f sum : %.02f \n",float(Nin)/float(Npts),float(Nin2)/float(Npts),float(Nout)/float(Npts),float(Nin+Nin2+Nout)/float(Npts));

    if(float(Nin+Nin2)/float(Npts) > tmp && float(Nin)/float(Npts) < tmp)
    {
      v[0] = 2;
      v[1] = cluster[0].size();
      return v;
    }
    else if(float(Nin)/float(Npts) < tmp)
    {
      v[0] = 0;
      v[1] = cluster[0].size();
      return v;
    }

    v[0] = 1;
    v[1] = cluster[0].size();
    return v;

}

std::vector<Cluster> pdvdana::SingleHit::GetCluster( int n_point , int n_cluster , point p , 
                                                     std::vector<float> vEInd1PointByEvent , std::vector<float> vEInd2PointByEvent , 
                                                     std::vector<int> vChInd1PointByEvent , std::vector<int> vChInd2PointByEvent ,  
                                                     std::vector<float> vEnergyColByEvent , std::vector<float> vPeakTimeColByEvent , std::vector<int> vChannelColByEvent ,  
                                                     std::vector<int> vMCPDGByEvent , std::vector<int> vMCMOMpdgByEvent ,std::vector<float> vMCWeightByEvent ,  
                                                     std::vector<std::string> vGeneratorTagByEvent ,  
                                                     std::vector<float> vMCXByEvent ,  std::vector<float> vMCYByEvent , std::vector<float> vMCZByEvent , 
                                                     std::vector<int> vNoFByEvent )
{

  int k;
  point vp;
  TempCluster NullCluster;
  NullCluster.Sumz = 0;
  NullCluster.Sumy = 0;
  NullCluster.Npoint = 0;
  NullCluster.ECol = 0;
  NullCluster.PeakTime = 0;
  NullCluster.NCol = 0;
  NullCluster.vNOF.clear();
  NullCluster.vMCPDG.clear();
  NullCluster.vMCMOMpdg.clear();
  NullCluster.vMCWEI.clear();
  NullCluster.vMCGenTag.clear();
  NullCluster.vMCX.clear();
  NullCluster.vMCY.clear();
  NullCluster.vMCZ.clear();
  NullCluster.lChannelCol.clear();
  NullCluster.lChannelInd1.clear();
  NullCluster.EInd1 = 0;
  NullCluster.NInd1 = 0;
  NullCluster.lChannelInd2.clear();
  NullCluster.EInd2 = 0;
  NullCluster.NInd2= 0;
    
  std::vector<TempCluster> vTempCluster(n_cluster, NullCluster);
 
  int out = 0;

  if( LogLevel > 5) std::cout << "there are " << n_cluster << " clusters to check" << std::endl;
  for( k = 0 , vp=p ; k<n_point ; k++ , vp++)
  {

    int index = vp->index;
    int ClusterID   = vp->group;


    int NoF         = vNoFByEvent[index];
    int ChannelCol  = vChannelColByEvent[index];
    int ChannelInd1 = vChInd1PointByEvent[index];
    int ChannelInd2 = vChInd2PointByEvent[index];
    float ECol      = vEnergyColByEvent[index];
    float PeakTime  = vPeakTimeColByEvent[index];

    int MCPart_pdg = vMCPDGByEvent[index];
    int MCPart_mompdg = vMCMOMpdgByEvent[index];
    float MCPart_weight = vMCWeightByEvent[index];
    float MCPart_x = vMCXByEvent[index];
    float MCPart_y = vMCYByEvent[index];
    float MCPart_z = vMCZByEvent[index];

    std::string Generator_tag = vGeneratorTagByEvent[index];

    if (ClusterID >=  n_cluster)
    {
      out++;
      continue;
    }

    if (ClusterID ==  -1)
    {
      out++;
      continue;
    }
    vTempCluster[ClusterID].Npoint += 1;

    if ( !Inside( ChannelCol , vTempCluster[ClusterID].lChannelCol ) )
    {
      if (NoF == -1) vTempCluster[ClusterID].Sumz += ECol * ( -1*(vp->z)-fgeoZmax );
      else vTempCluster[ClusterID].Sumz += ECol * ( vp->z );

      vTempCluster[ClusterID].Sumy += ECol * ( vp->y );
      vTempCluster[ClusterID].ECol += ECol;
      vTempCluster[ClusterID].PeakTime += PeakTime;
      vTempCluster[ClusterID].NCol += 1;
      (vTempCluster[ClusterID].vMCPDG).push_back( MCPart_pdg );
      (vTempCluster[ClusterID].vMCMOMpdg).push_back( MCPart_mompdg );
      (vTempCluster[ClusterID].vMCWEI).push_back( MCPart_weight );
      (vTempCluster[ClusterID].vMCGenTag).push_back( Generator_tag );
      (vTempCluster[ClusterID].vMCX).push_back( MCPart_x );
      (vTempCluster[ClusterID].vMCY).push_back( MCPart_y );
      (vTempCluster[ClusterID].vMCZ).push_back( MCPart_z );
      (vTempCluster[ClusterID].lChannelCol).push_back( ChannelCol );
      (vTempCluster[ClusterID].vNOF).push_back( NoF );
    }
    if ( (ChannelInd1 != -1) && (!Inside( ChannelInd1 , vTempCluster[ClusterID].lChannelInd1 ) ) )
    {
      vTempCluster[ClusterID].EInd1 += vEInd1PointByEvent[index];
      vTempCluster[ClusterID].NInd1 += 1;
      (vTempCluster[ClusterID].lChannelInd1).push_back( ChannelInd1 );
    }
    if ( (ChannelInd2 != -1) && (!Inside( ChannelInd2 , vTempCluster[ClusterID].lChannelInd2 ) ) )
    {
      vTempCluster[ClusterID].EInd2 += vEInd2PointByEvent[index];
      vTempCluster[ClusterID].NInd2 += 1;
      (vTempCluster[ClusterID].lChannelInd2).push_back( ChannelInd2 );
    }

  }
  if ( LogLevel > 3) std::cout << "temporary clusters done" << std::endl;

  std::vector<Cluster> vCluster(n_cluster);

  for( int j = 0 ; j < n_cluster ; j++ )
  {

    if (vTempCluster[j].ECol) 
    {
      vCluster[j].z = ( vTempCluster[j].Sumz )/(vTempCluster[j].ECol);
      vCluster[j].y = ( vTempCluster[j].Sumy )/(vTempCluster[j].ECol);
    }
    else 
    {
      vCluster[j].z = -999.;
      vCluster[j].y = -999.;
      if ( LogLevel > 3) std::cout << "cluster with null energy on colection" << std::endl;
    }

    vCluster[j].Npoint = vTempCluster[j].Npoint;
    vCluster[j].NCol   = vTempCluster[j].NCol;
    vCluster[j].NInd1  = vTempCluster[j].NInd1;
    vCluster[j].NInd2  = vTempCluster[j].NInd2;

    std::vector<int> vtemp_NOF = vTempCluster[j].vNOF;
    if ( AllSame( vtemp_NOF ) )  vCluster[j].NOF = vTempCluster[j].vNOF[0];
    else 
    {
      vCluster[j].NOF = -999;
      if (LogLevel > 2) std::cout << "CLUSTER WITH HITS FROM MULTIPLES TPC" << std::endl;
    }
	  
    vCluster[j].ECol     = vTempCluster[j].ECol;
    vCluster[j].EInd1    = vTempCluster[j].EInd1;
    vCluster[j].EInd2    = vTempCluster[j].EInd2;
    if (vTempCluster[j].NCol) vCluster[j].PeakTime = vTempCluster[j].PeakTime/vTempCluster[j].NCol;
    else
    {
      vCluster[j].PeakTime = -999.;
      if ( LogLevel > 3) std::cout << "cluster with no hits on colection" << std::endl;
    }

    vCluster[j].vMCPDG = vTempCluster[j].vMCPDG;
    vCluster[j].vMCMOMpdg = vTempCluster[j].vMCMOMpdg;
    vCluster[j].vMCWEI = vTempCluster[j].vMCWEI;
    vCluster[j].vMCGenTag = vTempCluster[j].vMCGenTag;

    vCluster[j].vMCX = vTempCluster[j].vMCX;
    vCluster[j].vMCY = vTempCluster[j].vMCY;
    vCluster[j].vMCZ = vTempCluster[j].vMCZ;
  }

  if ( LogLevel > 2) std::cout << "WE HAVE " << (float) out/n_point << " POINTS OUTSIDE OF CLUSTERS" << std::endl;
  return vCluster;
}



std::vector<std::string> pdvdana::SingleHit::GetGeneratorTag( art::Event const &e , std::string fG4Label , int LogLevel , art::ServiceHandle<cheat::BackTrackerService> bt_serv )
{
    //int MCPartcounter = 0;
    std::vector<std::pair<int, std::string>> vTrackIdToLabelPair;

    // get all MC truth object
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mcTruths;
    mcTruths = e.getMany<std::vector<simb::MCTruth>>();

    for (auto const &mcTruth : mcTruths)
    {
      // get generator tag (ie name)
      const std::string &sModuleLabel = mcTruth.provenance()->moduleLabel();

      // MC truth (gen) Mc particle association (g4)
      art::FindManyP<simb::MCParticle> MCTruthsToMCParticles( mcTruth , e , fG4Label );
      std::vector<art::Ptr<simb::MCParticle>> mcParts = MCTruthsToMCParticles.at(0);

      //MCPartcounter += (int) mcParts.size();

      for (const art::Ptr<simb::MCParticle> ptr : mcParts)
      {
        int track_id = ptr->TrackId();
        //creation trackID -> Label association
        vTrackIdToLabelPair.push_back(std::make_pair(track_id, sModuleLabel));
      }

      if( LogLevel > 3) std::cout << "THERE ARE " << (int) mcParts.size() << " MCPARTICLES FROM GENERATOR " << sModuleLabel << std::endl;

    }// end for MCtruth

    //sort to but greates trackID in first position
    std::sort(vTrackIdToLabelPair.begin(), vTrackIdToLabelPair.end(), [](std::pair<int,std::string> a, std::pair<int,std::string> b){ return a.first > b.first;});

    // reassociation for quick access
    std::string noTrackID = "no association";
    std::vector<std::string> vGeneratorLabels( vTrackIdToLabelPair[0].first +1 , noTrackID);

    for(int j = 0 ; j < (int) vTrackIdToLabelPair.size() ; j++)
    {
      if (vGeneratorLabels[vTrackIdToLabelPair[j].first] == noTrackID )
      {
        vGeneratorLabels[vTrackIdToLabelPair[j].first] = vTrackIdToLabelPair[j].second;
      }
      else
      {
        if ( LogLevel > 0) std::cout << "EXCEPTION ERROR : ISSUE WITH ASSOCIATION " << vTrackIdToLabelPair[j].first << std::endl;
        vGeneratorLabels[vTrackIdToLabelPair[j].first] = vTrackIdToLabelPair[j].second;
      }

    }// end for pair(trackID,tag)

    return vGeneratorLabels;
}


DEFINE_ART_MODULE(pdvdana::SingleHit)
