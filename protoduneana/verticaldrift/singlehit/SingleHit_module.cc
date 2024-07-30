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
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

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

using std::vector;
using std::string;

typedef struct { float z, y, ECol, EInd1 ,EInd2 , PeakTime; int group, ch_Col, ch_Ind1, ch_Ind2; } point_t, *point;

typedef struct                                                                
{                                                                                 
  float Sumz , Sumy , SumECol , ECol , EInd1 , EInd2 , PeakTime;                  
  int Npoint, NCol , NInd1 , NInd2;
  std::list<int> lChannelCol , lChannelInd1 , lChannelInd2;
} TempCluster ;                            
                              
typedef struct
{ 
  float z , y , ECol , EInd1 , EInd2 , PeakTime;
  int Npoint, NCol , NInd1 , NInd2;
} Cluster ; 

namespace pdvdana {
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

  int fHitNumber   = -999;
  int fPlane       = -999;
  int fChannel     = -999;
  int fHitWidth    = -999;
  int fTimeIsolation      = -999;
  int fCoincidence = -999;

  float fEnergy         = -999.;
  float fPeakTime       = -999.;
  float fSigmaPeakTime  = -999.;
  float fRMS            = -999.;
  float fAmplitude      = -999.;
  float fSigmaAmplitude = -999.;
  float fGoodnessOfFit  = -999.;
  float fIntegral       = -999.;
  float fSigmaIntegral  = -999.;

  std::list<geo::WireID>   lWireInd1; //working variable
  std::list<int>           lChannelInd1;
  std::list<float>         lEnergyInd1;
  std::list<float>         lPeakTimeInd1;
  std::list<float>         lYInd1;
  std::list<float>         lZInd1;
  std::list<int>           lChIntersectInd1;

  std::list<geo::WireID>   lWireInd2; //working variable
  std::list<int>           lChannelInd2;
  std::list<float>         lEnergyInd2;
  std::list<float>         lPeakTimeInd2;
  std::list<float>         lYInd2;
  std::list<float>         lZInd2;
  std::list<int>           lChIntersectInd2;

  std::list<float> lYPoint;
  std::list<float> lZPoint;
  std::list<float> lEInd1Point;
  std::list<float> lEInd2Point;
  std::list<int>   lChInd1Point;
  std::list<int>   lChInd2Point;

  // Cluster Tree Variables
  int NCluster;

  std::vector<float> vZCluster;
  std::vector<float> vYCluster;
  std::vector<float> vEColCluster;
  std::vector<float> vEInd1Cluster;
  std::vector<float> vEInd2Cluster;
  std::vector<float> vPTCluster;

  std::vector<int>   vNPointCluster;
  std::vector<int>   vNColCluster; 
  std::vector<int>   vNInd1Cluster; 
  std::vector<int>   vNInd2Cluster;

  //Input variables
  std::string fSpacePointLabel;
  std::string fClusterLabel;
  std::string fTrackLabel;
  std::string fHitLabel;

  int fMultiplicity;
  float fTickTimeInMus;

  float fRadiusInt;
  float fRadiusExt;
  float fElectronVelocity;
  float fCoincidenceWd;
  float fTimePlane1ToPlane2; 
  float fPitch;
  float fPitchMultiplier;
  //int   fTagHDVD;  // commented out to make clang happy

  float fNumberInitClusters;
  float fMaxSizeCluster;
  float fMinSizeCluster;
  float fClusterSizeMulti;
  float fNumberConvStep;
  float fCovering;

  // geometry
  const geo::Geometry* fGeom;

  // working variables
  int fHitCounter = 0;
  int fNHits      = 0;

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

  //function needed
  bool Inside( int k , std::list<int> liste);

  float GetDist( float x0 , float y0 , float z0 , float x1 , float y1 , float z1 );

  void GetSingle(art::Event const & ev, std::string HitLabel, std::list<int> & index_list_single, int const Multiplicity);

  void GetTimeIsolation(art::Event const & ev, std::string HitLabel, float const PeakTimeWdInt, float const PeakTimeWdExt, std::list<int> & index_list_single, std::list<int> & index_listIsolatedlated);

  void GetListOfTimeCoincidenceHit(art::Event const & ev, std::string HitLabel, const float CoincidenceWd, float const TimeInd1ToInd2, const recob::Hit & HitCol,
                                                                                  std::list<geo::WireID> & WireInd1,
                                                                                  std::list<geo::WireID> & WireInd2,
                                                                                  std::list<int>   & ChannelInd1,
                                                                                  std::list<int>   & ChannelInd2,
                                                                                  std::list<float> & EInd1,
                                                                                  std::list<float> & EInd2,
                                                                                  std::list<float> & PTInd1,
                                                                                  std::list<float> & PTInd2);

  void GetListOfCrossingChannel( geo::WireID & WireCol , std::list<geo::WireID> & WireInd1 , std::list<geo::WireID> & WireInd2 ,
                              std::list<int>  & ChInd1 , std::list<float> & YInd1 , std::list<float> & ZInd1 , std::list<int>  & ChIntersectInd1 ,
                              std::list<int>  & ChInd2 , std::list<float> & YInd2 , std::list<float> & ZInd2 , std::list<int>  & ChIntersectInd2 );

  void GetListOf3ViewsPoint( float pitch , float alpha ,
                             std::list<int> & ChIntersectInd1 , std::list<float> YInd1 , std::list<float> ZInd1 , std::list<float> EInd1,
                             std::list<int> & ChIntersectInd2 , std::list<float> YInd2 , std::list<float> ZInd2 , std::list<float> EInd2 ,
                             std::list<float> & listYSP       , std::list<float> & listZSP ,
                             std::list<float> & listEind1SP   , std::list<float> & listEind2SP ,
                             std::list<int> & listCh1SP       , std::list<int> & listCh2SP );

  int GetXYZIsolatedPoint( std::vector<float> vYPoint      , std::vector<float> vZPoint      ,
                           std::vector<float> vEInd1Point  , std::vector<float> vEInd2Point  ,
                           std::vector<int>   vChInd1Point , std::vector<int>   vChInd2Point ,
                           std::vector<float> vEnergyCol   , std::vector<float> vPeakTimeCol ,
                           std::vector<int>   vChannelCol  ,
                           float fElectronVelocity , float fTickToMus ,
                           float radiusInt , float radiusExt   ,
                           std::vector<float> &vZIsolated      , std::vector<float> &vYIsolated ,
                           std::vector<float> &vEColIsolated   , std::vector<float> &vPTIsolated ,
                           std::vector<float> &vEInd1Isolated  , std::vector<float> &vEInd2Isolated ,
                           std::vector<int>   &vChColIsolated  , std::vector<int>   &vChInd1Isolated ,
                           std::vector<int>   &vChInd2Isolated );

  // CLUSTER FUNCTIONS
  point gen_zy(int size , std::vector<float> vZ , std::vector<float> vY , 
                          std::vector<float> vECol , std::vector<float> vPT, 
                          std::vector<float> vEInd1 , std::vector<float> vEInd2 , 
                          std::vector<int> vChCol , std::vector<int> vChInd1 , std::vector<int> vChInd2);
  
  float dist2(point a, point b);

  float randf(float m);

  int nearest(point pt, point cent, int n_cluster, float *d2);

  int reallocate(point pt, std::vector<std::vector<float>> ClusterPosition , float threshold);

  float GetDist2D(float z0,float y0,float z1,float y1);

  float mean(float z,float y);
 
  void kpp(point pts, int len, point cent, int n_cent);

  std::vector<std::vector<float> > lloyd(point pts, int len, int n_cluster);

  std::vector<std::vector<float> > GetData(int len,point data);

  std::vector<int> CheckCompletude(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult );

  std::vector<int> CheckClusters(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult , float tmp);

  std::vector<Cluster> GetCluster( int n_point , int n_cluster , point p );
};


pdvdana::SingleHit::SingleHit(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSpacePointLabel(p.get<std::string>("SpacePointLabel")),
    fClusterLabel(p.get<std::string>("ClusterLabel")),
    fTrackLabel(p.get<std::string>("TrackLabel")),
    fHitLabel(p.get<std::string>("HitLabel")),

    fMultiplicity(p.get<int>("HitMultiplicity")),
   
    fRadiusInt(p.get<float>("RadiusInt")),
    fRadiusExt(p.get<float>("RadiusExt")),

    fCoincidenceWd(p.get<float>("CoincidenceWindow")), //in mus,
    fTimePlane1ToPlane2(p.get<float>("TimePlane1ToPlane2")), //in mus,

    fPitch(p.get<float>("Pitch")),
    fPitchMultiplier(p.get<float>("PitchMultiplier")),    

    fNumberInitClusters(p.get<int>("NumberInitClusters")),
    fMaxSizeCluster(p.get<float>("MaxSizeCluster")),
    fMinSizeCluster(p.get<float>("MinSizeCluster")),
    fClusterSizeMulti(p.get<float>("ClusterSizeMulti")),
    fNumberConvStep(p.get<int>("NumberConvStep")),
    fCovering(p.get<float>("Covering"))

    //fTagHDVD(p.get<int>("tagPD"))
    // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
  
  fTickTimeInMus      = sampling_rate(clockData);
  fCoincidenceWd      = fCoincidenceWd/fTickTimeInMus;      // in tt
  fTimePlane1ToPlane2 = fTimePlane1ToPlane2/fTickTimeInMus; // in tt

  fElectronVelocity   = detProp.DriftVelocity();
}

void pdvdana::SingleHit::analyze(art::Event const& e)
{
  //Set event ID
  fEventID = e.id().event();

  // Initializing colection Hit counter
  fHitCounter = 0;

  //Retreive hit list
  //art::InputTag hittag(fHitLabel);
  auto const HitList = e.getValidHandle<vector<recob::Hit>>(fHitLabel);
  fNHits = HitList->size();
 
  if (fNHits == 0)
  {
    fChannel        = -999;
    fPlane          = -999;
    fHitWidth       = -999;

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

    tHitTree->Fill();
    return;
  } 

  if( !lSingleIndex.empty()   ) lSingleIndex.clear();
  if( !lIsolatedIndex.empty() ) lIsolatedIndex.clear();

  // Single Isolated concideration
  GetSingle(   e, fHitLabel, lSingleIndex, fMultiplicity );

  float fPeakTimeWdInt = ( fRadiusInt / fElectronVelocity )/fTickTimeInMus;
  float fPeakTimeWdExt = ( fRadiusExt / fElectronVelocity )/fTickTimeInMus;

  GetTimeIsolation( e, fHitLabel, fPeakTimeWdInt, fPeakTimeWdExt, lSingleIndex, lIsolatedIndex);

  for(int index =0 ; index<fNHits; index++)
  {
    const recob::Hit& hit = HitList->at(index);
    fWire                 = hit.WireID();
    fChannel              = hit.Channel();
    fPlane                = hit.WireID().Plane;
    fHitWidth             = hit.EndTick() - hit.StartTick();

    fEnergy         = hit.SummedADC();///fADCtoEl;
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
      //GetListOfTimeCoincidenceHit( e, fHitLabel, fCoincidenceWd, fTimePlane1ToPlane2 , hit, lWireInd1, lWireInd2, lChannelInd1, lChannelInd2, lEnergyInd1, lEnergyInd2, lPeakTimeInd1, lPeakTimeInd2);
      GetListOfTimeCoincidenceHit( e, fHitLabel, 10 , 0 , hit, lWireInd1, lWireInd2, lChannelInd1, lChannelInd2, lEnergyInd1, lEnergyInd2, lPeakTimeInd1, lPeakTimeInd2);
      fCoincidence = 0;
      if ( !lWireInd1.empty() || !lWireInd2.empty() ) fCoincidence += 1;
      if ( !lWireInd1.empty() && !lWireInd2.empty() ) fCoincidence += 1;

      if ( fCoincidence > 0 )
      {
        GetListOfCrossingChannel( fWire , lWireInd1 , lWireInd2 , lChannelInd1 , lYInd1 , lZInd1 , lChIntersectInd1 , lChannelInd2 , lYInd2 , lZInd2 , lChIntersectInd2 ); 
        GetListOf3ViewsPoint( fPitch , fPitchMultiplier , lChIntersectInd1 , lYInd1 , lZInd1 , lEnergyInd1 , lChIntersectInd2 , lYInd2 , lZInd2 , lEnergyInd2 , lYPoint , lZPoint , lEInd1Point , lEInd2Point , lChInd1Point , lChInd2Point);
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
    }
  
    fTimeIsolation = 0;
    if ( Inside(index, lSingleIndex)   ) fTimeIsolation+=1;
    if ( Inside(index, lIsolatedIndex) ) fTimeIsolation+=1;
  
    //Filling Hit tree
    tHitTree->Fill();

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

      vEnergyColByEvent.push_back( fEnergy   );
      vPeakTimeColByEvent.push_back(  fPeakTime );
      //vPeakTimeColByEvent.push_back( 0 );
      vChannelColByEvent.push_back( fChannel );

      z++;
      ch1++;
      ch2++;
      e1++;
      e2++;

    }
  }// end hit loop

  std::vector<float> vZIsolated;
  std::vector<float> vYIsolated;
  std::vector<float> vEColIsolated;
  std::vector<float> vPTIsolated;
  std::vector<float> vEInd1Isolated;
  std::vector<float> vEInd2Isolated;
  std::vector<int>   vChColIsolated;
  std::vector<int>   vChInd1Isolated;
  std::vector<int>   vChInd2Isolated;

  int nIso = GetXYZIsolatedPoint( vYPointByEvent , vZPointByEvent , vEInd1PointByEvent , vEInd2PointByEvent , vChInd1PointByEvent , vChInd2PointByEvent , vEnergyColByEvent , vPeakTimeColByEvent , vChannelColByEvent , fElectronVelocity , fTickTimeInMus , fRadiusInt , fRadiusExt , vZIsolated , vYIsolated ,vEColIsolated , vPTIsolated , vEInd1Isolated , vEInd2Isolated , vChColIsolated , vChInd1Isolated , vChInd2Isolated );

  if (nIso == 0)
  {
    std::cout << " THERE IS NO ISOLATED POINT IN EVENT " << fEventID << std::endl;

    NCluster = 0;
    vZCluster.push_back(-999);
    vYCluster.push_back(-999);
    vEColCluster.push_back(-999);
    vEInd1Cluster.push_back(-999);
    vEInd2Cluster.push_back(-999);
    vPTCluster.push_back(-999);

    vNPointCluster.push_back(-999);
    vNColCluster.push_back(-999);
    vNInd1Cluster.push_back(-999);
    vNInd2Cluster.push_back(-999);

    tClusterTree->Fill();
    return;
  }

  std::cout << " THERE IS " << vYPointByEvent.size() << " POINTS IN EVENT " << fEventID << std::endl;
  int PTSIsolated = vZIsolated.size();
  std::cout << " THERE IS " << PTSIsolated << " ISOLATED POINTS IN EVENT " << fEventID << std::endl;

  point v = gen_zy( PTSIsolated , vZIsolated , vYIsolated , vEColIsolated , vPTIsolated , vEInd1Isolated , vEInd2Isolated ,  vChColIsolated , vChInd1Isolated , vChInd2Isolated);

  std::vector<std::vector<float> > dataPos = GetData(PTSIsolated,v);
  std::vector<std::vector<float> > clustersPos;

  std::vector<int> vchecks(2,0);

  int K = fNumberInitClusters;

  int check = 0;
  float threshold = fMinSizeCluster;

  for( int i = 0 ; i < fNumberConvStep ; i++ )
  {
    printf("%d %d %d ",i,K,check);

    clustersPos = lloyd(v, PTSIsolated, K);
    vchecks = CheckClusters( dataPos , clustersPos , threshold , fClusterSizeMulti, fCovering);
    check = vchecks[0];

    if(check == 1)  break;
    if(check == 2)
    {
      if (threshold < fMaxSizeCluster)
      {
        threshold *= fClusterSizeMulti;
        printf("Threshold increased \n");
        K = vchecks[1];
      }
      else
      {
        K = vchecks[1] + 1;
        printf("Threshold Max reached \n");
      }
    }
    else K = vchecks[1] + 5;
  }

  printf("Data size : %lu x %lu \n",dataPos.size(),dataPos[0].size());
  printf("Cluster size : %lu x %lu \n",clustersPos.size(),clustersPos[0].size());

  printf("Data clustering ended successfully \n");

  int j = 0;
  point p;
  for (j = 0, p = v; j < PTSIsolated ; j++, p++)
  {
    p->group = reallocate( p , clustersPos , threshold);
  }

  std::vector<Cluster> vCluster = GetCluster( PTSIsolated , clustersPos[0].size() , v );
  NCluster = vCluster.size();

  for( int j = 0 ; j < NCluster ; j++ )
  {
    // Filling tree variables
    vZCluster.push_back(vCluster[j].z);
    vYCluster.push_back(vCluster[j].y);
    vEColCluster.push_back(vCluster[j].ECol);
    vEInd1Cluster.push_back(vCluster[j].EInd1);
    vEInd2Cluster.push_back(vCluster[j].EInd2);
    vPTCluster.push_back(vCluster[j].PeakTime);

    vNPointCluster.push_back(vCluster[j].Npoint);
    vNColCluster.push_back(vCluster[j].NCol);
    vNInd1Cluster.push_back(vCluster[j].NInd1);
    vNInd2Cluster.push_back(vCluster[j].NInd2);

  }
  tClusterTree->Fill();
   

  // setting variables values to initiale values
  fChannel     = -999;
  fPlane       = -999;
  fHitWidth    = -999;
  fHitNumber   = -999;
  fCoincidence = -999;
  fTimeIsolation      = -999;

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


  NCluster = -999;

  vZCluster.clear();
  vYCluster.clear();
  vEColCluster.clear();
  vEInd1Cluster.clear();
  vEInd2Cluster.clear();
  vPTCluster.clear();

  vNPointCluster.clear();
  vNColCluster.clear();
  vNInd1Cluster.clear();
  vNInd2Cluster.clear();

  vYPointByEvent.clear();
  vZPointByEvent.clear();

  vEnergyColByEvent.clear();
  vEInd1PointByEvent.clear();
  vEInd2PointByEvent.clear();

  vPeakTimeColByEvent.clear();

  vChannelColByEvent.clear();
  vChInd1PointByEvent.clear();
  vChInd2PointByEvent.clear();

  // Implementation of required member function here.
}

void pdvdana::SingleHit::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  
  // HIT TREE
  tHitTree = tfs->make<TTree>("hitTree", "Output Tree");

  //add branches
  tHitTree->Branch("eventID"       , &fEventID       );
  tHitTree->Branch("wireID"        , &fWire          );
  tHitTree->Branch("hitNumber"     , &fHitNumber     );
  tHitTree->Branch("plane"         , &fPlane         );
  tHitTree->Branch("channel"       , &fChannel       );
  tHitTree->Branch("timeIsolation" , &fTimeIsolation );
  tHitTree->Branch("coincidence"   , &fCoincidence   );

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
  tHitTree->Branch("listYIntersecPointInd1" , &lYInd1           );
  tHitTree->Branch("listZIntersecPointInd1" , &lZInd1           );
  tHitTree->Branch("listChIntersecInd1"     , &lChIntersectInd1 );

  tHitTree->Branch("listChannelInd2"        , &lChannelInd2     );
  tHitTree->Branch("listEnergyInd2"         , &lEnergyInd2      );
  tHitTree->Branch("listPeakTimeInd2"       , &lPeakTimeInd2    );
  tHitTree->Branch("listYIntersecPointInd2" , &lYInd2           );
  tHitTree->Branch("listZIntersecPointInd2" , &lZInd2           );
  tHitTree->Branch("listChIntersecInd2"     , &lChIntersectInd2 );

  tHitTree->Branch("listOfYPoint"           , &lYPoint      );
  tHitTree->Branch("listOfZPoint"           , &lZPoint      );
  tHitTree->Branch("listEInd1Point"         , &lEInd1Point  );
  tHitTree->Branch("listEInd2Point"         , &lEInd2Point  );
  tHitTree->Branch("listChInd1Point"        , &lChInd1Point );
  tHitTree->Branch("listChInd2Point"        , &lChInd2Point );

  // CLUSTER TREE
  tClusterTree = tfs->make<TTree>("ClusterTree","ClusterTree");

  tClusterTree->Branch("eventID"            , &fEventID       );
  tClusterTree->Branch("NumberOfCluster"    , &NCluster       );
  tClusterTree->Branch("Z"                  , &vZCluster      );
  tClusterTree->Branch("Y"                  , &vYCluster      );
  tClusterTree->Branch("EnergyCollection"   , &vEColCluster   );
  tClusterTree->Branch("EnergyPlane0"       , &vEInd1Cluster  );
  tClusterTree->Branch("EnergyPlane1"       , &vEInd2Cluster  );
  tClusterTree->Branch("NumberOfPoint"      , &vNPointCluster );
  tClusterTree->Branch("NumberOfCollection" , &vNColCluster   );
  tClusterTree->Branch("NumberOfPlane0"     , &vNInd1Cluster  );
  tClusterTree->Branch("NumberOfPlane1"     , &vNInd2Cluster  );
  tClusterTree->Branch("PeakTime"           , &vPTCluster     );
}


void pdvdana::SingleHit::endJob()
{
  // Implementation of optional member function here.
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

void pdvdana::SingleHit::GetTimeIsolation(art::Event const & ev, std::string HitLabel, float const PeakTimeWdInt, float const PeakTimeWdExt,  std::list<int> & index_list_single, std::list<int> & index_listIsolatedlated)
{
  auto const hitlist = ev.getValidHandle<vector<recob::Hit>>(HitLabel);
  index_listIsolatedlated = index_list_single;

  recob::Hit hit = hitlist->at(0);
  float PeakTime = -999;

  float PeakTimeSingle  = -999;

  float PeakTimeMinInt     = -999;
  float PeakTimeMaxInt     = -999;
  float PeakTimeMinExt     = -999;
  float PeakTimeMaxExt     = -999;
   
  if(not( index_listIsolatedlated.empty()))
  {
    for (int i=0, sz=hitlist->size(); i!=sz; ++i)
    {
      hit      = hitlist->at(i);
      PeakTime = hit.PeakTime();

      PeakTimeSingle  = -999;
      PeakTimeMinInt     = -999;
      PeakTimeMaxInt     = -999;
      PeakTimeMinExt     = -999;
      PeakTimeMaxExt     = -999;

      std::list<int>::iterator elem = index_listIsolatedlated.begin();

      while( elem != index_listIsolatedlated.end())
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
            elem = index_listIsolatedlated.erase(elem);
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

void pdvdana::SingleHit::GetListOfTimeCoincidenceHit(art::Event const & ev, std::string HitLabel, float const CoincidenceWd, float const TimeInd1ToInd2, const recob::Hit & HitCol, 
                                                                                  std::list<geo::WireID> & WireInd1,
                                                                                  std::list<geo::WireID> & WireInd2,
                                                                                  std::list<int>   & ChannelInd1,
                                                                                  std::list<int>   & ChannelInd2,
                                                                                  std::list<float> & EInd1,
                                                                                  std::list<float> & EInd2,
                                                                                  std::list<float> & PTInd1,
                                                                                  std::list<float> & PTInd2)
{
  auto const hitlist = ev.getValidHandle<vector<recob::Hit>>(HitLabel);

  recob::Hit hit = hitlist->at(0);

  float PeakTimeCol    = HitCol.PeakTime();
  float RMSPeakTimeCol = HitCol.RMS();

  float EndTime    = PeakTimeCol + RMSPeakTimeCol/2;
  float StartTime1 = PeakTimeCol - CoincidenceWd;
  float StartTime2 = PeakTimeCol - CoincidenceWd + TimeInd1ToInd2;

  float PeakTime = -999;
  int   Plane    = -999;
 
  for (int i=0, sz=hitlist->size(); i!=sz; ++i)
  { 
    hit   = hitlist->at(i);
    Plane = hit.WireID().Plane;
    if (Plane == 2) continue;

    PeakTime = hit.PeakTime();
    if (Plane == 0)
    {
      if ((PeakTime < StartTime1)||(PeakTime > EndTime)) continue;

      WireInd1.push_back(hit.WireID());
      ChannelInd1.push_back(hit.Channel());
      EInd1.push_back(hit.SummedADC());
      PTInd1.push_back(PeakTime);
      continue;
    }
    if (Plane == 1)
    {
      if ((PeakTime < StartTime2)||(PeakTime > EndTime)) continue;

      WireInd2.push_back(hit.WireID());
      ChannelInd2.push_back(hit.Channel());
      EInd2.push_back(hit.SummedADC());
      PTInd2.push_back(PeakTime);
    }
  }
}

void pdvdana::SingleHit::GetListOfCrossingChannel( geo::WireID & WireCol , std::list<geo::WireID> & WireInd1 , std::list<geo::WireID> & WireInd2 , std::list<int>  & ChInd1,
                                                                                                                                                std::list<float> & YInd1 , 
                                                                                                                                                std::list<float> & ZInd1 ,
                                                                                                                                                std::list<int>  & ChIntersectInd1 , 
                                                                                                                                                std::list<int>  & ChInd2 ,
                                                                                                                                                std::list<float> & YInd2 , 
                                                                                                                                                std::list<float> & ZInd2 , 
                                                                                                                                                std::list<int>  & ChIntersectInd2 )
{
  geo::Point_t point = geo::Point_t(-999,-999,-999);
  bool drap;

  std::list<int>::iterator ch1  = ChInd1.begin();
  for (auto const elementInd1 : WireInd1)
  {
    if (WireCol.TPC != elementInd1.TPC )
    {
      ++ch1;
      continue ;
    }
    drap = fGeom->WireIDsIntersect( WireCol , elementInd1 , point);
    if ( drap )
    {
      YInd1.push_back(point.Y());
      ZInd1.push_back(point.Z());
      ChIntersectInd1.push_back(*ch1);
    }
    ++ch1;
  }
  std::list<int>::iterator ch2  = ChInd2.begin();
  for (auto const elementInd2 : WireInd2)
  { 
    if (WireCol.TPC != elementInd2.TPC )
    { 
      ++ch2; 
      continue ;
    }
    drap = fGeom->WireIDsIntersect( WireCol , elementInd2 , point);
    if ( drap ) 
    { 
      YInd2.push_back(point.Y());
      ZInd2.push_back(point.Z());
      ChIntersectInd2.push_back(*ch2);
    }
    ++ch2;
  }
}

void pdvdana::SingleHit::GetListOf3ViewsPoint( float pitch , float alpha , 
                                        std::list<int> & ChIntersectInd1 , std::list<float> YInd1 , std::list<float> ZInd1 , std::list<float> EInd1, 
                                        std::list<int> & ChIntersectInd2 , std::list<float> YInd2 , std::list<float> ZInd2 , std::list<float> EInd2 , 
                                        std::list<float> & listYSP       , std::list<float> & listZSP , 
                                        std::list<float> & listEind1SP   , std::list<float> & listEind2SP , 
                                        std::list<int> & listCh1SP       , std::list<int> & listCh2SP)
{

  std::list<int>::iterator  ch1t = ChIntersectInd1.begin();
  std::list<float>::iterator z1t = ZInd1.begin();
  std::list<float>::iterator e1t = EInd1.begin();

  float dy, dz, dr;

  for( auto const yind1 : YInd1)
  {
    std::list<int>::iterator  ch2t = ChIntersectInd2.begin();
    std::list<float>::iterator z2t = ZInd2.begin();
    std::list<float>::iterator e2t = EInd2.begin();

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


int pdvdana::SingleHit::GetXYZIsolatedPoint( std::vector<float> vYPoint      , std::vector<float> vZPoint      ,
		                            std::vector<float> vEInd1Point  , std::vector<float> vEInd2Point  ,
					    std::vector<int>   vChInd1Point , std::vector<int>   vChInd2Point ,
					    std::vector<float> vEnergyCol   , std::vector<float> vPeakTimeCol , 
					    std::vector<int>   vChannelCol  , 
					    float fElectronVelocity , float fTickToMus ,
					    float radiusInt , float radiusExt   ,
					    std::vector<float> &vZIsolated      , std::vector<float> &vYIsolated ,
                                            std::vector<float> &vEColIsolated   , std::vector<float> &vPTIsolated ,
                                            std::vector<float> &vEInd1Isolated  , std::vector<float> &vEInd2Isolated ,
                                            std::vector<int>   &vChColIsolated  , std::vector<int>   &vChInd1Isolated ,
                                            std::vector<int>   &vChInd2Isolated )
{

  if (vYPoint.size() != vZPoint.size())  throw std::invalid_argument( "BIG PROBLEM" );


  vZIsolated.clear();
  vYIsolated.clear();
  vEColIsolated.clear();
  vPTIsolated.clear();
  vEInd1Isolated.clear();
  vEInd2Isolated.clear();
  vChColIsolated.clear();
  vChInd1Isolated.clear();
  vChInd2Isolated.clear();

  int nIso = 0;
  int npoint = vYPoint.size();
  std::vector<int> vIsIsolated( npoint , -1 );

  for( int k = 0 ; k<npoint ; k++)
  {
    float zIs = vZPoint[k];
    float yIs = vYPoint[k];
    float xIs = vPeakTimeCol[k]/fElectronVelocity/fTickToMus;

    int i = 0;
    bool flag = true;

    int indic = 0;

    if (( yIs == -999) || (zIs == -999))
    {
      vIsIsolated[k] = 0;
      continue;
    }
    if (vIsIsolated[k] == 0)
    {  
      continue;
    }

    while((flag)&&(i<npoint))
    {
      float y = vYPoint[i];
      float z = vZPoint[i];
      float x = vPeakTimeCol[i]/fElectronVelocity/fTickToMus;

      float d = GetDist( xIs , yIs , zIs , x , y , z );

      if (( d > radiusInt)&&( d < radiusInt + radiusExt))
      {
        flag = false;
        indic = i;
      }
      i++;
    }

    if (flag)
    {
      vIsIsolated[k] = 1;
      nIso += 1;

      vZIsolated.push_back(zIs);
      vYIsolated.push_back(yIs);
      vEColIsolated.push_back(vEnergyCol[k]);
      vPTIsolated.push_back(vPeakTimeCol[k]);
      vEInd1Isolated.push_back(vEInd1Point[k]);
      vEInd2Isolated.push_back(vEInd2Point[k]);
      vChColIsolated.push_back(vChannelCol[k]);
      vChInd1Isolated.push_back(vChInd1Point[k]);
      vChInd2Isolated.push_back(vChInd1Point[k]);

      continue;
    }
    else
    {
      vIsIsolated[k] = 0;
      vIsIsolated[indic] = 0;

      continue;
    }
  }

  if ( (int) vZIsolated.size() != nIso ) throw std::invalid_argument( "BIG PROBLEM" );
  return nIso;
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

point pdvdana::SingleHit::gen_zy(int size , std::vector<float> vZ , std::vector<float> vY , std::vector<float> vECol , std::vector<float> vPT, std::vector<float> vEInd1 , std::vector<float> vEInd2 , std::vector<int> vChCol , std::vector<int> vChInd1 , std::vector<int> vChInd2)
{
  int i = 0;
  point p, pt = (point) malloc(sizeof(point_t) * size);

  for (p = pt + size; p-- > pt;)
  {
    p->z = vZ[i];
    p->y = vY[i];
    p->ECol = vECol[i];
    p->EInd1 = vEInd1[i];
    p->EInd2 = vEInd2[i];
    p->ch_Col = vChCol[i];
    p->ch_Ind1 = vChInd1[i];
    p->ch_Ind2 = vChInd2[i];
    p->PeakTime = vPT[i];
    i++;
  }

  return pt;
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


float pdvdana::SingleHit::GetDist2D(float z0,float y0,float z1,float y1){
    float z = z0-z1;
    float y = y0-y1;
    return z*z+y*y;
}

float pdvdana::SingleHit::mean(float z,float y){
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
    std::vector<float> clusterPosZ,clusterPosY;

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

    printf("Counting : %.03f %.03f %.03f sum : %.02f \n",float(Nin)/float(Npts),float(Nin2)/float(Npts),float(Nout)/float(Npts),float(Nin+Nin2+Nout)/float(Npts));


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

    if(overlap>0)
    {
      vector<float> newclusterZ, newclusterY;
      float meanZ, meanY;

      for(int i = 0;i<Ncls;i++)
      {
        meanZ = 0.;
        meanY = 0.;
        overlap = 0;

        for(int j = i;j<Ncls;j++)
        {
          if(IDoverlap[i][j] == 1 && cluster[0][j] != 666)
          {
            meanZ += mean(cluster[0][i],cluster[0][j]);
            meanY += mean(cluster[1][i],cluster[1][j]);
            cluster[0][j] = 666;
            cluster[1][j] = 666;
            overlap++;
          }
        }
        if(overlap == 0 && cluster[0][i] != 666)
        {
          newclusterZ.push_back(cluster[0][i]);
          newclusterY.push_back(cluster[1][i]);
        }
        else if(cluster[0][i] != 666)
        {
          newclusterZ.push_back(meanZ/float(overlap));
          newclusterY.push_back(meanY/float(overlap));
        }
      }

      cluster.clear();
      cluster.push_back(newclusterZ);
      cluster.push_back(newclusterY);

      printf("%lu clusters has been removed \n",Ncls-cluster[0].size());

    }

    NComp = CheckCompletude( data, cluster , RMS , mult );
    Nin  = NComp[0];
    Nin2 = NComp[1];
    Nout = NComp[2];

    printf("Counting : %.03f %.03f %.03f sum : %.02f \n",float(Nin)/float(Npts),float(Nin2)/float(Npts),float(Nout)/float(Npts),float(Nin+Nin2+Nout)/float(Npts));

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

std::vector<Cluster> pdvdana::SingleHit::GetCluster( int n_point , int n_cluster , point p )
{

  int k;
  point vp;
  std::vector<TempCluster> vTempCluster(n_cluster);
  int out = 0;

  for( k = 0 , vp=p ; k<n_point ; k++ , vp++)
  {
    int ClusterID   = vp->group;
    int ChannelCol  = vp->ch_Col;
    int ChannelInd1 = vp->ch_Ind1;
    int ChannelInd2 = vp->ch_Ind2;
    float ECol     = vp->ECol;
    float PeakTime = vp->PeakTime;

    //std::cout << " energy de colection hit " << k << " = " << ECol << std::endl;
    if (ClusterID >=  n_cluster)
    {
      //std::cout << " issue with cluster size " << std::endl;
      //std::cout << ClusterID << std::endl;
      out++;
      continue;
    }

    if (ClusterID ==  -1)
    {
      //std::cout << " issue with cluster size " << std::endl;
      //std::cout << ClusterID << std::endl;
      out++;
      continue;
    }
    vTempCluster[ClusterID].Sumz += ECol * ( vp->z );
    vTempCluster[ClusterID].Sumy += ECol * ( vp->y );
    vTempCluster[ClusterID].SumECol += ECol;
    vTempCluster[ClusterID].Npoint += 1;

    if ( !Inside( ChannelCol , vTempCluster[ClusterID].lChannelCol ) )
    {
      vTempCluster[ClusterID].ECol += ECol;
      vTempCluster[ClusterID].PeakTime += PeakTime;
      vTempCluster[ClusterID].NCol += 1;
      (vTempCluster[ClusterID].lChannelCol).push_back( ChannelCol );
    }
    if ( !Inside( ChannelInd1 , vTempCluster[ClusterID].lChannelInd1 ) )
    {
      vTempCluster[ClusterID].EInd1 += vp->EInd1;
      vTempCluster[ClusterID].NInd1 += 1;
      (vTempCluster[ClusterID].lChannelInd1).push_back( ChannelInd1 );
    }
    if ( !Inside( ChannelInd2 , vTempCluster[ClusterID].lChannelInd2 ) )
    {
      vTempCluster[ClusterID].EInd2 += vp->EInd2;
      vTempCluster[ClusterID].NInd2 += 1;
      (vTempCluster[ClusterID].lChannelInd2).push_back( ChannelInd2 );
    }

  }
  std::cout << "temporary clusters done" << std::endl;

  std::vector<Cluster> vCluster(n_cluster);

  for( int j = 0 ; j < n_cluster ; j++ )
  {

    vCluster[j].z = ( vTempCluster[j].Sumz )/(vTempCluster[j].SumECol);
    vCluster[j].y = ( vTempCluster[j].Sumy )/(vTempCluster[j].SumECol);

    vCluster[j].Npoint = vTempCluster[j].Npoint;
    vCluster[j].NCol   = vTempCluster[j].NCol;
    vCluster[j].NInd1  = vTempCluster[j].NInd1;
    vCluster[j].NInd2  = vTempCluster[j].NInd2;

    vCluster[j].ECol  = vTempCluster[j].ECol;
    vCluster[j].EInd1 = vTempCluster[j].EInd1;
    vCluster[j].EInd2 = vTempCluster[j].EInd2;
    vCluster[j].PeakTime = vTempCluster[j].PeakTime/vTempCluster[j].NCol;

  }

  std::cout << "WE HAVE " << (float) out/n_point << " POINTS OUTSIDE OF CLUSTERS" << std::endl;
  return vCluster;
}












DEFINE_ART_MODULE(pdvdana::SingleHit)
