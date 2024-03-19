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
  float fTimeInd1ToInd2 = 0.0;   // in ticktime
  float ftick_in_mus = 0.4;

private:

  // Declare member data here.
  TTree *fAnaTree;

  // Tree variables
  unsigned int fEventID = 0;

  int fHitNumber   = -999;
  int fPlane       = -999;
  int fChannel     = -999;
  int fHitWidth    = -999;
  int fCutHit      = -999;
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

  std::list<int>   lChannelInd1;
  std::list<float> lEnergyInd1;
  std::list<float> lPeakTimeInd1;

  std::list<int>   lChannelInd2;
  std::list<float> lEnergyInd2;
  std::list<float> lPeakTimeInd2;

  //Input variables
  std::string fSpacePointLabel;
  std::string fClusterLabel;
  std::string fTrackLabel;
  std::string fHitLabel;

  int fChannelWd;

  float fPeakTimeWd;
  float fCoincidenceWd;

  // working variables
  int fHitCounter = 0;
  int fNHits      = 0;

  std::list<int> lSingleIndex;
  std::list<int> lIsolatedIndex;

  //function needed
  bool Inside( int k , std::list<int> liste);

  void GetSingle(art::Event const & ev, std::string HitLabel, std::list<int> & index_list_single);

  void GetIsolated(art::Event const & ev, std::string HitLabel, int const ChannelWd, float const PeakTimeWd, std::list<int> & index_list_single, std::list<int> & index_list_isolated);

  void GetTimeCoincidence(art::Event const & ev, std::string HitLabel, const float CoincidenceWd, float const TimeInd1ToInd2, float const ADCtoEl, float const tick_in_mus, const recob::Hit & HitCol,
                                                                                  std::list<int> & ChInd1,
                                                                                  std::list<int> & ChInd2,
                                                                                  std::list<float> & EInd1,
                                                                                  std::list<float> & EInd2,
                                                                                  std::list<float> & PTInd1,
                                                                                  std::list<float> & PTInd2);
};


pdvdana::SingleHit::SingleHit(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSpacePointLabel(p.get<std::string>("SpacePointLabel")),
    fClusterLabel(p.get<std::string>("ClusterLabel")),
    fTrackLabel(p.get<std::string>("TrackLabel")),
    fHitLabel(p.get<std::string>("HitLabel")),

    fChannelWd(p.get<int>("ChannelWindow")),
    fPeakTimeWd(p.get<float>("PeakTimeWindow")), //in ticktime
    fCoincidenceWd(p.get<float>("CoincidenceWindow")) //in ticktime,
    
    // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
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
  
  if( !lSingleIndex.empty()   ) lSingleIndex.clear();
  if( !lIsolatedIndex.empty() ) lIsolatedIndex.clear();

  // Single Isolated concideration
  GetSingle(   e, fHitLabel, lSingleIndex );
  GetIsolated( e, fHitLabel, fChannelWd, fPeakTimeWd, lSingleIndex, lIsolatedIndex);

  for(int index =0 ; index<fNHits; index++)
  {
    const recob::Hit& hit = HitList->at(index);
    fChannel              = hit.Channel();
    fPlane                = hit.WireID().Plane;
    fHitWidth             = hit.EndTick() - hit.StartTick();

    fEnergy         = hit.SummedADC()/fADCtoEl;
    fPeakTime       = hit.PeakTime()*ftick_in_mus;
    fSigmaPeakTime  = hit.SigmaPeakTime()*ftick_in_mus;
    fRMS            = hit.RMS();
    fAmplitude      = hit.PeakAmplitude();
    fSigmaAmplitude = hit.SigmaPeakAmplitude();
    fGoodnessOfFit  = hit.GoodnessOfFit();
    fIntegral       = hit.Integral()/fADCtoEl;
    fSigmaIntegral  = hit.SigmaIntegral();

    if( !lChannelInd1.empty()  ) lChannelInd1.clear();
    if( !lChannelInd2.empty()  ) lChannelInd2.clear();
    if( !lEnergyInd1.empty()   ) lEnergyInd1.clear();
    if( !lEnergyInd2.empty()   ) lEnergyInd2.clear();
    if( !lPeakTimeInd1.empty() ) lPeakTimeInd1.clear();
    if( !lPeakTimeInd2.empty() ) lPeakTimeInd2.clear();

    if (fPlane == 2)
    {
      fHitNumber = fHitCounter;
      ++fHitCounter;

      // Coincidence research
      GetTimeCoincidence( e, fHitLabel, fCoincidenceWd, fTimeInd1ToInd2, fADCtoEl, ftick_in_mus, hit, lChannelInd1, lChannelInd2, lEnergyInd1, lEnergyInd2, lPeakTimeInd1, lPeakTimeInd2);
      fCoincidence = 0;
      std::cout<<"channel in coinc ind1 [ ";
      for (auto ch : lChannelInd1)
      {
        std::cout<< ch <<" ";
      }
      std::cout<<" ]"<<std::endl;
      if ( !lChannelInd1.empty() || !lChannelInd2.empty() ) fCoincidence += 1;
      if ( !lChannelInd1.empty() && !lChannelInd2.empty() ) fCoincidence += 1;
      std::cout<< "col hit num. "<< fHitNumber << " with "<< fCoincidence << " coincidences"<< std::endl;
    }

    if (fPlane != 2)
    {
      fHitNumber   = -1;
      fCoincidence = -999;

      lChannelInd1.clear();
      lChannelInd1.push_back(-1);
      lChannelInd2.clear();
      lChannelInd2.push_back(-1);
      lEnergyInd1.clear();
      lEnergyInd1.push_back(-1);
      lEnergyInd2.clear();
      lEnergyInd2.push_back(-1);
      lPeakTimeInd1.clear();
      lPeakTimeInd1.push_back(-1);
      lPeakTimeInd2.clear();
      lPeakTimeInd2.push_back(-1);
    }
  
    fCutHit = 0;
    if ( Inside(index, lSingleIndex)   ) fCutHit+=1;
    if ( Inside(index, lIsolatedIndex) ) fCutHit+=1;
  
    //Filling tree
    fAnaTree->Fill();

  }
  // setting variables values to initiale values
  fChannel     = -999;
  fPlane       = -999;
  fHitWidth    = -999;
  fHitNumber   = -999;
  fCoincidence = -999;
  fCutHit      = -999;

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

  lChannelInd1.clear();
  lChannelInd2.clear();
  lEnergyInd1.clear();
  lEnergyInd2.clear();
  lPeakTimeInd1.clear();
  lPeakTimeInd2.clear();

  // Implementation of required member function here.
}

void pdvdana::SingleHit::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fAnaTree = tfs->make<TTree>("anaTree", "Output Tree");

  //add branches
  fAnaTree->Branch("eventID"    , &fEventID    );
  fAnaTree->Branch("hitNumber"  , &fHitNumber  );
  fAnaTree->Branch("plane"      , &fPlane      );
  fAnaTree->Branch("channel"    , &fChannel    );
  fAnaTree->Branch("cutHit"     , &fCutHit     );
  fAnaTree->Branch("coincidence", &fCoincidence);

  fAnaTree->Branch("energy"           , &fEnergy        );
  fAnaTree->Branch("peakTime"         , &fPeakTime      );
  fAnaTree->Branch("hitWidth"         , &fHitWidth      );
  fAnaTree->Branch("sigmaPeakTime"    , &fSigmaPeakTime );
  fAnaTree->Branch("amplitudePeaktime", &fAmplitude     );
  fAnaTree->Branch("sigmaAmplitude"   , &fSigmaAmplitude);
  fAnaTree->Branch("rms"              , &fRMS           );
  fAnaTree->Branch("goodnessOfFit"    , &fGoodnessOfFit );
  fAnaTree->Branch("integral"         , &fIntegral      );
  fAnaTree->Branch("sigmaIntegral"    , &fSigmaIntegral );

  fAnaTree->Branch("listChannelInd1"      , &lChannelInd1   );
  fAnaTree->Branch("listEnergyInd1"       , &lEnergyInd1    );
  fAnaTree->Branch("listPeakTimeInd1"     , &lPeakTimeInd1  );
  fAnaTree->Branch("listChannelInd2"      , &lChannelInd2   );
  fAnaTree->Branch("listEnergyInd2"       , &lEnergyInd2    );
  fAnaTree->Branch("listPeakTimeInd2"     , &lPeakTimeInd2  );

}

void pdvdana::SingleHit::endJob()
{
  // Implementation of optional member function here.
}

bool pdvdana::SingleHit::Inside( int k , std::list<int> list){
  return (std::find(list.begin(), list.end(), k) != list.end());
}

void pdvdana::SingleHit::GetSingle(art::Event const & ev, std::string HitLabel, std::list<int> & index_list_single)
{
  auto const hitlist = ev.getValidHandle<vector<recob::Hit>>(HitLabel);
  recob::Hit hit = hitlist->at(0);

  for (int i=0, sz=hitlist->size(); i!=sz; ++i)
  {
    hit = hitlist->at(i);
    if(hit.Multiplicity()!=1) continue;
    index_list_single.push_front(i);    
  } 
}

void pdvdana::SingleHit::GetIsolated(art::Event const & ev, std::string HitLabel, int const ChannelWd, float const PeakTimeWd, std::list<int> & index_list_single, std::list<int> & index_list_isolated)
{
  auto const hitlist = ev.getValidHandle<vector<recob::Hit>>(HitLabel);
  index_list_isolated = index_list_single;

  recob::Hit hit = hitlist->at(0);
  float PeakTime = -999;
  int   Channel  = -999;

  int   ChannelSingle   = -999;
  float PeakTimeSingle  = -999;
  int   ChannelMin      = -999;
  float PeakTimeMin     = -999;
  int   ChannelMax      = -999;
  float PeakTimeMax     = -999;
   
  if(not( index_list_isolated.empty()))
  {
    for (int i=0, sz=hitlist->size(); i!=sz; ++i)
    {
      hit      = hitlist->at(i);
      PeakTime = hit.PeakTime();
      Channel  = hit.Channel();

      ChannelSingle   = -999;
      PeakTimeSingle  = -999;
      ChannelMin      = -999;
      PeakTimeMin     = -999;
      ChannelMax      = -999;
      PeakTimeMax     = -999;

      std::list<int>::iterator elem = index_list_isolated.begin();

      while( elem !=index_list_isolated.end())
      {
        ChannelSingle   = (hitlist->at(*elem)).Channel();
        PeakTimeSingle  = (hitlist->at(*elem)).PeakTime();
      
        ChannelMin      = ChannelSingle  - ChannelWd;
        PeakTimeMin     = PeakTimeSingle - PeakTimeWd;
        ChannelMax      = ChannelSingle  + ChannelWd;
        PeakTimeMax     = PeakTimeSingle + PeakTimeWd;
        
        if ((Channel >= ChannelMin)&&(Channel <= ChannelMax)&&(PeakTime >= PeakTimeMin)&&(PeakTime <= PeakTimeMax))
        {
          if (i != *elem)
          {
            elem = index_list_isolated.erase(elem);
          }
        }
        ++elem;
        
        ChannelSingle   = -999;
        PeakTimeSingle  = -999;
        ChannelMin      = -999;
        PeakTimeMin     = -999;
        ChannelMax      = -999;
        PeakTimeMax     = -999;
      }
      
      PeakTime        = -999;
      Channel         = -999;
    }
  }
}

void pdvdana::SingleHit::GetTimeCoincidence(art::Event const & ev, std::string HitLabel, float const CoincidenceWd, float const TimeInd1ToInd2, float const ADCtoEl, float const tick_in_mus, const recob::Hit & HitCol, 
                                                                                  std::list<int> & ChInd1,
                                                                                  std::list<int> & ChInd2,
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
      //std::cout<<"Ind 1 : ch= "<< hit.Channel() << " peaktime = "<< PeakTime <<std::endl;
      ChInd1.push_back(hit.Channel());
      EInd1.push_back(hit.SummedADC()/ADCtoEl);
      PTInd1.push_back(PeakTime*ftick_in_mus);
      continue;
    }
    if (Plane == 1)
    {
      if ((PeakTime < StartTime2)||(PeakTime > EndTime)) continue;
      //std::cout<<"Ind 2 : ch= "<< hit.Channel() << " peaktime = "<< PeakTime <<std::endl;
      ChInd2.push_back(hit.Channel());
      EInd2.push_back(hit.SummedADC()/fADCtoEl);
      PTInd2.push_back(PeakTime*tick_in_mus);
    }
  }
}


 

DEFINE_ART_MODULE(pdvdana::SingleHit)
