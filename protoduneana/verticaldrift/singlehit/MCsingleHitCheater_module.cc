////////////////////////////////////////////////////////////////////////
// Class:       MCsingleHitCheater
// Plugin Type: analyzer (Unknown Unknown)
// File:        MCsingleHitCheater_module.cc
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
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

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

namespace pdvdana 
{
  using std::vector;
  using std::string;
  class MCsingleHitCheater;
}


class pdvdana::MCsingleHitCheater : public art::EDAnalyzer {
public:
  explicit MCsingleHitCheater(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCsingleHitCheater(MCsingleHitCheater const&) = delete;
  MCsingleHitCheater(MCsingleHitCheater&&) = delete;
  MCsingleHitCheater& operator=(MCsingleHitCheater const&) = delete;
  MCsingleHitCheater& operator=(MCsingleHitCheater&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  bool Inside( std::string k , std::vector<std::string> v);
  std::vector<std::string> HowManyLabels( std::vector<std::pair<int, std::string>> vTrackIdToLabelPair);

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  unsigned int fEventID = 0;
  const geo::Geometry* fGeom;

  //Input variables
  std::string fHitLabel , fG4Label;
  int LogLevel;
  // working variables
 

};


pdvdana::MCsingleHitCheater::MCsingleHitCheater(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fHitLabel(p.get<std::string>("HitLabel")),
  fG4Label(p.get<std::string>("G4Label")),
  LogLevel(p.get<int>("LogLevel"))

{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
}

void pdvdana::MCsingleHitCheater::analyze(art::Event const& e)
{
  //Set event ID
  fEventID = e.id().event();

  //clock 
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  // back-tracker service
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  art::ValidHandle<std::vector<simb::MCParticle>> MCPartList_handle = e.getValidHandle<vector<simb::MCParticle>>(fG4Label);

  //Retreive MCTruths info list

  std::vector<std::pair<int, std::string>> vTrackIdToLabelPair;

  int MCPartcounter = 0;

  if( !e.isRealData() )
  {

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

      MCPartcounter += (int) mcParts.size();

      for (const art::Ptr<simb::MCParticle> ptr : mcParts) 
      {
        int track_id = ptr->TrackId();
        //creation trackID -> Label association
	vTrackIdToLabelPair.push_back(std::make_pair(track_id, sModuleLabel));
      }
	
      if( LogLevel > 2) std::cout << "THERE ARE " << (int) mcParts.size() << " MCPARTICLES FROM GENERATOR " << sModuleLabel << std::endl;
    
    }// end for MCtruth

    if( LogLevel > 2) 
    { 
      if( MCPartcounter == (int) MCPartList_handle->size()) std::cout << "ALL G4-PARTICLES ARE ASSOCIATED TO A GENERATOR TAG " << std::endl;
      else std::cout << "NOT ALL G4-PARTICLES ARE ASSOCIATED TO A GENERATOR TAG --> ABORT" << std::endl;
    }

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
        std::cout << "ISSUE WITH ASSOCIATION " << vTrackIdToLabelPair[j].first << std::endl;
        vGeneratorLabels[vTrackIdToLabelPair[j].first] = vTrackIdToLabelPair[j].second;
      }

    }// end for pair(trackID,tag)
      

    // get vector of generator tag (ie labels)
    std::vector<std::string> vLabels = HowManyLabels(vTrackIdToLabelPair);

    std::cout << "ASSOCIATION TRACKID GENERATOR TAG MADE ---> HIT ASSOCIATION NOW" <<  std::endl;

    // Hit association

    art::ValidHandle<std::vector<recob::Hit>> HitList_handle = e.getValidHandle<vector<recob::Hit>>(fHitLabel);
    int fNHits = HitList_handle->size();

    if( LogLevel > 2) std::cout << "THERE ARE " << fNHits << " HITS IN EVENT " << fEventID << std::endl;
    std::cout << "STARTING GENRATOR HIT ASSOCIATION " << std::endl;


    std::vector<int> vNumberOfHits((int) vLabels.size(),0);
    int NumberOfElecNoiseHit = 0;

    for(int index = 0 ; index < fNHits ; index++)
    { 

      const recob::Hit& hit = HitList_handle->at(index);
      std::vector<std::pair<int , float>> tempIDPair;

      for (const auto & ide : bt_serv->HitToEveTrackIDEs(clockData, hit))
      {
        //get pair(trckID,energy)
      	tempIDPair.push_back(std::make_pair(ide.trackID,ide.energy));
      }

      // sort by energy to have the main MCparticle for each hit
      std::sort(tempIDPair.begin(), tempIDPair.end(), [](std::pair<int,float> a, std::pair<int,float> b){ return a.second > b.second;});

      if(tempIDPair.size() == 0)
      {
        NumberOfElecNoiseHit += 1;
        continue;
      }
      //get tag associated to TrackID 
      std::string tag = vGeneratorLabels[tempIDPair[0].first];

      for( int k = 0 ; k < (int) vLabels.size() ; k++ )
      {
        if( tag == vLabels[k] ) 
        {
          vNumberOfHits[k] += 1;
          break;
        }
      }
    }// end loop on hits


    if( LogLevel > 2) std::cout << "THERE ARE " << NumberOfElecNoiseHit << " WITHOUT MCTRUTH  (ie Electronic NOISE)" << std::endl;
 
    int sum = 0;
    for( int k = 0 ; k < (int) vLabels.size() ; k++ )
    {  
      if( LogLevel > 2) std::cout << "THERE ARE " << vNumberOfHits[k] << " HITS ASSOCIATED WITH GENERATOR " << vLabels[k] << std::endl;
      sum += vNumberOfHits[k];
    }
    sum += NumberOfElecNoiseHit;
    
    if(sum == fNHits) std::cout << " ALL HITS ARE WELL ASSOCIATED" << std::endl;
      
  }// end if e = simu
}


void pdvdana::MCsingleHitCheater::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
 
}

void pdvdana::MCsingleHitCheater::endJob()
{
  // Implementation of optional member function here.
}
 
bool pdvdana::MCsingleHitCheater::Inside( std::string k , std::vector<std::string> v){
  return (std::find(v.begin(), v.end(), k) != v.end());
}

std::vector<std::string> pdvdana::MCsingleHitCheater::HowManyLabels( std::vector<std::pair<int, std::string>> vTrackIdToLabelPair) 
{
  std::vector<std::string> v;
  for( std::pair<int, std::string> pair : vTrackIdToLabelPair)
  {
    if( !Inside(pair.second , v) )
    {
      v.push_back( pair.second );
    }
  }

  return v;
}

DEFINE_ART_MODULE(pdvdana::MCsingleHitCheater)

