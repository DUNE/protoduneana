////////////////////////////////////////////////////////////////////////
// Class:       EMshower
// Plugin Type: analyzer (Unknown Unknown)
// File:        EMshower_module.cc
//
// Generated at Fri Jan 23 11:56:46 2026 by Yoann Kermaidic using cetskelgen
// from cetlib version 3.18.02.
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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"

#include "TTree.h"

namespace pdvd {
  class EMshower;
}


class pdvd::EMshower : public art::EDAnalyzer {
public:
  explicit EMshower(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EMshower(EMshower const&) = delete;
  EMshower(EMshower&&) = delete;
  EMshower& operator=(EMshower const&) = delete;
  EMshower& operator=(EMshower&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::string fBeamEventLabel;
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fSpacePointLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;

  TTree* fTree;

  unsigned int fEventID;
  float fBeamTof;
  short fBeamCKov0;
  short fBeamCKov1;


  int fNPFParticles;
  int fNPrimaryChildren;
  int fPrimaryPdg;
  int fPrimaryKey;

  std::vector<int> fChildBeamPFPKey;
  std::vector<int>  fChildBeamPFPPdg;
  std::vector<int>  fChildBeamPFPChildren;
  std::vector<int>  fChildBeamTrackFlag;

  std::vector<float> fChildBeamPFPTrackScore;
  std::vector<float> fChildBeamLength;
  std::vector<float> fChildBeamEnergy;
  
  std::vector<std::vector<float>> fChildBeamX;
  std::vector<std::vector<float>> fChildBeamY;
  std::vector<std::vector<float>> fChildBeamZ;

};


pdvd::EMshower::EMshower(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fBeamEventLabel(p.get<std::string>("BeamEventLabel")), 
  fSliceLabel(p.get<std::string>("SliceLabel")), 
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")), 
  fSpacePointLabel(p.get<std::string>("SpacePointLabel")), 
  fTrackLabel(p.get<std::string>("TrackLabel")), 
  fShowerLabel(p.get<std::string>("ShowerLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdvd::EMshower::analyze(art::Event const& e)
{
  fEventID = e.id().event();
  fBeamTof = 0.0;

  // Get the beam instrumentation related objects
  art::ValidHandle<std::vector<beam::ProtoDUNEBeamEvent>> beamHandle = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventLabel);
  std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;

  if(!beamHandle.isValid()) return;

  art::fill_ptr_vector(beaminfo, beamHandle);

  if ( beaminfo.size() == 0 ) { std::cout << "WARNING: Beam event vector is empty." << std::endl; return;}
  if ( beaminfo.size() >  1 ) { std::cout << "WARNING: Beam event vector has size " << beaminfo.size() << std::endl;}

  // for beam particle selection strategy 
  // see p. 18 https://inspirehep.net/files/134bfcc0c74ca45b565e7e9514434c24
  int beamTrigFlag =  beaminfo[0]->GetTimingTrigger();
  std::cout << "Beam event trigger is " << beamTrigFlag << std::endl;

  // beam inst. Time Of Flight
  fBeamTof = beaminfo[0]->GetTOF();
  std::cout << "Beam event TOF is " << fBeamTof << std::endl;

  // beam inst. Cerenkov trigger
  fBeamCKov0 = beaminfo[0]->GetCKov0Status();
  fBeamCKov1 = beaminfo[0]->GetCKov1Status();
  std::cout << "Beam event CKovs are " << fBeamCKov0 << " / " << fBeamCKov1 << std::endl;

  // Get the pandora related objects
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle    = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  // TODO: commented out below to avoid unused variable warning in c14:prof qualifier
  //art::ValidHandle<std::vector<recob::Track>> trackHandle    = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
  //art::ValidHandle<std::vector<recob::Shower>> showerHandle  = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);
  //art::ValidHandle<std::vector<recob::SpacePoint>> spHandle  = e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);

  // Associate PFPs to tracks, showers and space points
  art::FindManyP<recob::Track>  pfpTrackAssoc(                   pfpHandle, e, fTrackLabel);
  art::FindManyP<recob::Shower> pfpShowerAssoc(                  pfpHandle, e, fShowerLabel);
  art::FindManyP<recob::SpacePoint> pfpSPAssoc(                  pfpHandle, e, fSpacePointLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetaAssoc(pfpHandle, e, fPFParticleLabel);

  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if(!sliceHandle.isValid()) return;

  art::fill_ptr_vector(sliceVector, sliceHandle);
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fPFParticleLabel);

  for(const art::Ptr<recob::Slice> &slice: sliceVector){
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));
   
    //std::cout << "slice: " << slice.key() << std::endl;
 
    for(const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs){
      const bool isPrimary(slicePFP->IsPrimary());

      fPrimaryPdg       = 0;
      fPrimaryKey       = 0;
      fNPFParticles     = 0;
      fNPrimaryChildren = 0;

      fChildBeamPFPKey.clear();
      fChildBeamPFPPdg.clear();
      fChildBeamPFPChildren.clear();
     
      fChildBeamPFPTrackScore.clear();
      fChildBeamTrackFlag.clear();
      fChildBeamLength.clear();
      fChildBeamEnergy.clear();

      fChildBeamX.clear();
      fChildBeamY.clear();
      fChildBeamZ.clear();

      // Only consider primary particles to idenfy the incoming beam particle
      if(!isPrimary)                          continue;
      // Skip primaires idendified as muons for now - will better consider the beam trigger information later
      // if(std::abs(slicePFP->PdgCode()) == 13) continue;
	
      std::cout << "  found primary particle at slice " << slice.key() << " with PDG code: " << slicePFP->PdgCode() << std::endl;
     
      fPrimaryKey       = slice.key();
      fPrimaryPdg       = slicePFP->PdgCode(); 
      fNPFParticles     = slicePFPs.size();
      fNPrimaryChildren = slicePFP->NumDaughters();
 
      std::vector<art::Ptr<recob::PFParticle>> beamPFPs(slicePFPAssoc.at(fPrimaryKey));

      // Loop over the PFPs of the primary beam particle
      for(const art::Ptr<recob::PFParticle> &beamPFP : beamPFPs){

	fChildBeamPFPKey.push_back(beamPFP.key());
	fChildBeamPFPPdg.push_back(beamPFP->PdgCode());
	fChildBeamPFPChildren.push_back(beamPFP->NumDaughters());

	// space points distribution
	std::vector<art::Ptr<recob::SpacePoint>> spacepoints = pfpSPAssoc.at(beamPFP.key());
	
	std::vector<float> tmpX, tmpY, tmpZ;
	for (const auto& sp : spacepoints){
	  // remove badly reconstructed space points	  
	  if(sp->XYZ()[2] < 0) continue;

	  tmpX.push_back(sp->XYZ()[0]);
	  tmpY.push_back(sp->XYZ()[1]);
	  tmpZ.push_back(sp->XYZ()[2]);
	} // end of shower space points loop
	
	fChildBeamX.push_back(tmpX);
	fChildBeamY.push_back(tmpY);
	fChildBeamZ.push_back(tmpZ);

	// Retrieve the track score for PFP metadata
	if (pfpMetaAssoc.isValid()){
          if (pfpMetaAssoc.at(beamPFP->Self()).size() > 0){
            const larpandoraobj::PFParticleMetadata metaData = *((pfpMetaAssoc.at(beamPFP->Self())).at(0));
            const std::map<std::string, float> &propMap = metaData.GetPropertiesMap();
            if(propMap.count("TrackScore") > 0) fChildBeamPFPTrackScore.push_back(propMap.at("TrackScore"));
	    else fChildBeamPFPTrackScore.push_back(-1);
          }
	  else fChildBeamPFPTrackScore.push_back(-1);
        }
	else fChildBeamPFPTrackScore.push_back(-1);

	// Retrieve tracks from PFPs list
	std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(beamPFP.key());
	// Retrieve shower from PFPs list
	std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(beamPFP.key());

	if(tracks.size() == 1){
	  art::Ptr<recob::Track> track = tracks.at(0);
	  fChildBeamTrackFlag.push_back(1);
	  fChildBeamLength.push_back(track->Length());
	  fChildBeamEnergy.push_back(-1);
	} // end of track condition - avoid double counting of tracks and shower
	else if(showers.size() == 1){
	  art::Ptr<recob::Shower> shower = showers.at(0);
	  fChildBeamTrackFlag.push_back(0);
	  fChildBeamLength.push_back(shower->Length());
	  fChildBeamEnergy.push_back(shower->Energy()[2]);
	} // end of shower condition
      } // end of beamPFP loop
    
      // Filling only primary particles
      fTree->Fill();

    } // end of PFP loop
  } //end of slice loop
}

void pdvd::EMshower::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree","Pandora PFParticles");

  // event
  fTree->Branch("eventID",               &fEventID);
  // beam instrumentation
  fTree->Branch("BeamTOF",               &fBeamTof);
  fTree->Branch("BeamCKov0",             &fBeamCKov0);
  fTree->Branch("BeamCKov1",             &fBeamCKov1);
  // pandora PFPs - primary
  fTree->Branch("PrimaryKey",            &fPrimaryKey);
  fTree->Branch("PrimaryPdg",            &fPrimaryPdg);
  fTree->Branch("NPFParticles",          &fNPFParticles);
  // pandora PFPs - children 
  fTree->Branch("NPrimaryChildren",      &fNPrimaryChildren);
  fTree->Branch("ChildBeamPFPKey",       &fChildBeamPFPKey);
  fTree->Branch("ChildBeamPFPPdg",       &fChildBeamPFPPdg);
  fTree->Branch("ChildBeamPFPChildren",  &fChildBeamPFPChildren);
  fTree->Branch("ChildBeamPFPTrackScore",&fChildBeamPFPTrackScore);
  fTree->Branch("ChildBeamPFPTrackFlag", &fChildBeamTrackFlag);
  fTree->Branch("ChildBeamPFPLength",    &fChildBeamLength);
  fTree->Branch("ChildBeamPFPEnergy",    &fChildBeamEnergy);
  fTree->Branch("ChildBeamPFPX",         &fChildBeamX);
  fTree->Branch("ChildBeamPFPY",         &fChildBeamY);
  fTree->Branch("ChildBeamPFPZ",         &fChildBeamZ);
}

void pdvd::EMshower::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pdvd::EMshower)
