////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEelectronWireAna
// Plugin Type: analyzer (art v3_03_01)
// File:        ProtoDUNEelectronWireAna_module.cc
//
// Generated at Mon Dec  2 16:17:42 2019 by Aaron Higuera Pichardo using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

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

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "TFile.h"
#include "TProfile.h"
#include "TH1.h"
#include "TTree.h"
#include <stdio.h>
#include <stdlib.h> 
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

const int MAXHits = 6000;

class ProtoDUNEelectronWireAna;


class ProtoDUNEelectronWireAna : public art::EDAnalyzer {
public:
  explicit ProtoDUNEelectronWireAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtoDUNEelectronWireAna(ProtoDUNEelectronWireAna const&) = delete;
  ProtoDUNEelectronWireAna(ProtoDUNEelectronWireAna&&) = delete;
  ProtoDUNEelectronWireAna& operator=(ProtoDUNEelectronWireAna const&) = delete;
  ProtoDUNEelectronWireAna& operator=(ProtoDUNEelectronWireAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;
  void beginJob() override;
private:

  void Initialize();
  // Declare member data here.
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
 
  geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fWireTag;

  TTree *fTree;
  unsigned int frun;
  unsigned int fsubrun;
  unsigned int fevent;
  std::vector<int>  fprimary_Shower_wire_w;
  std::vector<float>  fprimary_Shower_wire_ch; //charge 
  std::vector<float>fprimary_Shower_wire_X;  
  std::vector<float>fprimary_Shower_wire_Z;  
  std::vector<float>fprimary_Shower_wire_Y;  

  std::vector<float>  fprimary_Shower_MCwire_E;  //energy in MeV
  std::vector<int>  fprimary_Shower_MCwire_w;
  float fprimaryStartPosition[3];
};

void ProtoDUNEelectronWireAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("WireAna","Recob::Wire");
  fTree->Branch("run",&frun,"run/I");
  fTree->Branch("subrun",&fsubrun,"subrun/I");
  fTree->Branch("event",&fevent,"event/I");
  
  fTree->Branch("primaryStartPosition",          &fprimaryStartPosition,         "primaryStartPosition[3]/f");
  fTree->Branch("primary_Shower_wire_w",&fprimary_Shower_wire_w);
  fTree->Branch("primary_Shower_wire_ch",&fprimary_Shower_wire_ch);
  fTree->Branch("primary_Shower_wire_X",&fprimary_Shower_wire_X);
  fTree->Branch("primary_Shower_wire_Z",&fprimary_Shower_wire_Z);
  fTree->Branch("primary_Shower_wire_Y",&fprimary_Shower_wire_Y);
  fTree->Branch("primary_Shower_MCwire_w",&fprimary_Shower_MCwire_w);
  fTree->Branch("primary_Shower_MCwire_E",&fprimary_Shower_MCwire_E);
 
}


ProtoDUNEelectronWireAna::ProtoDUNEelectronWireAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fWireTag(p.get<std::string>("WireTag"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ProtoDUNEelectronWireAna::analyze(art::Event const& evt)
{

  Initialize(); 
  frun = evt.run();
  fsubrun = evt.subRun();
  fevent  = evt.id().event(); 
 
  //check for reco pandora stuff
  auto recoParticleHandle = evt.getHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  if( !recoParticleHandle ) return;

  //std::vector<int> hit_w;
  std::map<int,std::vector<double>> hit_w_and_t1;
  std::map<int,std::vector<double>> hit_w_and_t2;
  std::map<int,std::vector<double>> hit_w_and_y;
  std::map<int,std::vector<double>> hit_w_and_x;
  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  const simb::MCParticle* mcparticle = NULL;
  bool doAna = false;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
  for(const recob::PFParticle* particle : pfParticles){
     const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
     if( thisShower == 0x0 ) continue;
     if( !evt.isRealData() ){
       mcparticle = truthUtil.GetMCParticleFromRecoShower(clockData, *thisShower, evt, "pandoraShower");
       if( abs(mcparticle->PdgCode()) != 11 ) return;
     }
     fprimaryStartPosition[0] = thisShower->ShowerStart().X();
     fprimaryStartPosition[1] = thisShower->ShowerStart().Y();
     fprimaryStartPosition[2] = thisShower->ShowerStart().Z();
     doAna = true; 
     const std::vector<const recob::Hit*> sh_hits = showerUtil.GetRecoShowerHits(*thisShower, evt, fShowerTag);
     art::FindManyP<recob::Wire> wFromHits(sh_hits,evt,"hitpdune");
     art::FindManyP<recob::SpacePoint> spFromShowerHits(sh_hits,evt,fPFParticleTag);
     for( size_t j=0; j<sh_hits.size() && j<MAXHits; ++j){
       if( sh_hits[j]->WireID().Plane != 2 ) continue;
       std::vector<art::Ptr<recob::Wire>> wires = wFromHits.at(j);
       //hit_w.push_back(wires[0]->Channel());
       hit_w_and_t1[wires[0]->Channel()].push_back(sh_hits[j]->PeakTimeMinusRMS(5.0));
       hit_w_and_t2[wires[0]->Channel()].push_back(sh_hits[j]->PeakTimePlusRMS(5.0));
       //std::cout<<wires[0]->Channel()<<" "<<sh_hits[j]->PeakTime()<<" "<<sh_hits[j]->PeakTimePlusRMS(5.0)<<" "<<sh_hits[j]->PeakTimeMinusRMS(5.0)<<std::endl;
       hit_w_and_x[wires[0]->Channel()].push_back(detProp.ConvertTicksToX(sh_hits[j]->PeakTime(),sh_hits[j]->WireID().Plane,sh_hits[j]->WireID().TPC,0));
       std::vector<art::Ptr<recob::SpacePoint>> sp = spFromShowerHits.at(j); 
       if(!sp.empty()){
          hit_w_and_y[wires[0]->Channel()].push_back(sp[0]->XYZ()[1]);
       }
       else hit_w_and_y[wires[0]->Channel()].push_back(fprimaryStartPosition[1]); //use vtx if no sp 
     }   
     break;
  }
 
  //one wire can have various hits so lets remove duplicate wires
  //std::sort(hit_w.begin(),hit_w.end());
  //hit_w.erase(std::unique(hit_w.begin(),hit_w.end()), hit_w.end()); 
  
  if( doAna ){
    auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >(fWireTag);
    //auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >("wclsdatasp:gauss"); //new reco 
    auto w1 = hit_w_and_t1.begin();
    auto w2 = hit_w_and_t2.begin();
    auto x  = hit_w_and_x.begin();
    auto y  = hit_w_and_y.begin(); 
    while( w1 != hit_w_and_t1.end()){
      int it_w = w1->first;
      int n_hits = w1->second.size();
      // std::cout<<it_w<<" "<<n_hits<<" "<<w1->second[n_hits-1]<<" "<<w2->second[n_hits-1]<<std::endl;
      double t1 =  w1->second[0]; //first hit
      double t2 =  w2->second[n_hits-1]; //last hit
      double x_w, y_w;
      if( x->second.size() < 1 ){
        x_w = (x->second[0]-x->second[n_hits-1])/2.0;
        y_w = (y->second[0]-y->second[n_hits-1])/2.0;
      }
      else {
        x_w = x->second[0];
        y_w = y->second[0];
      }
      for(auto & wire : * wires){
         int channel_no = wire.Channel();
         int plane = fWireReadoutGeom.View(channel_no);
         if( plane != 2 ) continue;
         std::vector< geo::WireID > wireID= fWireReadoutGeom.ChannelToWire(channel_no);
         const geo::WireGeo* pwire = fWireReadoutGeom.WirePtr(wireID.at(0)); //for collection plane there is one wire per channel
         auto xyzWire = pwire->GetCenter();
         if( it_w == channel_no ){  
           double charge =0.0;
           for(size_t i = 0; i < wire.Signal().size(); ++i){
              if( i > t1 && i < t2 ) charge += wire.Signal()[i];
           }   
           fprimary_Shower_wire_Z.push_back(xyzWire.Z()); 
           fprimary_Shower_wire_Y.push_back(y_w); 
           fprimary_Shower_wire_X.push_back(x_w); 
           fprimary_Shower_wire_ch.push_back(charge);
           fprimary_Shower_wire_w.push_back(channel_no);
           
           break;
         }   
      }// all recob::wire
      w1 ++; w2 ++;
      x ++; y ++;
    }//wire from shower 

    //look at MC info if available
    if(!evt.isRealData()){ 
      auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
      const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
      int trackid = mcparticle->TrackId();
      if( mcparticle->TrackId() == geantGoodParticle->TrackId() ){
        if( doAna ){   //we have a reco shower
          auto simchannelHandle = evt.getHandle< std::vector<sim::SimChannel> >("largeant");
          if (simchannelHandle){
            for(auto const& simchannel : (*simchannelHandle)){
               if(fWireReadoutGeom.View(simchannel.Channel()) != 2) continue;
               auto const& alltimeslices = simchannel.TDCIDEMap();
               double EperCh=0;
               for(auto const& tslice : alltimeslices){
	          auto const& simide = tslice.second;
	          // Loop over energy deposits
	          for(auto const& eDep : simide){
	             if(eDep.trackID == trackid || eDep.trackID == -trackid){
                       EperCh += eDep.energy;
                     }
	          }
               } 
               if( EperCh== 0 ) continue;//save only channeles with ID
               fprimary_Shower_MCwire_w.push_back(simchannel.Channel());
               fprimary_Shower_MCwire_E.push_back(EperCh);
            }
          }
        }//do analysis of MC
      }//is MC same as beam?
    }//is MC?
  }//shower


  fTree->Fill();
}
void ProtoDUNEelectronWireAna::Initialize(){
  frun =-999;
  fsubrun =-999;
  fevent =-999; 
  fprimary_Shower_wire_w.clear();
  fprimary_Shower_wire_Y.clear();
  fprimary_Shower_wire_X.clear();
  fprimary_Shower_wire_Z.clear();
  fprimary_Shower_wire_ch.clear();
  fprimary_Shower_MCwire_w.clear();
  fprimary_Shower_MCwire_E.clear();
  for(int k=0; k < 3; k++){
    fprimaryStartPosition[k] = -999.0;
  }
}
DEFINE_ART_MODULE(ProtoDUNEelectronWireAna)
