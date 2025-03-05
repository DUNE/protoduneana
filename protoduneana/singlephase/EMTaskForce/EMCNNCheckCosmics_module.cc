////////////////////////////////////////////////////////////////////////
// Class:       EMCNNCheckCosmics
// Plugin Type: analyzer (art v3_03_01)
// File:        EMCNNCheckCosmics_module.cc
//
// Generated at Wed Jan  8 21:50:23 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_08_00.
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

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "TTree.h"

#include <iostream>

namespace pdsp {
  class EMCNNCheckCosmics;
}


class pdsp::EMCNNCheckCosmics : public art::EDAnalyzer {
public:
  explicit EMCNNCheckCosmics(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EMCNNCheckCosmics(EMCNNCheckCosmics const&) = delete;
  EMCNNCheckCosmics(EMCNNCheckCosmics&&) = delete;
  EMCNNCheckCosmics& operator=(EMCNNCheckCosmics const&) = delete;
  EMCNNCheckCosmics& operator=(EMCNNCheckCosmics&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;

private:

  void ResetVars();
  
  //int fselectpdg;
  std::string fGeneratorTag;
  std::string fCNNTag;

  TTree *ftree;
  int run;
  int subrun;
  int event;
  int beampdg;
  double average_score_em;
  double average_score_trk;
  double average_score_mic;
  double track_endz;
  int ndaughterhits;
  double average_daughter_score_mic;
  double vtxx, vtxy, vtxz;
  double endx, endy, endz;
  double dirx, diry, dirz;
  std::vector<short> channel;
  std::vector<short> tpc;
  std::vector<short> plane;
  std::vector<short> wire;
  std::vector<double> charge;
  std::vector<double> peakt;
  std::vector<double> score_em;
  std::vector<double> score_trk;
  std::vector<double> score_mic;

  std::vector<short> daughter_channel;
  std::vector<short> daughter_tpc;
  std::vector<short> daughter_plane;
  std::vector<short> daughter_wire;
  std::vector<double> daughter_charge;
  std::vector<double> daughter_peakt;
  std::vector<double> daughter_score_em;
  std::vector<double> daughter_score_trk;
  std::vector<double> daughter_score_mic;

  std::vector<int> pdg;
  std::vector<int> origin;
  std::vector<std::string> process;

};


pdsp::EMCNNCheckCosmics::EMCNNCheckCosmics(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
//fselectpdg(p.get<int>("selectpdg")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fCNNTag(p.get<std::string>("CNNTag","emtrkmichelid:emtrkmichel"))
  {
  }

void pdsp::EMCNNCheckCosmics::ResetVars()
{
  beampdg = 0;
  average_score_em  = 0.;
  average_score_trk = 0.;
  average_score_mic = 0.;
  track_endz = -1;
  ndaughterhits = 0;
  average_daughter_score_mic = 0.;
  vtxx = -9999;
  vtxy = -9999;
  vtxz = -9999;
  endx = -9999;
  endy = -9999;
  endz = -9999;
  dirx = -9999;
  diry = -9999;
  dirz = -9999;
  channel.clear();
  tpc.clear();
  plane.clear();
  wire.clear();
  charge.clear();
  peakt.clear();
  score_em.clear();
  score_trk.clear();
  score_mic.clear();

  daughter_channel.clear();
  daughter_tpc.clear();
  daughter_plane.clear();
  daughter_wire.clear();

  daughter_charge.clear();
  daughter_peakt.clear();
  daughter_score_em.clear();
  daughter_score_trk.clear();
  daughter_score_mic.clear();

  pdg.clear();
  origin.clear();
  process.clear();

}

void pdsp::EMCNNCheckCosmics::analyze(art::Event const& e)
{
  run = e.run();
  subrun = e.subRun();
  event = e.id().event();

  //Services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);
  art::ServiceHandle<geo::Geometry const> geom;

  std::vector < art::Ptr < recob::Slice > > sliceList;
  auto sliceListHandle = e.getHandle < std::vector < recob::Slice > >("pandora");
  if (sliceListHandle) {
    art::fill_ptr_vector(sliceList, sliceListHandle);
  }
  else return;

  // Get all pfparticles
  std::vector < art::Ptr < recob::PFParticle > > pfpList;
  auto pfpListHandle = e.getHandle < std::vector < recob::PFParticle > >("pandora");
  if (pfpListHandle) {
    art::fill_ptr_vector(pfpList, pfpListHandle);
  }

  // Get all clusters
  std::vector < art::Ptr < recob::Cluster > > cluList;
  auto cluListHandle = e.getHandle < std::vector < recob::Cluster > >("pandora");
  if (cluListHandle) {
    art::fill_ptr_vector(cluList, cluListHandle);
  }

  // Get all tracks
  std::vector < art::Ptr < recob::Track > > trkList;
  auto trkListHandle = e.getHandle < std::vector < recob::Track > >("pandoraTrack");
  if (trkListHandle) {
    art::fill_ptr_vector(trkList, trkListHandle);
  }

  // Get all hits
  std::vector < art::Ptr < recob::Hit > > hitList;
  auto hitListHandle = e.getHandle < std::vector < recob::Hit > >("hitpdune");
  if (hitListHandle) {
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  // Get hits-track association
  art::FindManyP<recob::Cluster> fmcpfp(pfpListHandle, e, "pandora");

  // Get vertex-PFParticle association
  art::FindManyP<recob::Vertex> fmvpfp(pfpListHandle, e, "pandora");

  // Get hit-cluster association
  art::FindManyP<recob::Hit> fmhc(cluListHandle, e, "pandora");

  art::FindManyP <recob::Hit> hitsFromSlice(sliceListHandle, e, "pandora");

  // Get track-hit association
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trkListHandle, e,"pandoraTrack"); // to associate tracks and hits

  // Get hit-track association
  art::FindManyP<recob::Hit> thass(trkListHandle, e, "pandoraTrack"); //to associate hit just trying

  anab::MVAReader<recob::Hit,4> hitResults(e, fCNNTag);

  // Get the PFParticle utility
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  const unsigned int beamSlice{pfpUtil.GetBeamSlice(e,"pandora")};

  const std::map<unsigned int,std::vector<const recob::PFParticle*>> sliceToPFParticleMap = pfpUtil.GetPFParticleSliceMap(e,"pandora");
  std::vector<unsigned int> tracksToUse;
  for (auto const & element : sliceToPFParticleMap)
  {
    if (element.first == beamSlice)
      continue;

    for (const recob::PFParticle *part : element.second)
    {
      // We only want track-like particles
      const recob::Track *track{pfpUtil.GetPFParticleTrack(*part,e,"pandora","pandoraTrack")};
      if (track == nullptr)
        continue;
  
      // Check that the track is long enough to be considered a cosmic muon
      if (track->Length() < 100.)
        continue;
  
      // Ignore tracks that start or end near the front face of the TPC
      if (track->Start<TVector3>().Z() < 50.0 || track->End<TVector3>().Z() < 50.0)
        continue;
 
      // Angular cut for very steep tracks
      if (track->VertexDirection<TVector3>().Y() < TMath::Cos(165. * TMath::Pi() / 180.))
        continue;
 
      // Now we should be in the situation where we have a good track
      tracksToUse.push_back(track->ID());
    }
  }

  // We can now look at these particles
  int trackid = -1;
//  int endwire = -1;
//  int endpeakt = -1;
//  int endtpc = -1;
  std::vector<int> wirekeys;
  for(unsigned int t = 0; t < tracksToUse.size(); ++t){

    const art::Ptr<recob::Track> thisTrack = trkList.at(tracksToUse.at(t));
    if (thisTrack){
      this->ResetVars();
      //if (!beam_cuts.IsBeamlike(*thisTrack, e, "1")) return;
      // Track ID
      trackid = thisTrack->ID();
      // Track end point z
      track_endz = thisTrack->End().Z();
      vtxx = thisTrack->Vertex().X();
      vtxy = thisTrack->Vertex().Y();
      vtxz = thisTrack->Vertex().Z();
      if (!geom->FindTPCAtPosition(geo::Point_t(vtxx, vtxy, vtxz)).isValid) return;
      auto offset = SCE->GetCalPosOffsets(geo::Point_t(vtxx, vtxy, vtxz), (geom->FindTPCAtPosition(geo::Point_t(vtxx, vtxy, vtxz))).TPC);
//      std::cout<<"track "<<offset.X()<<" "<<offset.Y()<<" "<<offset.Z()<<std::endl;
      vtxx -= offset.X();
      vtxy += offset.Y();
      vtxz += offset.Z();
      endx = thisTrack->End().X();
      endy = thisTrack->End().Y();
      endz = thisTrack->End().Z();
      if (!geom->FindTPCAtPosition(geo::Point_t(endx, endy, endz)).isValid) return;
      offset = SCE->GetCalPosOffsets(geo::Point_t(endx, endy, endz), (geom->FindTPCAtPosition(geo::Point_t(endx, endy, endz))).TPC);
      endx -= offset.X();
      endy += offset.Y();
      endz += offset.Z();
      TVector3 dir(endx-vtxx, endy-vtxy, endz-vtxz);
      dir = dir.Unit();
      dirx = dir.X();
      diry = dir.Y();
      dirz = dir.Z();
      // Find the last wire number and peak time on the track
      if (fmthm.isValid()){
        float zlast0=-99999;
        auto vhit=fmthm.at(trackid);
        auto vmeta=fmthm.data(trackid);
        for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
          bool fBadhit = false;
          if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
            fBadhit = true;
            continue;
          }
          if (vmeta[ii]->Index()>=thisTrack->NumberTrajectoryPoints()){
            throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<thisTrack->NumberTrajectoryPoints()<<" for track index "<<trackid<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
          }
          if (!thisTrack->HasValidPoint(vmeta[ii]->Index())){
            fBadhit = true;
            continue;
          }
          auto loc = thisTrack->LocationAtPoint(vmeta[ii]->Index());
          if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
          if (loc.Z()<-100) continue; //hit not on track
          if(vhit[ii]->WireID().Plane==2){
            wirekeys.push_back(vhit[ii].key());
            float zlast=loc.Z();
            if(zlast>zlast0){
              zlast0=zlast;
//              endwire=vhit[ii]->WireID().Wire;
//              endpeakt=vhit[ii]->PeakTime();
//              endtpc=vhit[ii]->WireID().TPC;
            }
          }
        }
      }
      
      // Now we can loop over the hits
      auto const &hits = thass.at(trackid);
      for (auto & hit : hits){
        std::array<float,4> cnn_out = hitResults.getOutput(hit);
        if (hit->WireID().Plane == 2){
          channel.push_back(hit->Channel());
          tpc.push_back(hit->WireID().TPC);
          plane.push_back(hit->WireID().Plane);
          wire.push_back(hit->WireID().Wire);
          charge.push_back(hit->Integral());
          peakt.push_back(hit->PeakTime());
          score_em.push_back(cnn_out[hitResults.getIndex("em")]);
          score_trk.push_back(cnn_out[hitResults.getIndex("track")]);
          score_mic.push_back(cnn_out[hitResults.getIndex("michel")]);
        }
      }
    }
  
    // Get the average of the collection plane scores
    unsigned int nCollectionHits = 0;
    for(unsigned int h = 0; h < plane.size(); ++h){
      if(plane.at(h) == 2){
        ++nCollectionHits;
        average_score_em += score_em.at(h);
        average_score_trk += score_trk.at(h);
        average_score_mic += score_mic.at(h);
      }
    }
  
    if(nCollectionHits > 0){
      average_score_em /= static_cast<double>(nCollectionHits);
      average_score_trk /= static_cast<double>(nCollectionHits);
      average_score_mic /= static_cast<double>(nCollectionHits);
    }
  
    if (!channel.empty())
    {
      std::cout << "Filling output tree" << std::endl;
      ftree->Fill();
    }
  }
}

void pdsp::EMCNNCheckCosmics::beginJob(){

  art::ServiceHandle<art::TFileService> fileServiceHandle;
  ftree = fileServiceHandle->make<TTree>("ftree", "hit info");
  ftree->Branch("run", &run, "run/I");
  ftree->Branch("event", &event, "event/I");
  ftree->Branch("beampdg", &beampdg, "beampdg/I");
  ftree->Branch("average_score_em" , &average_score_em , "average_score_em/D");
  ftree->Branch("average_score_trk", &average_score_trk, "average_score_trk/D");
  ftree->Branch("average_score_mic", &average_score_mic, "average_score_mic/D");
  ftree->Branch("track_endz", &track_endz, "track_endz/D");
  ftree->Branch("ndaughterhits", &ndaughterhits, "ndaughterhits/I");
  ftree->Branch("average_daughter_score_mic", &average_daughter_score_mic, "average_daughter_score_mic/D");
  ftree->Branch("vtxx", &vtxx, "vtxx/D");
  ftree->Branch("vtxy", &vtxy, "vtxy/D");
  ftree->Branch("vtxz", &vtxz, "vtxz/D");
  ftree->Branch("endx", &endx, "endx/D");
  ftree->Branch("endy", &endy, "endy/D");
  ftree->Branch("endz", &endz, "endz/D");
  ftree->Branch("dirx", &dirx, "dirx/D");
  ftree->Branch("diry", &diry, "diry/D");
  ftree->Branch("dirz", &dirz, "dirz/D");
  ftree->Branch("channel", &channel);
  ftree->Branch("tpc", &tpc);
  ftree->Branch("plane", &plane);
  ftree->Branch("wire", &wire);
  ftree->Branch("charge", &charge);
  ftree->Branch("peakt", &peakt);
  ftree->Branch("score_em", &score_em);
  ftree->Branch("score_trk", &score_trk);
  ftree->Branch("score_mic", &score_mic);

  ftree->Branch("daughter_channel", &daughter_channel);
  ftree->Branch("daughter_tpc", &daughter_tpc);
  ftree->Branch("daughter_plane", &daughter_plane);
  ftree->Branch("daughter_wire", &daughter_wire);
  ftree->Branch("daughter_charge", &daughter_charge);
  ftree->Branch("daughter_peakt", &daughter_peakt);
  ftree->Branch("daughter_score_em", &daughter_score_em);
  ftree->Branch("daughter_score_trk", &daughter_score_trk);
  ftree->Branch("daughter_score_mic", &daughter_score_mic);

  ftree->Branch("pdg", &pdg);
  ftree->Branch("origin", &origin);
  ftree->Branch("process", &process);

}


DEFINE_ART_MODULE(pdsp::EMCNNCheckCosmics)
