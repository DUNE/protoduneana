////////////////////////////////////////////////////////////////////////
// Class:       HadronHitsRemoval
// Plugin Type: producer (art v3_06_03)
// File:        HadronHitsRemoval_module.cc
//
// Generated at Tue Jul 13 22:36:12 2021 by Tingjun Yang using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "TMath.h"
#include "TVector3.h"
#include <memory>
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <iostream>

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace pdsp {
  class HadronHitsRemoval;
}
using namespace std;


class pdsp::HadronHitsRemoval : public art::EDProducer {
public:
  explicit HadronHitsRemoval(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HadronHitsRemoval(HadronHitsRemoval const&) = delete;
  HadronHitsRemoval(HadronHitsRemoval&&) = delete;
  HadronHitsRemoval& operator=(HadronHitsRemoval const&) = delete;
  HadronHitsRemoval& operator=(HadronHitsRemoval&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;
  void endJob() override;
  void reset();
  
  bool PassBeamQualityCut() const;
  void TrueBeamInfo(const art::Event & evt, const simb::MCParticle* true_beam_particle);

private:
  // Declare member data here.
  TTree *fTree;
  geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
  double reco_beam_len_sce;
  int selected_track;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  // for beam quality cut of data
  const double beam_startX_data = -27.911;
  const double beam_startY_data = 424.364;
  const double beam_startZ_data = 3.77836;
  const double beam_startX_rms_data = 4.71128;
  const double beam_startY_rms_data = 5.16472;
  const double beam_startZ_rms_data = 1.10265;
  const double beam_angleX_data = 100.454;
  const double beam_angleY_data = 103.523;
  const double beam_angleZ_data = 17.8288;
  
  const double beam_startX_mc = -30.8075;
  const double beam_startY_mc = 422.41;
  const double beam_startZ_mc = 0.11171;
  const double beam_startX_rms_mc = 5.01719;
  const double beam_startY_rms_mc = 4.50862;
  const double beam_startZ_rms_mc = 0.217733;
  const double beam_angleX_mc = 101.578;
  const double beam_angleY_mc = 101.189;
  const double beam_angleZ_mc = 16.5942;

  double reco_beam_calo_startX, reco_beam_calo_startY, reco_beam_calo_startZ;
  double reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ;
  double beam_dx, beam_dy, beam_dz, beam_dxy, beam_costh;
  
  int true_beam_PDG;
  std::string true_beam_endProcess;
  double true_beam_len;
  double true_beam_len_sce;
  // Declare member data here.
  std::string fHitModuleLabel;
  std::string fPFParticleTag;
  std::string fCalorimetryTagSCE;
  std::string fTrackerTag;
  std::string fGeneratorTag;

  double fLimit1;
  double fLimit2Lower, fLimit2Upper;
  double fTruncationLength;
};


pdsp::HadronHitsRemoval::HadronHitsRemoval(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fHitModuleLabel(p.get<std::string>("HitModuleLabel", "hitpdune")),
    fPFParticleTag(p.get<std::string>("PFParticleTag", "pandora")),
    fCalorimetryTagSCE(p.get<std::string>("CalorimetryTagSCE", "pandoracalo")),
    fTrackerTag(p.get<std::string>("TrackerTag", "pandoraTrack")),
    fGeneratorTag(p.get<std::string>("GeneratorTag", "generator")),
    fLimit1(p.get<double>("Limit1", 100.)),
    fLimit2Lower(p.get<double>("Limit2Lower", 230.)),
    fLimit2Upper(p.get<double>("Limit2Upper", 500.)),
    fTruncationLength(p.get<double>("TruncationLength", 10.)) {
  cout<<"$$$HadronHitsRemoval$$$"<<endl;
  // Call appropriate produces<>() functions here.
  recob::HitCollectionCreator::declare_products(producesCollector(), "", true, false);
  produces<std::vector<recob::SpacePoint>>();
  produces<art::Assns<recob::Hit, recob::SpacePoint>>(); // this space point will only be used as a tag
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::HadronHitsRemoval::beginJob(){
  cout<<"$$$beginJob"<<endl;
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");
  fTree->Branch("reco_beam_len_sce", &reco_beam_len_sce);
  fTree->Branch("beam_dx", &beam_dx);
  fTree->Branch("beam_dy", &beam_dy);
  fTree->Branch("beam_dz", &beam_dz);
  fTree->Branch("beam_dxy", &beam_dxy);
  fTree->Branch("beam_costh", &beam_costh);
  fTree->Branch("selected_track", &selected_track);
  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_endProcess", &true_beam_endProcess);
  fTree->Branch("true_beam_len", &true_beam_len);
  fTree->Branch("true_beam_len_sce", &true_beam_len_sce);
}

void pdsp::HadronHitsRemoval::produce(art::Event& evt)
{
  cout<<"#####EvtNo."<<evt.id().event()<<endl;
  reset();
  // Implementation of required member function here.
  // Add code to select beam tracks using Pandora information
  
  recob::HitCollectionCreator hcol(evt, true, false);
  std::unique_ptr<std::vector<recob::SpacePoint>> spvp(new std::vector<recob::SpacePoint>);
  auto assns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
  
  if (!evt.isRealData()) {// MC
    const simb::MCParticle* true_beam_particle = 0x0;
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    if( !true_beam_particle ){
      cout << "No true beam particle" << endl;
    }
    else{
      TrueBeamInfo(evt, true_beam_particle);
      cout<<"true_beam_PDG"<<true_beam_PDG<<endl;
      cout<<"true_beam_endProcess"<<true_beam_endProcess<<endl;
      cout<<"true_beam_len"<<true_beam_len<<endl;
      cout<<"true_beam_len_sce"<<true_beam_len_sce<<endl;
    }
  }
  
  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap
      = pfpUtil.GetPFParticleSliceMap(evt, fPFParticleTag);
  std::vector<std::vector<const recob::PFParticle*>> beam_slices;
  for(auto slice : sliceMap){
    for(auto particle : slice.second){
      bool added = false;
      if(pfpUtil.IsBeamParticle(*particle,evt, fPFParticleTag)){
        if (!added) {
          beam_slices.push_back(slice.second);
          added = true;
        }
      }
      if (!added) {continue;}
      //else {
      //  std::cout << "Not beam particle" << std::endl;
      //}
    }
  }

  std::cout << "Got beam slices: " << beam_slices.size() << std::endl;

  if (beam_slices.size() > 0){
    const recob::PFParticle* particle = beam_slices.at(0).at(0);
    
    const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);

    //Primary Track Calorimetry
    if( thisTrack ){
      auto calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag,
                                                    fCalorimetryTagSCE);
      bool found_calo = false;
      size_t index = 0;
      for ( index = 0; index < calo.size(); ++index) {
        if (calo[index].PlaneID().Plane == 2) {
          found_calo = true;
          break;
        }
      }
      if (found_calo) {
        reco_beam_len_sce = calo[index].Range();
        auto theXYZPoints = calo[index].XYZ();
        std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b)
            {return (a.Z() < b.Z());}); // sort according to Z position
        if (theXYZPoints.size()) {
          reco_beam_calo_startX = theXYZPoints[0].X();
          reco_beam_calo_startY = theXYZPoints[0].Y();
          reco_beam_calo_startZ = theXYZPoints[0].Z();
          reco_beam_calo_endX = theXYZPoints.back().X();
          reco_beam_calo_endY = theXYZPoints.back().Y();
          reco_beam_calo_endZ = theXYZPoints.back().Z();
        }
        
        TVector3 beamdir;
        if (evt.isRealData()) { // beam quality requirement (for data)
          cout<<"@@@Is real data"<<endl;
          beam_dx = (reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
          beam_dy = (reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
          beam_dz = (reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
          beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
          
          TVector3 beamdir_data(cos(beam_angleX_data*TMath::Pi()/180),
                           cos(beam_angleY_data*TMath::Pi()/180),
                           cos(beam_angleZ_data*TMath::Pi()/180));
          beamdir = beamdir_data.Unit();
        }
        else { // beam quality requirement (for MC)
          cout<<"@@@Is MC"<<endl;
          beam_dx = (reco_beam_calo_startX - beam_startX_mc)/beam_startX_rms_mc;
          beam_dy = (reco_beam_calo_startY - beam_startY_mc)/beam_startY_rms_mc;
          beam_dz = (reco_beam_calo_startZ - beam_startZ_mc)/beam_startZ_rms_mc;
          beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
          
          TVector3 beamdir_mc(cos(beam_angleX_mc*TMath::Pi()/180),
                           cos(beam_angleY_mc*TMath::Pi()/180),
                           cos(beam_angleZ_mc*TMath::Pi()/180));
          beamdir = beamdir_mc.Unit();
        }
        
        TVector3 pt0(reco_beam_calo_startX,
                     reco_beam_calo_startY,
                     reco_beam_calo_startZ);
        TVector3 pt1(reco_beam_calo_endX,
                     reco_beam_calo_endY,
                     reco_beam_calo_endZ);
        TVector3 dir = pt1 - pt0;
        dir = dir.Unit();
        beam_costh = dir.Dot(beamdir);

        std::cout << "Length: " << reco_beam_len_sce << std::endl;
        if (reco_beam_len_sce>fLimit1 && PassBeamQualityCut()){
          selected_track = 1;
          if (reco_beam_len_sce>fLimit2Lower && reco_beam_len_sce<fLimit2Upper)
            selected_track = 2;
        }
        std::cout << "Selected: " << selected_track << std::endl;
      }
      
      //string fHitModuleLabel = "hitpdune";
      //string fSpModuleLabel = "reco3d";
      
      //double sel_len = 10.; // shorten to how long (cm)
      
      if (selected_track == 2) { // selected long tracks (not smaller than 5 m)
        auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
        //auto spHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fSpModuleLabel);
        art::FindOneP<recob::Wire>   channelHitWires    (hitsHandle, evt, fHitModuleLabel);
        //art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fSpModuleLabel);
        //recob::HitCollectionCreator hcol(evt, true, true);
        std::vector< art::Ptr<recob::Hit> > eventHits;
        art::fill_ptr_vector(eventHits, hitsHandle);
        //std::unordered_map< size_t, geo::WireID > hitToWire;
        //hitToWire.reserve(eventHits.size());
        
        const std::vector< art::Ptr< recob::Hit > > beamPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *particle, evt, fPFParticleTag );
        auto calo_dQdX = calo[index].dQdx();
        double reco_beam_len = thisTrack->Length();
        cout<<"$$Length"<<reco_beam_len<<endl;
        size_t i = 0;
        for (; thisTrack->Length(i) > reco_beam_len - fTruncationLength/*sel_len*/; ++i ){
          ;
        }
        auto const pos = thisTrack->LocationAtPoint(i);
        cout<<"$$"<<"\tX"<<pos.X()<<"\tY"<<pos.Y()<<"\tZ"<<pos.Z()<<endl; // position of cut point
        
        auto TpIndices = calo[index].TpIndices();
        
        double wirecoord_U = fWireReadoutGeom.Plane(geo::PlaneID{0, 1, 0}).WireCoordinate(pos);
        double wirecoord_V = fWireReadoutGeom.Plane(geo::PlaneID{0, 1, 1}).WireCoordinate(pos);
        double wirecoord_X = fWireReadoutGeom.Plane(geo::PlaneID{0, 1, 2}).WireCoordinate(pos);
        cout<<"$$$WireCoord: U "<<wirecoord_U<<"\tV "<<wirecoord_V<<"\tX "<<wirecoord_X<<endl; // Wire ID of cut point
        
        std::vector< art::Ptr< recob::Hit > > remove_hits;
        std::vector< art::Ptr< recob::Hit > > tag_hits;
        //art::Ptr< recob::SpacePoint > SpPtr;
        std::vector< recob::SpacePoint > spv;
        recob::SpacePoint sp;
        spv.push_back(sp);
        spvp->insert(spvp->end(), spv.begin(), spv.end());
        for (size_t kk = 0; kk < beamPFP_hits.size(); ++kk){
          auto hit = beamPFP_hits[kk];
          
          geo::WireID hitid = hit->WireID();
          //std::vector<geo::WireID> cwids = fWireReadoutGeom.ChannelToWire(hit->Channel());
          //cout<<"@@@"<<hitid.Plane<<"\t"<<hitid.Wire<<"\t"<<hitid.TPC<<"\t"<<hitid.Cryostat<<endl;
          
          if (hitid.TPC == 1 && hitid.Cryostat == 0) {
            if (hitid.Plane == 0) { // plane U
              if (hitid.Wire > wirecoord_U) {
                remove_hits.push_back(hit);
              }
              else {
                tag_hits.push_back(hit);
              }
            }
            else if (hitid.Plane == 1) { // plane V
              if (hitid.Wire < wirecoord_V) {
                remove_hits.push_back(hit);
              }
              else {
                tag_hits.push_back(hit);
              }
            }
            else if (hitid.Plane == 2) { // plane X
              if (hitid.Wire > wirecoord_X) {
                remove_hits.push_back(hit);
              }
              else {
                tag_hits.push_back(hit);
              }
            }
          }
        }
        
        auto const hitPtrMaker = art::PtrMaker<recob::Hit>(evt);
        // fill hits
        for (size_t kk = 0; kk < eventHits.size(); ++kk){
          auto hit = eventHits[kk];
          vector<art::Ptr<recob::Hit>>::iterator itr;
          itr = find(remove_hits.begin(), remove_hits.end(), hit);
          if (itr == remove_hits.end()) { // not be to removed
            hcol.emplace_back(std::move(*hit), channelHitWires.at(kk));
            
            vector<art::Ptr<recob::Hit>>::iterator itr_tag;
            itr_tag = find(tag_hits.begin(), tag_hits.end(), hit);
            if (itr_tag != tag_hits.end()) { // to be tagged (associated)
              auto hitPtr = hitPtrMaker(hcol.size() - 1);
              util::CreateAssn(evt, *spvp, hitPtr, *assns, 0);
            }
          }
        }
        
        cout<<"$$$reco_beam_len_sce "<<reco_beam_len_sce<<endl;
      }
    }
  }
  
  
  hcol.put_into(evt);
  evt.put(std::move(spvp));
  evt.put(std::move(assns));
  fTree->Fill();
}

void pdsp::HadronHitsRemoval::endJob(){
  cout<<"$$$endJob"<<endl;
}

void pdsp::HadronHitsRemoval::reset(){
  cout<<"$$$reset"<<endl;
  reco_beam_len_sce = -999;
  beam_dx = -999;
  beam_dy = -999;
  beam_dz = -999;
  beam_dxy = -999;
  beam_costh = -999;
  reco_beam_calo_startX = -999;
  reco_beam_calo_startY = -999;
  reco_beam_calo_startZ = -999;
  reco_beam_calo_endX = -999;
  reco_beam_calo_endY = -999;
  reco_beam_calo_endZ = -999;
  selected_track = 0;
  
  true_beam_PDG = -999;
  true_beam_endProcess = "";
  true_beam_len = -999.;
  true_beam_len_sce = -999.;
}

bool pdsp::HadronHitsRemoval::PassBeamQualityCut() const{ // cut on beam entrance location and beam angle

  if (beam_dxy<-1) return false;
  if (beam_dxy>3) return false;
  
  if (beam_dz<-3) return false;
  if (beam_dz>3) return false;
  
  if (beam_costh<0.95) return false;
  if (beam_costh>2) return false;
  
  return true;
}

void pdsp::HadronHitsRemoval::TrueBeamInfo( const art::Event & evt, const simb::MCParticle* true_beam_particle)
{
  true_beam_PDG = true_beam_particle->PdgCode();
  true_beam_endProcess = true_beam_particle->EndProcess();
  const simb::MCTrajectory & true_beam_trajectory = true_beam_particle->Trajectory();
  //auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  
  vector<double> true_beam_traj_X;
  vector<double> true_beam_traj_Y;
  vector<double> true_beam_traj_Z;
  vector<double> true_beam_traj_X_SCE;
  vector<double> true_beam_traj_Y_SCE;
  vector<double> true_beam_traj_Z_SCE;
  for (size_t i = 0; i < true_beam_trajectory.size(); ++i) {
    double trajX = true_beam_trajectory.X(i);
    double trajY = true_beam_trajectory.Y(i);
    double trajZ = true_beam_trajectory.Z(i);
    true_beam_traj_X.push_back(trajX);
    true_beam_traj_Y.push_back(trajY);
    true_beam_traj_Z.push_back(trajZ);
    /*auto offset = sce->GetPosOffsets( {trajX, trajY, trajZ} );
    cout<<"@@@ "<<i<<endl;
    true_beam_traj_X_SCE.push_back(trajX - offset.X());
    true_beam_traj_Y_SCE.push_back(trajY + offset.Y());
    true_beam_traj_Z_SCE.push_back(trajZ + offset.Z());*/
  }
  true_beam_len = 0;
  true_beam_len_sce = 0;
  for (size_t i = 1; i < true_beam_trajectory.size(); ++i) {
    true_beam_len += sqrt( pow( true_beam_traj_X.at(i)-true_beam_traj_X.at(i-1), 2)
                       + pow( true_beam_traj_Y.at(i)-true_beam_traj_Y.at(i-1), 2)
                       + pow( true_beam_traj_Z.at(i)-true_beam_traj_Z.at(i-1), 2));
    /*true_beam_len_sce += sqrt( pow( true_beam_traj_X_SCE.at(i)-true_beam_traj_X_SCE.at(i-1), 2)
                       + pow( true_beam_traj_Y_SCE.at(i)-true_beam_traj_Y_SCE.at(i-1), 2)
                       + pow( true_beam_traj_Z_SCE.at(i)-true_beam_traj_Z_SCE.at(i-1), 2));*/
  }
}

DEFINE_ART_MODULE(pdsp::HadronHitsRemoval)
