////////////////////////////////////////////////////////////////////////
// Class:       CosmicsdQdx
// Plugin Type: analyzer (art v3_05_01)
// File:        CosmicsdQdx_module.cc
//
// Generate simple summary tree from track dQ/dx measurements for each
// collection view to check track calorimetry and effecit CRP gains
// Follow the same logic as in TrackCalorimetryAlg
//
// Generated at Tue May 26 16:28:38 2020 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <utility>

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

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"

//
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEDPCRPGeo.h"

// ROOT
#include "TTree.h"
//#include "TMath.h"
//#include "TH1F.h"
//#include "TFile.h"

//
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

//
namespace pddpana {
  class CosmicsdQdx;
}

//
namespace {
  //from calo::TrackCalorimetryAlg
  class dist_projected{
  public:
    dist_projected(recob::Hit const* h, geo::WireReadoutGeom const* cm):
      hit(h), wireReadout(cm){}
    bool operator() (std::pair<geo::WireID,float> const& i, std::pair<geo::WireID,float> const& j)
    {
      float dw_i = ((int)(i.first.Wire) - (int)(hit->WireID().Wire))*wireReadout->Plane(i.first).WirePitch();
      float dw_j = ((int)(j.first.Wire) - (int)(hit->WireID().Wire))*wireReadout->Plane(j.first).WirePitch();
      float dt_i = i.second - hit->PeakTime();
      float dt_j = j.second - hit->PeakTime();
      return (std::sqrt(dw_i*dw_i + dt_i*dt_i) < std::sqrt(dw_j*dw_j + dt_j*dt_j));
    }
  private:
    recob::Hit const* hit;
    geo::WireReadoutGeom const* wireReadout;
  };
}


class pddpana::CosmicsdQdx : public art::EDAnalyzer {
public:
  explicit CosmicsdQdx(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicsdQdx(CosmicsdQdx const&) = delete;
  CosmicsdQdx(CosmicsdQdx&&) = delete;
  CosmicsdQdx& operator=(CosmicsdQdx const&) = delete;
  CosmicsdQdx& operator=(CosmicsdQdx&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  int      fLogLevel;
  string   fTrackModuleLabel;
  float    fTrackMinLen;
  float    fTrackDriftCut;
  float    fTrackWallCut;
  float    fTrackLemCut;
  unsigned fMaxHitMultiplicity;

  unsigned fTotalTracks;
  unsigned fSelectedTracks;
  unsigned fHitCutTracks;
  unsigned fGeoCutTracks;

  // summary tree
  TTree   *fTree;

  //
  unsigned fEventNum;
  unsigned fTrackId;
  unsigned fTrajPoints;
  float    fTrackLen;
  float    fDriftOffset;
  float    fStartX;
  float    fStartY;
  float    fStartZ;
  float    fStartNx;
  float    fStartNy;
  float    fStartNz;
  float    fEndX;
  float    fEndY;
  float    fEndZ;
  float    fEndNx;
  float    fEndNy;
  float    fEndNz;

  //
  unsigned fPlane;
  float    fHitAdcSum;
  float    fHitIntegral;
  float    fHitPeak;
  float    fHitTime;
  float    fPitch;
  float    fdQdx;
  float    fX;
  float    fY;
  float    fZ;
  float    fQAsym;

  // plane hits
  vector< vector< const recob::Hit* > > fPlaneHits;

  // track utils
  protoana::ProtoDUNETrackUtils trackUtil;

  // CRP Geo util
  protoana::ProtoDUNEDPCRPGeo crpGeoUtil;

  // detector geometry
  const geo::Geometry* fGeom;
  const geo::WireReadoutGeom* fWireReadoutGeom;

  bool checkCutsAndGetT0( const recob::Track& track, float &T0 );
  bool checkTrackHitInfo( const recob::Track& track, art::Event const& e);
};


pddpana::CosmicsdQdx::CosmicsdQdx(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fTrackModuleLabel( p.get< std::string  >("TrackModuleLabel") ),
  fTrackMinLen( p.get< float  >("TrackMinLen") ),
  fTrackDriftCut( p.get< float  >("TrackDriftCut") ),
  fTrackWallCut( p.get< float  >("TrackWallCut") ),
  fTrackLemCut( p.get< float  >("TrackLemCut") ),
  fMaxHitMultiplicity( p.get< float  >("MaxHitMultiplicity") )
  {
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
    fWireReadoutGeom = &art::ServiceHandle<geo::WireReadout>()->Get();

    //auto const &tpc = ;
    auto const [axis, _] = fGeom->TPC().DriftAxisWithSign();
    if( axis != geo::Coordinate::X ){
      throw cet::exception("CosmicsdQdx") << "Drift direction "<<axis
                                          << " is not support\n";
    }
  }

//
void pddpana::CosmicsdQdx::analyze(art::Event const& e)
{
  const string myname = "pddpana::CosmicsdQdx::analyze: ";

  // get tracks
  auto Tracks   = e.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
  fTotalTracks += Tracks->size();

  if( fLogLevel >= 2 ){
    cout<<myname<<"The event contains "<< Tracks->size() <<" tracks\n";
  }

  //
  fEventNum = e.id().event();

  // loop over tracks
  for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {
    const recob::Track& track = Tracks->at(itrk);

    fTrackId     = itrk;
    fTrajPoints  = track.NumberTrajectoryPoints();
    fTrackLen    = track.Length();

    //
    fStartX  = track.Vertex().X();
    fStartY  = track.Vertex().Y();
    fStartZ  = track.Vertex().Z();

    fEndX    = track.End().X();
    fEndY    = track.End().Y();
    fEndZ    = track.End().Z();

    fStartNx = track.VertexDirection().X();
    fStartNy = track.VertexDirection().Y();
    fStartNz = track.VertexDirection().Z();

    fEndNx   = track.EndDirection().X();
    fEndNy   = track.EndDirection().Y();
    fEndNz   = track.EndDirection().Z();

    fQAsym   = 0;
    if( !checkTrackHitInfo( track, e ) ){
      if( fLogLevel >= 2 ){
        cout<<myname<<"Skipping track ID "<<fTrackId
            <<" due to hit plane mismatch "<<endl;
      }
      fHitCutTracks++;
      continue;
    }

    fDriftOffset = 0;
    if( !checkCutsAndGetT0( track, fDriftOffset) ){
      if( fLogLevel >= 2 ){
        cout<<myname<<"Skipping track ID "<<fTrackId<<" has "
            <<fTrajPoints<<" points and is "<<fTrackLen<<" cm long\n";
      }
      fGeoCutTracks++;
      continue;
    }

    //
    if( fLogLevel >= 3 ){
      cout<<myname<<"Track ID "<<fTrackId<<" has "
          <<fTrajPoints<<" points and is "<<fTrackLen<<" cm long"
          << "\n  start at: ( " << track.Vertex().X()
          << " ; " << track.Vertex().Y()
          << " ; " << track.Vertex().Z()
          << "\n  end at:   ( " << track.End().X()
          << " ; " << track.End().Y()
          << " ; " << track.End().Z()<<" )"<<endl;
    }

    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
    // loop over the planes
    for(size_t i_plane=0; i_plane<fWireReadoutGeom->Nplanes(); i_plane++) {
      fPlane = i_plane;

      // get hits in this plane
      //auto hits = trackUtil.GetRecoTrackHitsFromPlane( track, e, fTrackModuleLabel, i_plane );
      auto hits = fPlaneHits[ i_plane ]; // loaded in checkTrackHitInfo
      if( fLogLevel >= 3 ){
        cout<<myname<<"Hits in plane "<<i_plane<<" "<<hits.size()<<endl;
      }
      //

      //project down the track into wire/tick space for this plane
      vector< unsigned > traj_points_idx;
      vector< pair<geo::WireID,float> > traj_points_in_plane; //(track.NumberTrajectoryPoints());
      for(size_t i_trjpt=0; i_trjpt<track.NumberTrajectoryPoints(); i_trjpt++){
        auto pnt        = track.LocationAtPoint(i_trjpt);
        if( !track.HasValidPoint( i_trjpt ) ){
          if( fLogLevel>=4 ){
            cerr<<"track point "<<i_trjpt<<" has position ("
                <<pnt.X()<<", "<<pnt.Y()<<", "<<pnt.Z()<<") \n"
                <<" and is not a valid point "<<track.HasValidPoint( i_trjpt )<<endl;
          }
          continue;
        }

        double x_pos    = pnt.X();
        float tick      = detProp.ConvertXToTicks(x_pos,(int)i_plane,0,0);
        try{
          auto tpcid = fGeom->PositionToTPCID(pnt);
          // if( tpcid.TPC >= fGeom->NTPC() ){
          //   if( fLogLevel>=3 ){
          //     cerr<<myname<<"tpc no "<<tpcid.TPC<<" is not valid \n"
          //      <<" track id "<<fTrackId<<" at ("
          //      <<pnt.X()<<", "<<pnt.Y()<<", "<<pnt.Z()<<") \n"
          //      <<" has valid point "<<track.HasValidPoint( i_trjpt )<<endl;
          //   }
          //   continue;
          // }
          geo::WireID wid = fWireReadoutGeom->Plane(geo::PlaneID(tpcid, i_plane)).NearestWireID(pnt);
          //traj_points_in_plane[i_trjpt] = std::make_pair(wid, tick);
          traj_points_in_plane.push_back( std::make_pair(wid, tick) );
          traj_points_idx.push_back( i_trjpt );
        }
        catch(cet::exception &e){
          if( fLogLevel>=4 ){
            cout<<e;
          }
          continue;
        }
      }

      //
      // from calo::TrackCalorimetryAlg::AnalyzeHit
      for(auto const &hit: hits ){
        //skip high mulitplicity hits
        if(hit->Multiplicity() > (int)fMaxHitMultiplicity) continue;

        //
        size_t traj_iter = std::distance( traj_points_in_plane.begin(),
                                          std::min_element( traj_points_in_plane.begin(),
                                                            traj_points_in_plane.end(),
                                                            dist_projected(hit, fWireReadoutGeom) ) );
        size_t traj_pnt_idx = traj_points_idx[traj_iter];
        try{
          // fGeom->View(hit->WireID().Plane) always returns 3 (geo::kY) for some reason
          //fPitch = lar::util::TrackPitchInView(track, fGeom->View(hit->WireID().Plane), traj_pnt_idx);
          fPitch = lar::util::TrackPitchInView(track, hit->View(), traj_pnt_idx);
        }
        catch(cet::exception &e){
          if( fLogLevel>=3 ){
            cout<<e;
          }
          continue;
        }

        // hit properties
        fHitAdcSum   = hit->SummedADC();
        fHitIntegral = hit->Integral();
        fHitPeak     = hit->PeakAmplitude();
        fHitTime     = hit->PeakTime();
        fdQdx        = fHitAdcSum / fPitch;

        auto pnt = track.LocationAtPoint(traj_iter);
        fX = pnt.X();
        fY = pnt.Y();
        fZ = pnt.Z();

        // if( fEventNum == 11540 && fTrackId == 0 ){
        //   auto dir1 = track.DirectionAtPoint(traj_iter);
        //   cout<<i_plane<<" "<<fPitch<<" "<<fX<<" "<<fY<<" "<<fZ<<" "
        //       <<dir1.X()<<" "<<dir1.Y()<<" "<<dir1.Z()<<endl;
        //}

        fTree->Fill();
      }// loop over track plane hits

    }// end plane loop

    //
    fSelectedTracks++;
  }// end track loop

} // end analyze()

//
void pddpana::CosmicsdQdx::beginJob()
{
  fPlaneHits.resize( fWireReadoutGeom->Nplanes() );

  fTotalTracks    = 0;
  fSelectedTracks = 0;
  fHitCutTracks   = 0;
  fGeoCutTracks   = 0;

  // init summary tree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("cranaTree","Check cosmics dQdx");
  fTree->Branch("EventNum", &fEventNum, "EventNum/i");
  fTree->Branch("TrackId", &fTrackId, "TrackId/i");
  fTree->Branch("TrajPoints", &fTrajPoints, "TrajPoints/i");
  fTree->Branch("TrackLen",   &fTrackLen,   "TrackLen/F");
  fTree->Branch("DriftOffset",  &fDriftOffset,  "DriftOffset/F");
  fTree->Branch("StartX",  &fStartX,  "StartX/F");
  fTree->Branch("StartY",  &fStartY,  "StartY/F");
  fTree->Branch("StartZ",  &fStartZ,  "StartZ/F");
  fTree->Branch("StartNx",  &fStartNx,  "StartNx/F");
  fTree->Branch("StartNy",  &fStartNy,  "StartNy/F");
  fTree->Branch("StartNz",  &fStartNz,  "StartNz/F");
  fTree->Branch("EndX",  &fEndX,  "EndX/F");
  fTree->Branch("EndY",  &fEndY,  "EndY/F");
  fTree->Branch("EndZ",  &fEndZ,  "EndZ/F");
  fTree->Branch("EndNx",  &fEndNx,  "EndNx/F");
  fTree->Branch("EndNy",  &fEndNy,  "EndNy/F");
  fTree->Branch("EndNz",  &fEndNz,  "EndNz/F");

  fTree->Branch("Plane", &fPlane, "Plane/i");
  fTree->Branch("HitAdcSum", &fHitAdcSum, "HitAdcSum/F");
  fTree->Branch("HitIntegral", &fHitIntegral, "HitIntegral/F");
  fTree->Branch("HitPeak", &fHitPeak, "HitPeak/F");
  fTree->Branch("HitTime", &fHitTime, "HitTime/F");
  fTree->Branch("Pitch", &fPitch, "Pitch/F");
  fTree->Branch("dQdx", &fdQdx, "dQdx/F");
  fTree->Branch("QAsym", &fQAsym, "QAsym/F");
  fTree->Branch("X", &fX, "X/F");
  fTree->Branch("Y", &fY, "Y/F");
  fTree->Branch("Z", &fZ, "Z/F");
}

//
void pddpana::CosmicsdQdx::endJob()
{
  const string myname = "pddpana::CosmicsdQdx::endJob: ";
  cout<<myname<<"tracks processed total     : "<<fTotalTracks<<endl;
  cout<<myname<<"tracks processed selected  : "<<fSelectedTracks<<endl;
  cout<<myname<<"tracks cut due to hit info : "<<fHitCutTracks<<endl;
  cout<<myname<<"tracks cut due to geo info : "<<fGeoCutTracks<<endl;

}

//
//
bool pddpana::CosmicsdQdx::checkCutsAndGetT0( const recob::Track& track, float &T0 )
{
  const string myname = "pddpana::CosmicsdQdx::checkCutsAndGetT0: ";

  //
  T0 = 0; // track distance in cm from the top

  if( track.Length() < fTrackMinLen ){
    if( fLogLevel>=3 ){
      cout<<myname<<"Failed to pass length cut "<<track.Length()<<endl;
    }
    return false;
  }

  auto tstart = track.Vertex();
  auto tend   = track.End();

  // cosmics are downward going and anode is at +300cm
  if( tstart.X() < tend.X() ){
    tstart = tend;
  }

  auto crp_geo_info = crpGeoUtil.GetCRPGeoInfo( tstart );
  if( not crp_geo_info.valid ) return false;

  // drift distance to the anode
  T0 = crp_geo_info.danode;
  if( T0 < fTrackDriftCut ){
    if( fLogLevel>=3 ){
      cout<<myname<<"Failed to pass T0 cut "<<T0<<endl;
    }
    return false;
  }

  // distance to CRP border
  if( crp_geo_info.dedge < fTrackWallCut ){
    if( fLogLevel>=3 ){
      cout<<myname<<"Failed to pass wall cut ("
          <<tstart.X()<<", "<<tstart.Y()<<", "<<tstart.Z()<<")"<<endl;
    }

    return false;
  }

  // distance to LEM border
  if( crp_geo_info.dlem < fTrackLemCut ){
    if( fLogLevel>=3 ){
      cout<<myname<<"Failed to pass LEM edge cut ("
          <<tstart.X()<<", "<<tstart.Y()<<", "<<tstart.Z()<<")"<<endl;
    }

    return false;
  }

  return true;
}

//
// check that hits associated to a track belong to the same CRP
bool pddpana::CosmicsdQdx::checkTrackHitInfo( const recob::Track& track, art::Event const& e ){
  const string myname = "pddpana::CosmicsdQdx::checkTrackHitInfo: ";

  //
  unsigned EventId = e.id().event();


  // loop over the planes and get all the hits
  for(size_t i_plane = 0; i_plane<fWireReadoutGeom->Nplanes(); i_plane++) {
    // get hits in this plane
    auto hits = trackUtil.GetRecoTrackHitsFromPlane( track, e, fTrackModuleLabel, i_plane );
    fPlaneHits[i_plane] = hits;
  }

  //
  // calculate charge asymmetry
  //
  // charge in view 0
  double v0 = 0;
  for( auto const &h: fPlaneHits[0] ){
    v0 += h->SummedADC();
  }
  // charge in view 1
  double v1 = 0;
  for( auto const &h: fPlaneHits[1] ){
    v1 += h->SummedADC();
  }
  // charge asymmetry
  fQAsym = (v0 - v1)/(v0 + v1);

  //
  //
  // check for TPC mismatches
  vector<unsigned> hitsTpcId( fWireReadoutGeom->Nplanes() );

  for(size_t i_plane = 0; i_plane<fWireReadoutGeom->Nplanes(); i_plane++) {
    auto hits = fPlaneHits[i_plane];

    vector<unsigned> plane_hits_tpcid( hits.size() );
    // loop over hits
    for(size_t i_hit = 0; i_hit<hits.size(); i_hit++ ){
      plane_hits_tpcid[i_hit] = hits[i_hit]->WireID().TPC;
    }

    int this_tpcid = -1;
    auto start = plane_hits_tpcid.begin();
    auto end   = plane_hits_tpcid.end();
    if( std::equal(start + 1, end, start)) {
      this_tpcid = (int)(*start);
    } //
    else {
      if( fLogLevel>= 3){
        cout<<myname<<"hits for the same plane have mixed TPC IDs. Skipping...\n";
      }
    }
    if( this_tpcid < 0 ){
      hitsTpcId.clear();
      break;
    }
    else{
      hitsTpcId[i_plane] = (unsigned)this_tpcid;
    }
  }// end plane loop

  if( hitsTpcId.empty() ) return true;

  auto start = hitsTpcId.begin();
  auto end   = hitsTpcId.end();
  if( !std::equal(start + 1, end, start)) {
    if( fLogLevel >= 3){
      cout<<myname<<"mismatch in TPC ID for the initial hits: ";
      cout<<" event ID "<<EventId<<"; hit TPC IDs ";
      for( auto const &v: hitsTpcId ){ cout<<v<<" ";}
      cout<<endl;
    }
    return false;
  } //

  return true;
}


DEFINE_ART_MODULE(pddpana::CosmicsdQdx)
