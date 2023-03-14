////////////////////////////////////////////////////////////////////////
// Class:       CheckHitsAndTracks
// Plugin Type: analyzer (art v3_05_01)
// File:        CheckHitsAndTracks_module.cc
//
// Generate simple summary tree from tracks and hits to check sim and reco
// Follow the same logic as in CosmicsdQdx
//
// Generated at 31/01/2023 by Thibaut Houdy
//////////////////////////////////////////////////////////////////////Check//
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
#include "canvas/Persistency/Common/FindManyP.h"

#include "art_root_io/TFileService.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"

#include "TPolyMarker3D.h"
#include "TPolyMarker.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TLine.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

// ROOT
#include "TTree.h"
#include "TLegend.h"

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

//
namespace pdvdana {
  class CheckHitsAndTracks;
}

class pdvdana::CheckHitsAndTracks : public art::EDAnalyzer {
public:
  explicit CheckHitsAndTracks(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CheckHitsAndTracks(CheckHitsAndTracks const&) = delete;
  CheckHitsAndTracks(CheckHitsAndTracks&&) = delete;
  CheckHitsAndTracks& operator=(CheckHitsAndTracks const&) = delete;
  CheckHitsAndTracks& operator=(CheckHitsAndTracks&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  int      fLogLevel;
  int      fFlagBasics = 1;
  int      fFlagInfos = 2;
  int      fFlagWarning = 3;
  int      fFlagDetails = 4;
  //int      fFlagDebug = 5;

  string   fHitModuleLabel;
  string   fTrackModuleLabel;
  float    fTrackMinLen;
  int      ev_num;

  unsigned fTotalTracks;
  unsigned fSelectedTracks;

  // summary tree
  TTree   *ftrackTree;
  TTree   *fhitTree;

  //
  unsigned fEventNum;
  unsigned fTrackId;
  unsigned fNtpcs, fNplanes;
  unsigned ftrackStartTick;
  unsigned ftrackEndTick;

  float ftrackLen;
  float ftrackDx;
  float ftrackDy;
  float ftrackDz;
  float ftrackStartX;
  float ftrackStartY;
  float ftrackStartZ;
  float ftrackEndX;
  float ftrackEndY;
  float ftrackEndZ;

  float ftracktheta;
  float ftrackphi;
  float ftracknorm;

  float fhitX;
  float fhitY;
  float fhitZ;
  unsigned fhitPlane, fhitTPC, fhitWire, fhitChannel;
  int fhitTime;

  vector<int> fOffsetWireID_u;
  vector<int> fOffsetWireID_v;
  vector<int> fOffsetWireID_z;

  float fgeoXmin = 1e6;
  float fgeoXmax =-1e6;
  float fgeoYmin = 1e6;
  float fgeoYmax =-1e6;
  float fgeoZmin = 1e6;
  float fgeoZmax =-1e6;

  // detector geometry
  const geo::Geometry* fGeom;
  const geo::WireReadoutGeom* fWireReadoutGeom;

  void checkTrackCharacs( float Dx, float Dy, float Dz, vector<float> &track_char );
  void DrawRectangle(TPolyLine *Rectangle, float a_1, float a_2, float b_1, float b_2);
  void DrawCube(TPolyLine3D *Cube, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max );
  void Drawing3D_HitsAndTracks(TCanvas* l_c, vector<vector<float>> l_t, vector<vector<float>> l_v);
  void DrawingTPC_HitsAndTracks(TCanvas* l_c, vector<vector<int>> l_a, vector<vector<int>> l_b, vector<vector<int>> l_d,  vector<vector<int>> l_e);
  void Drawing2D_HitsAndTracks(TCanvas* l_c, vector<vector<float>> l_t, vector<vector<float>> l_v);
  void SetHitStyle(TPolyMarker *l_hits);
  unsigned GetWireOffset(unsigned plane_id, unsigned tpc_id);

};


pdvdana::CheckHitsAndTracks::CheckHitsAndTracks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fHitModuleLabel( p.get< std::string >("HitModuleLabel") ),
  fTrackModuleLabel( p.get< std::string >("TrackModuleLabel") ),
  fTrackMinLen( p.get< float  >("TrackMinLen") )
  {
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
    fWireReadoutGeom = &art::ServiceHandle<geo::WireReadout>()->Get();
  }

//
void pdvdana::CheckHitsAndTracks::analyze(art::Event const& ev)
{

  const string myname = "pdvdana::CheckHitsAndTracks::analyze: ";

  // get tracks

  string hittag = "gaushit";
  art::InputTag hit_tag(fHitModuleLabel);
  string trcktag = "pandoraTrack";
  art::InputTag trck_tag(trcktag);
  string hittag_dis = "hitpdune";
  art::InputTag hit_tag_dis(hittag_dis);
  string calotag = "pandoraGnocalo";//"calo";
  art::InputTag calo_tag(calotag);
  string clutag = "pandora";
  art::InputTag cluster_tag(clutag);

  vector<vector<float>> trackSpacePoints;
  vector<vector<float>> vertexSpacePoints;

  vector<vector<int>> vector_Hits_tpc_channel(fNplanes);
  vector<vector<int>> vector_Hits_tpc_time(fNplanes);
  vector<vector<int> > fTrackStartWireID(fNplanes);
  vector<vector<int> > fTrackEndWireID(fNplanes);

  ev_num++;

  auto Tracks = ev.getValidHandle<vector<recob::Track>>(trck_tag);
  auto Hits = ev.getValidHandle<vector<recob::Hit>>(hit_tag);

  const art::FindManyP<recob::Track> fmp(Hits ,ev ,trcktag);

  fTotalTracks += Tracks->size();
  fEventNum = ev_num;

  vector<float> track_charac(6);

  auto const hitHandle = ev.getValidHandle<std::vector<recob::Hit>>(hit_tag);
  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits, hitHandle);

  if (fLogLevel>=fFlagInfos){
    cout<<myname<<Hits->size()<<" hits found"<<endl;
    cout<<myname<<Tracks->size()<<" tracks found"<<endl;
  }

  // Get the hits associated with the space points
  const art::FindManyP<recob::SpacePoint> fmsph(hitHandle, ev, cluster_tag);
  if (!fmsph.isValid()) {
    throw cet::exception("LArPandoraShowerCheatingAlg")
      << "Spacepoint and hit association not valid. Stopping.";
  }

  //----------------------------------------------------------------------------------------------
  //Loops on tracks
  for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {
    const recob::Track& track = Tracks->at(itrk);
    ftrackLen    = track.Length();

    if (ftrackLen<fTrackMinLen)
      continue;
    ftrackStartX  = track.Start().X();
    ftrackStartY  = track.Start().Y();
    ftrackStartZ  = track.Start().Z();
    ftrackStartTick  = 0;
    ftrackEndX    = track.End().X();
    ftrackEndY    = track.End().Y();
    ftrackEndZ    = track.End().Z();
    ftrackEndTick = 0;

    // Filling the track tree
    checkTrackCharacs(ftrackEndX-ftrackStartX, ftrackEndY-ftrackStartY, ftrackEndZ-ftrackStartZ, track_charac);
    ftrackDx = track_charac[0];
    ftrackDy = track_charac[1];
    ftrackDz = track_charac[2];
    ftracknorm = track_charac[3];
    ftracktheta = track_charac[4];
    ftrackphi = track_charac[5];
    fTrackId     = 0;
    ftrackTree->Fill();

    if (fLogLevel>=fFlagDetails){
      cout <<myname<< "-- Track: "<<itrk<<" x,y,z: (" <<ftrackStartX<<"; "<< ftrackStartY<<"; "<< ftrackStartZ<<" ) | x,y,z: ("<< ftrackEndX<<"; "<< ftrackEndY<<"; "<< ftrackEndZ<<")"<<endl;
    }

    geo::TPCID tpc_s = fGeom->FindTPCAtPosition(track.Start());
    geo::TPCID tpc_e = fGeom->FindTPCAtPosition(track.End());

    int wireID_s = -9999;
    int wireID_e = -9999;

    if(! (tpc_s.isValid && tpc_e.isValid)){
      if (fLogLevel>=fFlagWarning){
        cout<<"---!!--- TPCs are not valid : "<<endl;
        cout<<"---!!--- TPCs ID: "<<tpc_s.TPC <<" and "<<tpc_e.TPC << endl;
      }
      continue;
    }
    else{
      geo::WireID g_wireID_s;
      geo::WireID g_wireID_e;
      for(unsigned int id_pl = 0; id_pl<fNplanes; id_pl++){
        if ((tpc_s.TPC != tpc_e.TPC) && (fLogLevel>=fFlagWarning)){
          cout<<"---> Carefull, start and end are in different TPCs"<<endl;
        }
        try{
          g_wireID_s = fWireReadoutGeom->NearestWireID(track.Start(), geo::PlaneID(0, tpc_s.TPC, id_pl));
        }
        catch(geo::InvalidWireError const& e_1) {
          g_wireID_s = e_1.suggestedWireID(); // pick the closest valid wire
        }
        //Looking for the wire associated with the end of the track
        try{
          g_wireID_e = fWireReadoutGeom->NearestWireID(track.End(), geo::PlaneID(0, tpc_e.TPC, id_pl));
        }
        catch(geo::InvalidWireError const& e_2) {
          g_wireID_e = e_2.suggestedWireID(); // pick the closest valid wire
        }

        wireID_s = g_wireID_s.Wire;
        wireID_e = g_wireID_e.Wire;

        if ((wireID_e!=-9999) && (wireID_s!=-9999)){

          wireID_s += GetWireOffset(id_pl, tpc_s.TPC);
          wireID_e += GetWireOffset(id_pl, tpc_e.TPC);
          if (fLogLevel>=fFlagDetails)
            cout<<"---- START - TPC: "<<tpc_s.TPC <<" wire ID: "<< wireID_s<<" | END : "<<tpc_e.TPC << " wire ID: "<< wireID_e << endl;
          fTrackStartWireID[id_pl].push_back(g_wireID_s.Wire + GetWireOffset(id_pl, tpc_s.TPC));
          fTrackEndWireID[id_pl].push_back( g_wireID_e.Wire + GetWireOffset(id_pl, tpc_e.TPC));
          fSelectedTracks++;
        }
      }
    }
    vector<float> pos_track = {ftrackStartX, ftrackStartY, ftrackStartZ, ftrackEndX, ftrackEndY, ftrackEndZ};
    trackSpacePoints.push_back(pos_track);
  }// Tracks loop

  //----------------------------------------------------------------------------------------------
  //Loop on hits
  for (auto hit : hits) {
    std::vector<art::Ptr<recob::Track>> tracks = fmp.at(hit.key());
    std::vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(hit.key());

    if (sps.size() == 1 && tracks.size()==1) {
      art::Ptr<recob::SpacePoint> sp = sps.front();
      fhitX = sp->XYZ()[0];
      fhitY = sp->XYZ()[1];
      fhitZ = sp->XYZ()[2];
      fTrackId     = 0;

      geo::WireID hit_wireID = hit->WireID();
      fhitTPC = hit_wireID.TPC;
      fhitPlane = hit_wireID.Plane;
      fhitWire = hit_wireID.Wire+GetWireOffset(fhitPlane, fhitTPC);
      fhitTime = hit->PeakTime();
      fhitChannel = hit->Channel();
      if (fhitTime>10000)
        continue;

      std::vector<geo::WireID> cwids = fWireReadoutGeom->ChannelToWire(fhitChannel);
      vector_Hits_tpc_channel[fhitPlane].push_back(fhitChannel);
      vector_Hits_tpc_time[fhitPlane].push_back(fhitChannel);
      vector<float> pos_sp = {fhitX, fhitY, fhitZ};
      vertexSpacePoints.push_back(pos_sp);
      fhitTree->Fill();
      }
  }// Hits loop


  //----------------------------------------------------------------------------------------------
  // Drawing


  art::ServiceHandle<art::TFileService> tfs;

  gStyle->SetOptStat(0);

  //Drawing 3D tracks
  TString canvasName_3D = Form("canvas3D_%i", ev_num);
  TCanvas* canvas_3D = tfs->make<TCanvas>(canvasName_3D, canvasName_3D);
  Drawing3D_HitsAndTracks(canvas_3D, trackSpacePoints, vertexSpacePoints);
  canvas_3D->Write(canvasName_3D);


  TString canvasName_xy = Form("canvasxy_%i", ev_num);
  TCanvas* canvas_xy = tfs->make<TCanvas>(canvasName_xy, canvasName_xy);
  Drawing2D_HitsAndTracks(canvas_xy, trackSpacePoints, vertexSpacePoints);
  canvas_xy->Write(canvasName_xy);


  //Drawing TPC canvas -> comparison hits (wire, time) with start/end tracks (only wire, time is arbitrary)
  TString canvasName_tpc = Form("canvastpc_%i", ev_num);
  TCanvas* canvas_tpc = tfs->make<TCanvas>(canvasName_tpc, canvasName_tpc);
  DrawingTPC_HitsAndTracks(canvas_tpc, vector_Hits_tpc_channel, vector_Hits_tpc_time, fTrackStartWireID, fTrackEndWireID);
  canvas_tpc->Write(canvasName_tpc);


}
//
void pdvdana::CheckHitsAndTracks::beginJob()
{

  ev_num = 0;
  fTotalTracks    = 0;
  fSelectedTracks = 0;
  fNplanes = fWireReadoutGeom->Nplanes();
  fNtpcs = fGeom->NTPC();

  fOffsetWireID_u.clear();
  fOffsetWireID_v.clear();
  fOffsetWireID_z.clear();

  unsigned int mem_nwires_u = 0;
  unsigned int mem_nwires_v = 0;
  unsigned int mem_nwires_z = 0;

  if(fLogLevel>=fFlagBasics){
    cout<<"Geometry : "<<endl;
    cout<<"number of planes :  "<<fNplanes<<endl;
    cout<<"number of tpc :  "<<fNtpcs<<endl;
  }

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){

    geo::TPCID tpcid{0, t_tpc_id};
    geo::PlaneID const uplane_id{tpcid, geo::View_t::kU};
    geo::PlaneID const vplane_id{tpcid, geo::View_t::kV};
    geo::PlaneID const zplane_id{tpcid, geo::View_t::kZ};
    unsigned int nwires_u = fWireReadoutGeom->Nwires(uplane_id);
    unsigned int nwires_v = fWireReadoutGeom->Nwires(vplane_id);
    unsigned int nwires_z = fWireReadoutGeom->Nwires(zplane_id);

    fOffsetWireID_u.push_back(mem_nwires_u);
    fOffsetWireID_v.push_back(mem_nwires_v);
    fOffsetWireID_z.push_back(mem_nwires_z);

    mem_nwires_u += nwires_u;
    mem_nwires_v += nwires_v;
    mem_nwires_z += nwires_z;

    if(fLogLevel>=fFlagInfos){
      std::cout << "  TPC " << t_tpc_id << " center: ("<< fGeom->TPC(tpcid).GetCenter().X()      << "," << fGeom->TPC(tpcid).GetCenter().Y()      << ","<< fGeom->TPC(tpcid).GetCenter().Z() << ")"
                                 << " box:  ["  << fGeom->TPC(tpcid).BoundingBox().MinX() << "," << fGeom->TPC(tpcid).BoundingBox().MaxX() << "]"
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinY() << "," << fGeom->TPC(tpcid).BoundingBox().MaxY() << "]"
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinZ() << "," << fGeom->TPC(tpcid).BoundingBox().MaxZ() << "]" << std::endl;

      if(fgeoXmin > fGeom->TPC(tpcid).BoundingBox().MinX()) fgeoXmin = fGeom->TPC(tpcid).BoundingBox().MinX();
      if(fgeoXmax < fGeom->TPC(tpcid).BoundingBox().MaxX()) fgeoXmax = fGeom->TPC(tpcid).BoundingBox().MaxX();
      if(fgeoYmin > fGeom->TPC(tpcid).BoundingBox().MinY()) fgeoYmin = fGeom->TPC(tpcid).BoundingBox().MinY();
      if(fgeoYmax < fGeom->TPC(tpcid).BoundingBox().MaxY()) fgeoYmax = fGeom->TPC(tpcid).BoundingBox().MaxY();
      if(fgeoZmin > fGeom->TPC(tpcid).BoundingBox().MinZ()) fgeoZmin = fGeom->TPC(tpcid).BoundingBox().MinZ();
      if(fgeoZmax < fGeom->TPC(tpcid).BoundingBox().MaxZ()) fgeoZmax = fGeom->TPC(tpcid).BoundingBox().MaxZ();
      cout<< "Number of wires : "<<nwires_u<<"  "<<nwires_v <<"  "<<nwires_z<<endl;
    }

  }


  // init summary tree
  art::ServiceHandle<art::TFileService> tfs;

  ftrackTree = tfs->make<TTree>("tracks","Check tracks");
  ftrackTree->Branch("EventNum", &fEventNum, "EventNum/i");
  ftrackTree->Branch("TrackId", &fTrackId, "TrackId/i");
  ftrackTree->Branch("TrackLen",   &ftrackLen,   "TrackLen/F");
  ftrackTree->Branch("Dx",  &ftrackDx,  "Dx/F");
  ftrackTree->Branch("Dy",  &ftrackDy,  "Dy/F");
  ftrackTree->Branch("Dz",  &ftrackDz,  "Dz/F");
  ftrackTree->Branch("StartX",  &ftrackStartX,  "StartX/F");
  ftrackTree->Branch("StartY",  &ftrackStartY,  "StartY/F");
  ftrackTree->Branch("StartZ",  &ftrackStartZ,  "StartZ/F");
  ftrackTree->Branch("StartTick",  &ftrackStartTick,  "StartTick/i");

  ftrackTree->Branch("EndX",  &ftrackEndX,  "EndX/F");
  ftrackTree->Branch("EndY",  &ftrackEndY,  "EndY/F");
  ftrackTree->Branch("EndZ",  &ftrackEndZ,  "EndZ/F");
  ftrackTree->Branch("EndTick",  &ftrackEndTick,  "EndTick/i");

  ftrackTree->Branch("theta",&ftracktheta,  "theta/F");
  ftrackTree->Branch("phi",  &ftrackphi,  "phi/F");
  ftrackTree->Branch("norm", &ftracknorm,  "norm/F");


  fhitTree  = tfs->make<TTree>("hits","Check reconstruction");
  fhitTree->Branch("X", &fhitX, "X/F");
  fhitTree->Branch("Y", &fhitY, "Y/F");
  fhitTree->Branch("Z", &fhitZ, "Z/F");
  fhitTree->Branch("EventNum", &fEventNum, "EventNum/i");
  fhitTree->Branch("TrackId",  &fTrackId, "TrackId/i");
  fhitTree->Branch("Plane",  &fhitPlane, "Plane/i");
  fhitTree->Branch("TPC", &fhitTPC, "TPC/i");
  fhitTree->Branch("Wire",  &fhitWire, "Wire/i");
  fhitTree->Branch("Channel", &fhitChannel, "Channel/i");
  fhitTree->Branch("Time", &fhitTime, "Time/i");

}

//
void pdvdana::CheckHitsAndTracks::endJob()
{
  const string myname = "pdvdana::CheckHitsAndTracks::endJob: ";
  if(fLogLevel>=fFlagBasics){
    cout<<myname<<"tracks processed total     : "<<fTotalTracks<<endl;
    cout<<myname<<"tracks processed selected  : "<<fSelectedTracks<<endl;
  }
}

//
// check the characteristics of a track -> length in x, y, z, total length and theta/phi angles
void pdvdana::CheckHitsAndTracks::checkTrackCharacs( float Dx, float Dy, float Dz, vector<float> &TrackCharac ){

  double norm = sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
  if (Dx>0){
    Dx = -Dx;
    Dz = -Dz;
  }
  TrackCharac[0] = Dx;
  TrackCharac[1] = Dy;
  TrackCharac[2] = Dz;
  TrackCharac[3] = norm;
  TrackCharac[4] = 0;
  TrackCharac[5] = 0;

  if (norm<1)
    return;

  // Determination of the angle
  double theta = 180 - acos(abs(Dx)/norm)*180/3.1415;
  double phi = 0;

  if(Dy>=0){
    if (Dz>=0)
      phi = atan(abs(Dz)/abs(Dy))*180/3.1415 ;
    else{
      phi = 270+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  else{
    if (Dz<0){
      phi= atan(abs(Dz)/abs(Dy))*180/3.1415+180;
    }
    else{
      phi= 90+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  TrackCharac[4] = theta;
  TrackCharac[5] = phi;
  return;
}
void pdvdana::CheckHitsAndTracks::DrawRectangle(TPolyLine *Rectangle, float a_1, float a_2, float b_1, float b_2){
  Rectangle->SetPoint(0,a_1, b_1);
  Rectangle->SetPoint(1,a_1, b_2);
  Rectangle->SetPoint(2,a_2, b_2);
  Rectangle->SetPoint(3,a_2, b_1);
  Rectangle->SetPoint(4,a_1, b_1);
  Rectangle->SetLineColor(12);
  Rectangle->SetLineWidth(3);
  return;
}

void pdvdana::CheckHitsAndTracks::DrawCube(TPolyLine3D *Cube, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max ){
  Cube->SetPoint(0,x_min, y_min, z_min);
  Cube->SetPoint(1,x_max, y_min, z_min);
  Cube->SetPoint(2,x_max, y_max, z_min);
  Cube->SetPoint(3,x_min, y_max, z_min);
  Cube->SetPoint(4,x_min, y_min, z_min);
  Cube->SetPoint(5,x_min, y_min, z_max);
  Cube->SetPoint(6,x_max, y_min, z_max);
  Cube->SetPoint(7,x_max, y_min, z_min);
  Cube->SetPoint(8,x_max, y_max, z_min);
  Cube->SetPoint(9,x_max, y_max, z_max);
  Cube->SetPoint(10,x_min, y_max, z_max);
  Cube->SetPoint(11,x_min, y_min, z_max);
  Cube->SetPoint(12,x_min, y_max, z_max);
  Cube->SetPoint(13,x_min, y_max, z_min);
  Cube->SetPoint(14,x_min, y_max, z_max);
  Cube->SetPoint(15,x_max, y_max, z_max);
  Cube->SetPoint(16,x_max, y_min, z_max);
  Cube->SetLineColor(12);
  Cube->SetLineWidth(3);
  return;
}
void pdvdana::CheckHitsAndTracks::Drawing3D_HitsAndTracks(TCanvas* canvas_3D, vector<vector<float>> l_trackSpacePoints, vector<vector<float>> l_vertexSpacePoints){
  //Drawing 3D tracks

  TH3F *axes = new TH3F("axes", "3D Hits and tracks distribution; X; Y; Z", 1, -400, 400, 1, -400, 400, 1, -10, 400);
  axes->SetDirectory(0);
  canvas_3D->cd();
  axes->Draw();

  TPolyMarker3D* otherPoly = new TPolyMarker3D(l_vertexSpacePoints.size());
  int point=0;
  for (auto const& sp : l_vertexSpacePoints) {
    otherPoly->SetPoint(point, sp[0], sp[1], sp[2]);
    point++;
  }
  otherPoly->SetMarkerStyle(8);
  otherPoly->SetMarkerSize(1);
  otherPoly->SetMarkerColor(kBlue);
  TPolyLine3D* Detector3D = new TPolyLine3D(17);
  DrawCube(Detector3D, fgeoXmin, fgeoYmin, fgeoZmin, fgeoXmax, fgeoYmax, fgeoZmax);

  canvas_3D->cd();
  Detector3D->Draw();
  otherPoly->Draw();

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){
    geo::TPCID tpcid{0, t_tpc_id};
    TPolyLine3D* TPC_3D = new TPolyLine3D(17);
    DrawCube(TPC_3D, fGeom->TPC(tpcid).BoundingBox().MinX(), fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MinZ(),
            fGeom->TPC(tpcid).BoundingBox().MaxX(), fGeom->TPC(tpcid).BoundingBox().MaxY(), fGeom->TPC(tpcid).BoundingBox().MaxZ());
    TPC_3D->SetLineWidth(1);
    TPC_3D->SetLineColor(29);
    canvas_3D->cd();
    TPC_3D->Draw();
  }
  for (auto const& sp : l_trackSpacePoints) {
    TPolyLine3D* t_line = new TPolyLine3D(2);
    t_line->SetPoint(0, sp[0], sp[1], sp[2]);
    t_line->SetPoint(1, sp[3], sp[4], sp[5]);
    t_line->SetLineColor(2);
    t_line->SetLineWidth(2);
    t_line->Draw();
  }
  return;
}
  //Drawing 2D tracks
void pdvdana::CheckHitsAndTracks::Drawing2D_HitsAndTracks(TCanvas* l_canvas_xy, vector<vector<float>> l_trackSpacePoints, vector<vector<float>> l_vertexSpacePoints){

  TPolyMarker *HitsYZ = new TPolyMarker(l_vertexSpacePoints.size());
  TPolyMarker *HitsXY = new TPolyMarker(l_vertexSpacePoints.size());
  TPolyMarker *HitsXZ = new TPolyMarker(l_vertexSpacePoints.size());

  int point=0;

  for (auto const& sp : l_vertexSpacePoints) {
    HitsYZ->SetPoint(point, sp[1], sp[2]);
    HitsXY->SetPoint(point, sp[0], sp[1]);
    HitsXZ->SetPoint(point, sp[0], sp[2]);
    ++point;
  }

  SetHitStyle(HitsYZ);
  SetHitStyle(HitsXY);
  SetHitStyle(HitsXZ);


  TPolyLine *RectangleYZ = new TPolyLine(5);
  TPolyLine *RectangleXY = new TPolyLine(5);
  TPolyLine *RectangleXZ = new TPolyLine(5);

  DrawRectangle(RectangleXZ, fgeoXmin, fgeoXmax, fgeoZmin, fgeoZmax);
  DrawRectangle(RectangleYZ, fgeoYmin, fgeoYmax, fgeoZmin, fgeoZmax);
  DrawRectangle(RectangleXY, fgeoXmin, fgeoXmax, fgeoYmin, fgeoYmax);

  TH2F *axos_XY = new TH2F("axes_1", "Y(X); X; Y", 1, -400, 400, 1, -400, 400);
  TH2F *axos_XZ = new TH2F("axes_2", "Z(X); X; Z", 1, -400, 400, 1, -10, 400);
  TH2F *axos_YZ = new TH2F("axes_3", "Z(Y); Y; Z", 1, -400, 400, 1, -10, 400);
  l_canvas_xy->Divide(2,2);
  l_canvas_xy->cd(1);
  axos_XY->Draw();
  RectangleXY->Draw();
  HitsXY->Draw();
  l_canvas_xy->cd(2);
  axos_XZ->Draw();
  RectangleXZ->Draw();
  HitsXZ->Draw();
  l_canvas_xy->cd(3);
  axos_YZ->Draw();
  RectangleYZ->Draw();
  HitsYZ->Draw();

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){
    geo::TPCID tpcid{0, t_tpc_id};
    TPolyLine* Rectangle_tpc_xy = new TPolyLine(5);
    TPolyLine* Rectangle_tpc_yz = new TPolyLine(5);
    TPolyLine* Rectangle_tpc_xz = new TPolyLine(5);

    DrawRectangle(Rectangle_tpc_xy, fGeom->TPC(tpcid).BoundingBox().MinX(), fGeom->TPC(tpcid).BoundingBox().MaxX(),
      fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MaxY());

    DrawRectangle(Rectangle_tpc_yz, fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MaxY(),
      fGeom->TPC(tpcid).BoundingBox().MinZ(), fGeom->TPC(tpcid).BoundingBox().MaxZ());

    DrawRectangle(Rectangle_tpc_xz, fGeom->TPC(tpcid).BoundingBox().MinX(), fGeom->TPC(tpcid).BoundingBox().MaxX(),
      fGeom->TPC(tpcid).BoundingBox().MinZ(), fGeom->TPC(tpcid).BoundingBox().MaxZ());

    Rectangle_tpc_xy->SetLineWidth(1);
    Rectangle_tpc_yz->SetLineWidth(1);
    Rectangle_tpc_xz->SetLineWidth(1);
    Rectangle_tpc_xy->SetLineColor(29);
    Rectangle_tpc_yz->SetLineColor(29);
    Rectangle_tpc_xz->SetLineColor(29);

    l_canvas_xy->cd(1);
    Rectangle_tpc_xy->Draw();
    l_canvas_xy->cd(3);
    Rectangle_tpc_yz->Draw();
    l_canvas_xy->cd(2);
    Rectangle_tpc_xz->Draw();
  }

  TString vertexName = Form("space_points_%i", ev_num);
  for (auto const& sp : l_trackSpacePoints) {
    TLine* l1 = new TLine(sp[0], sp[1], sp[3], sp[4]);
    TLine* l2 = new TLine(sp[0], sp[2], sp[3], sp[5]);
    TLine* l3 = new TLine(sp[1], sp[2], sp[4], sp[5]);
    l_canvas_xy->cd(1);
    l1->SetLineColor(2);
    l1->SetLineWidth(1);
    l1->Draw();
    l_canvas_xy->cd(2);
    l2->SetLineColor(2);
    l2->SetLineWidth(1);
    l2->Draw();
    l_canvas_xy->cd(3);
    l3->SetLineColor(2);
    l3->SetLineWidth(1);
    l3->Draw();
  }
  return;
}
void pdvdana::CheckHitsAndTracks::DrawingTPC_HitsAndTracks(TCanvas* canvas_tpc, vector<vector<int>> vector_Hits_tpc_channel,
  vector<vector<int>> vector_Hits_tpc_time, vector<vector<int> > fTrackStartWireID,  vector<vector<int> > fTrackEndWireID){
  //Drawing TPC hits and tracks

  vector<TPolyMarker*> EndTracks_tpc(fNplanes);
  vector<TPolyMarker*> StartTracks_tpc(fNplanes);
  vector<TPolyMarker*> Hits_tpc(fNplanes);

  for(unsigned int id_pl =0; id_pl <fNplanes; id_pl++){

    StartTracks_tpc[id_pl] = new TPolyMarker(fTrackStartWireID[id_pl].size());
    EndTracks_tpc[id_pl] = new TPolyMarker(fTrackEndWireID[id_pl].size());
    Hits_tpc[id_pl] = new TPolyMarker(vector_Hits_tpc_channel[id_pl].size());

    //Setting the track polymarker
    for(unsigned int j=0; j<fTrackStartWireID[id_pl].size(); j++){
      StartTracks_tpc[id_pl]->SetPoint(j, fTrackStartWireID[id_pl][j], 2000);
      EndTracks_tpc[id_pl]->SetPoint(j, fTrackEndWireID[id_pl][j], 4000);
    }

    //Setting the hits polymarker
    for(unsigned int j=0; j<vector_Hits_tpc_channel[id_pl].size(); j++){
      Hits_tpc[id_pl]->SetPoint(j, vector_Hits_tpc_channel[id_pl][j], vector_Hits_tpc_time[id_pl][j]);
      if (fLogLevel >= fFlagDetails)
        cout<<"hits : "<<id_pl<<"  "<<j<<"  "<< vector_Hits_tpc_channel[id_pl][j]<<"  "<<vector_Hits_tpc_time[id_pl][j]<<endl;
    }
  }

  TH2F *axos_tpcU = new TH2F("U plane", "U plane; wire; time [ticks]", 1, 0, 4600, 1, 0, 6000);
  TH2F *axos_tpcV = new TH2F("V plane", "V plane; wire; time [ticks]", 1, 0, 4600, 1, 0, 6000);
  TH2F *axos_tpcZ = new TH2F("Z plane", "Z plane; wire; time [ticks]", 1, 0, 4600, 1, 0, 6000);

  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header

  auto legend_U = new TLegend(0.1,0.75,0.42,0.92);
  legend_U->AddEntry(Hits_tpc[0],"Hits","p");
  legend_U->AddEntry(StartTracks_tpc[0],"Start of the track","p");
  legend_U->AddEntry(EndTracks_tpc[0],"End of the track","p");

  auto legend_V = new TLegend(0.1,0.75,0.42,0.92);
  legend_V->AddEntry(Hits_tpc[1],"Hits","p");
  legend_V->AddEntry(StartTracks_tpc[1],"Start of the track","p");
  legend_V->AddEntry(EndTracks_tpc[1],"End of the track","p");

  auto legend_Z = new TLegend(0.1,0.75,0.42,0.92);
  legend_Z->AddEntry(Hits_tpc[2],"Hits","p");
  legend_Z->AddEntry(StartTracks_tpc[2],"Start of the track","p");
  legend_Z->AddEntry(EndTracks_tpc[2],"End of the track","p");

  for (unsigned i=0; i<fNplanes; i++ ){
    SetHitStyle(Hits_tpc[i]);
    StartTracks_tpc[i]->SetMarkerStyle(23);
    StartTracks_tpc[i]->SetMarkerSize(1);
    StartTracks_tpc[i]->SetMarkerColor(kRed);
    EndTracks_tpc[i]->SetMarkerStyle(22);
    EndTracks_tpc[i]->SetMarkerSize(1);
    EndTracks_tpc[i]->SetMarkerColor(kBlack);
  }

  canvas_tpc->Divide(2,2);
  canvas_tpc->cd(1);
  axos_tpcU->Draw();
  EndTracks_tpc[0]->Draw();
  StartTracks_tpc[0]->Draw();
  Hits_tpc[0]->Draw();
  for (unsigned i =0; i < fNtpcs; i++){
    TLine* lines_U = new TLine(fOffsetWireID_u[i],0,fOffsetWireID_u[i],6000);
    canvas_tpc->cd(1);
    lines_U->SetLineColor(29);
    lines_U->SetLineWidth(1);
    lines_U->Draw();
  }
  legend_U->Draw();

  canvas_tpc->cd(2);
  axos_tpcV->Draw();
  EndTracks_tpc[1]->Draw();
  StartTracks_tpc[1]->Draw();
  Hits_tpc[1]->Draw();
  for (unsigned i =0; i < fNtpcs; i++){
    TLine* lines_V = new TLine(fOffsetWireID_v[i],0,fOffsetWireID_v[i],6000);
    canvas_tpc->cd(2);
    lines_V->SetLineColor(29);
    lines_V->SetLineWidth(1);
    lines_V->Draw();
  }
  legend_V->Draw();

  canvas_tpc->cd(3);
  axos_tpcZ->Draw();
  EndTracks_tpc[2]->Draw();
  StartTracks_tpc[2]->Draw();
  Hits_tpc[2]->Draw();
  for (unsigned i =0; i < fNtpcs; i++){
    TLine* lines_Z = new TLine(fOffsetWireID_z[i],0,fOffsetWireID_z[i],6000);
    canvas_tpc->cd(3);
    lines_Z->SetLineColor(29);
    lines_Z->SetLineWidth(1);
    lines_Z->Draw();
  }
  legend_Z->Draw();

  return;
}

void pdvdana::CheckHitsAndTracks::SetHitStyle(TPolyMarker *l_hits){
  l_hits->SetMarkerStyle(8);
  l_hits->SetMarkerSize(0.5);
  l_hits->SetMarkerColor(kBlue);
  return;
}
unsigned pdvdana::CheckHitsAndTracks::GetWireOffset(unsigned plane_id, unsigned tpc_id){

  if (plane_id==0){
    return fOffsetWireID_u[tpc_id];
  }
  else if (plane_id==1){
    return fOffsetWireID_v[tpc_id];
  }
  else if(plane_id==2){
    return fOffsetWireID_z[tpc_id];
  }

  return 0;

}

DEFINE_ART_MODULE(pdvdana::CheckHitsAndTracks)
