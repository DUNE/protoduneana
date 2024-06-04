//////////////////////////////////////////////////////////////////////////////
// Class:       velocity                                                   ///
// File:        velocity_module.cc                                         ///   
//Description:                                                             /// 
//drift velocity calculation module                                        ///
//dumps all the infomrmation in TTrees which needs to analysed             ///
//contact person:apaudel@phys.ksu.edu                                      ///
//////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "larcore/Geometry/WireReadout.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

//Root and C++ include
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

//using namespace std;


namespace protoana{
  class velocity : public art::EDAnalyzer {
  public:
    explicit velocity(fhicl::ParameterSet const& pset);
    virtual ~velocity();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
  private:
    ProtoDUNEDataUtils fDataUtils;
    TTree* fEventTree;
    geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

    //These are the tree variables I will be using
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Double_t evttime; 
    int fNactivefembs[6];
    std::vector<float> trackthetaxz;
    std::vector<float>  trackthetayz;
    std::vector<float> trkstartx;
    std::vector<float> trkstarty;
    std::vector<float> trkstartz;
    std::vector< std::vector<float> > trkstartcosxyz;
    std::vector< std::vector<float> > trkendcosxyz;
    std::vector<float> trkendx;
    std::vector<float> trkendy;
    std::vector<float> trkendz;
    std::vector<float> trklen;
    std::vector<int> TrkID; 
    std::vector<float>  xprojectedlen;
    std::vector<double> T0_values;
    std::vector<int> tot_trks;
    std::vector< std::vector<float> > hit_peakT0;
    std::vector< std::vector<int> > hit_tpc0;
    std::vector< std::vector<int> > hit_wire0;
    std::vector< std::vector<int> > hit_channel0;
    std::vector< std::vector<float> > trkhitx0;
    std::vector< std::vector<float> > trkhity0;
    std::vector< std::vector<float> > trkhitz0;
    std::vector< std::vector<float> > trkdq_amp0;
    std::vector< std::vector<float> >trkdq_int0;
    std::vector< std::vector<float> > hit_peakT1;
    std::vector< std::vector<int> > hit_tpc1;
    std::vector< std::vector<int> > hit_wire1;
    std::vector< std::vector<int> > hit_channel1;
    std::vector< std::vector<float> > trkhitx1;
    std::vector< std::vector<float> > trkhity1;
    std::vector< std::vector<float> > trkhitz1;
    std::vector< std::vector<float> > trkdq_amp1;
    std::vector< std::vector<float> > trkdq_int1;
    std::vector< std::vector<float> > hit_peakT2;
    std::vector< std::vector<int> > hit_tpc2;
    std::vector< std::vector<int> > hit_wire2;
    std::vector< std::vector<int> > hit_channel2;
    std::vector< std::vector<float> > trkhitx2;
    std::vector< std::vector<float> > trkhity2;
    std::vector< std::vector<float> > trkhitz2;
    std::vector< std::vector<float> > trkhitz_wire2;
    std::vector< std::vector<float> > trkdq_amp2;
    std::vector< std::vector<float> >trkdq_int2;
    std::vector< std::vector<float> > trkhittruex0;
    std::vector< std::vector<float> > trkhittruey0;
    std::vector< std::vector<float> > trkhittruez0;
    std::vector< std::vector<float> > trkhittruex1;
    std::vector< std::vector<float> > trkhittruey1;
    std::vector< std::vector<float> > trkhittruez1;
    std::vector< std::vector<float> > trkhittruex2;
    std::vector< std::vector<float> > trkhittruey2;
    std::vector< std::vector<float> > trkhittruez2;
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveTrackInfo;
    bool  fSaveCaloInfo;
  };

  //========================================================================
  velocity::velocity(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fDataUtils                  (pset.get<fhicl::ParameterSet>("DataUtils")),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ), 
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false))

  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
 
  //========================================================================
  velocity::~velocity(){
  }
  //========================================================================

  //========================================================================
  void velocity::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("Nactivefembs",&fNactivefembs,"Nactivefembs[6]/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("T0_values",&T0_values);
    fEventTree->Branch("xprojectedlen",&xprojectedlen);
    fEventTree->Branch("trackthetaxz",&trackthetaxz);
    fEventTree->Branch("trackthetayz",&trackthetayz);
    fEventTree->Branch("trkstartx",&trkstartx);
    fEventTree->Branch("trkstarty",&trkstarty);
    fEventTree->Branch("trkstartz",&trkstartz);
    fEventTree->Branch("trkendx",&trkendx);
    fEventTree->Branch("trkendy",&trkendy);
    fEventTree->Branch("trkendz",&trkendz);
    fEventTree->Branch("trklen",&trklen);
    fEventTree->Branch("TrkID",&TrkID);
    fEventTree->Branch("tot_trks",&tot_trks);
    fEventTree->Branch("hit_peakT0",&hit_peakT0);
    fEventTree->Branch("hit_tpc0",&hit_tpc0);
    fEventTree->Branch("hit_wire0",&hit_wire0);
    fEventTree->Branch("hit_channel0",&hit_channel0);
    fEventTree->Branch("trkhitx0",&trkhitx0);
    fEventTree->Branch("trkhity0",&trkhity0);
    fEventTree->Branch("trkhitz0",&trkhitz0);
    fEventTree->Branch("trkdq_int0",&trkdq_int0);
    fEventTree->Branch("trkdq_amp0",&trkdq_amp0);
    fEventTree->Branch("hit_peakT1",&hit_peakT1);
    fEventTree->Branch("hit_tpc1",&hit_tpc1);
    fEventTree->Branch("hit_wire1",&hit_wire1);
    fEventTree->Branch("hit_channel1",&hit_channel1);
    fEventTree->Branch("trkhitx1",&trkhitx1);
    fEventTree->Branch("trkhity1",&trkhity1);
    fEventTree->Branch("trkhitz1",&trkhitz1);
    fEventTree->Branch("trkdq_int1",&trkdq_int1);
    fEventTree->Branch("trkdq_amp1",&trkdq_amp1);
    fEventTree->Branch("hit_peakT2",&hit_peakT2);
    fEventTree->Branch("hit_tpc2",&hit_tpc2);
    fEventTree->Branch("hit_wire2",&hit_wire2);
    fEventTree->Branch("hit_channel2",&hit_channel2);
    fEventTree->Branch("trkhitx2",&trkhitx2);
    fEventTree->Branch("trkhity2",&trkhity2);
    fEventTree->Branch("trkhitz2",&trkhitz2);
    fEventTree->Branch("trkhitz_wire2",&trkhitz_wire2);
    fEventTree->Branch("trkdq_int2",&trkdq_int2);
    fEventTree->Branch("trkdq_amp2",&trkdq_amp2);
    fEventTree->Branch("trkstartcosxyz",&trkstartcosxyz);
    fEventTree->Branch("trkendcosxyz",&trkendcosxyz);
    fEventTree->Branch("trkhittruex0",&trkhittruex0);
    fEventTree->Branch("trkhittruey0",&trkhittruey0);
    fEventTree->Branch("trkhittruez0",&trkhittruez0);
    fEventTree->Branch("trkhittruex1",&trkhittruex1);
    fEventTree->Branch("trkhittruey1",&trkhittruey1);
    fEventTree->Branch("trkhittruez1",&trkhittruez1);
    fEventTree->Branch("trkhittruex2",&trkhittruex2);
    fEventTree->Branch("trkhittruey2",&trkhittruey2);
    fEventTree->Branch("trkhittruez2",&trkhittruez2);

  }

  //========================================================================
  void velocity::endJob(){     

  }

  //========================================================================
  void velocity::beginRun(const art::Run&){
    mf::LogInfo("velocity")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void velocity::analyze( const art::Event& evt){//analyze
    reset();  

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
   
    std::vector<art::Ptr<recob::Track> > tracklist;
    auto trackListHandle = evt.getHandle< std::vector<recob::Track> >("pandoraTrack");
    if(trackListHandle){
      art::fill_ptr_vector(tracklist, trackListHandle);
    }
    else return;


    std::vector<art::Ptr<recob::PFParticle> > pfplist;
    auto PFPListHandle = evt.getHandle< std::vector<recob::PFParticle> >("pandora");
    if(PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);
  
   
      // get information about the hits
    std::vector<art::Ptr<recob::Hit>> hitlist;
    auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel); 
    if(hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
    art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,"pandoraTrack");
    art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");


    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime=tts.AsDouble();



    // Get number of active fembs
    if(!evt.isRealData()){
      for(int k=0; k < 6; k++)
	fNactivefembs[k] = 20;
    }
    else{
      for(int k=0; k < 6; k++)
	fNactivefembs[k] = fDataUtils.GetNActiveFembsForAPA(evt, k);
    }


    //defining the 1D temporary storage vectors
    std::vector< float> peakT_0;  std::vector< float> peakT_1;  std::vector< float> peakT_2; 
    std::vector< int > tpc_0;  std::vector< int > tpc_1;  std::vector<int> tpc_2; 
    std::vector< int > wire_0;  std::vector< int > wire_1;  std::vector<int> wire_2; 
    std::vector< float > int_0;  std::vector< float > int_1;  std::vector<float> int_2; 
    std::vector< float > amp_0;  std::vector< float > amp_1;  std::vector<float> amp_2; 
    std::vector<int> channel_0;   std::vector<int> channel_1; std::vector<int> channel_2;
    std::vector<float> hitx_0; std::vector<float> hitx_1; std::vector<float> hitx_2;
    std::vector<float> hity_0; std::vector<float> hity_1; std::vector<float> hity_2;
    std::vector<float> hitz_0; std::vector<float> hitz_1; std::vector<float> hitz_2;
    std::vector<float> hittruex_0; std::vector<float> hittruex_1; std::vector<float> hittruex_2;
    std::vector<float> hittruey_0; std::vector<float> hittruey_1; std::vector<float> hittruey_2;
    std::vector<float> hittruez_0; std::vector<float> hittruez_1; std::vector<float> hittruez_2;
    std::vector<float> hitz_wire2; std::vector<float> startcosxyz; std::vector<float> endcosxyz;
   
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    float max_value;
    float min_value;
    int trks=0;
    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      
      //clearing the 1D temporary storage vectors
      peakT_0.clear(); peakT_1.clear();  peakT_2.clear(); 
      tpc_0.clear();  tpc_1.clear();  tpc_2.clear(); 
      wire_0.clear();  wire_1.clear();  wire_2.clear(); 
      int_0.clear(); int_1.clear() ;   int_2.clear(); 
      amp_0.clear(); amp_1.clear() ;   amp_2.clear(); 
      channel_0.clear();   channel_1.clear(); channel_2.clear();
      hitx_0.clear();  hitx_1.clear();  hitx_2.clear();
      hity_0.clear();  hity_1.clear();  hity_2.clear();
      hitz_0.clear();  hitz_1.clear();  hitz_2.clear();
      hittruex_0.clear(); hittruey_0.clear(); hittruez_0.clear(); 
      hittruex_1.clear(); hittruey_1.clear(); hittruez_1.clear(); 
      hittruex_2.clear(); hittruey_2.clear(); hittruez_2.clear(); 
      hitz_wire2.clear(); startcosxyz.clear();endcosxyz.clear();



      art::Ptr<recob::Track> ptrack(trackListHandle, i);

      ///this block just saves the t0 values while I include entries with no T0s as well also this saves T0 coming from pandroaTrack alg only
      double t_zero=-999999;
      max_value=0.0;
      min_value=0.0;
      if(fTrackModuleLabel=="pandoraTrack"){ 
	std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(i);
	if(!pfps.size()) continue;
	std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	// if(!t0s.size()) continue;
	//auto t0 = t0s.at(0);
	// double t_zero=t0->Time();
	if(t0s.size()){ 
	  auto t0=t0s.at(0);
	  t_zero=t0->Time();
	}
      }
     
      /* 
	 if(fTrackModuleLabel=="pmtrack"){
	 std::vector<art::Ptr<anab::T0>> T0s=fmT0.at(i);
	 if(T0s.size()==0)
	 continue;
	 }
      */
      // all_trks++;
      const recob::Track& track = *ptrack;
      auto pos = track.Vertex();
      auto dir_start = track.VertexDirection();
      auto dir_end   = track.EndDirection();
      auto end = track.End();
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
     
      int planenum=999;
      float xpos=-9999;
      float ypos=-9999;
      float zpos=-9999;

      //hits and calorimetry loop
      if(fmthm.isValid()){
	auto vhit=fmthm.at(i);
	auto vmeta=fmthm.data(i);
	for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
	  bool fBadhit = false;
	  if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
	    fBadhit = true;
	    //cout<<"fBadHit"<<fBadhit<<endl;
	    continue;
	  }
	  if (vmeta[ii]->Index()>=tracklist[i]->NumberTrajectoryPoints()){
	    throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[i]->NumberTrajectoryPoints()<<" for track index "<<i<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
	  }
	  if (!tracklist[i]->HasValidPoint(vmeta[ii]->Index())){
	    fBadhit = true;
	    // cout<<"had valid point "<<fBadhit<<endl;
	    continue;
	  }
	
	  auto loc = tracklist[i]->LocationAtPoint(vmeta[ii]->Index());
	  xpos=loc.X();
	  ypos=loc.Y();
	  zpos=loc.Z();
	  //	cout<<"x, y, z "<<xpos<<"  "<<ypos<<"  "<<zpos<<endl;
	  //	cout<<"BadHit"<<fBadhit<<endl;
	  if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
	  if (zpos<-100) continue; //hit not on track

          float truexpos = -9999;
          float trueypos = -9999;
          float truezpos = -9999;
          
          if (!evt.isRealData()){
            auto simides = bt_serv->HitToSimIDEs_Ps(clockData, vhit[ii]);
            std::map<int, float> idemap;
            for (auto const & ide: simides){
              idemap[ide->trackID]+=ide->energy;
            }
            int trackid = 0;
            float maxe = -1;
            for (auto ii = idemap.begin(); ii!=idemap.end(); ++ii){
              if ((ii->second)>maxe){
                maxe = ii->second;
                trackid = ii->first;
              }
            }
            if (trackid){
              float xavg = 0;
              float yavg = 0;
              float zavg = 0;
              float tote = 0;
              for (auto const & ide: simides){
                if (ide->trackID != trackid) continue;
                tote += ide->energy;
                xavg += ide->x*ide->energy;
                yavg += ide->y*ide->energy;
                zavg += ide->z*ide->energy;
              }
              if (tote){
                truexpos = xavg/tote;
                trueypos = yavg/tote;
                truezpos = zavg/tote;
              }
            }
          }

	  planenum=vhit[ii]->WireID().Plane;
	  if(planenum==0){	 
	    peakT_0.push_back(vhit[ii]->PeakTime());
	    tpc_0.push_back(vhit[ii]->WireID().TPC);
	    wire_0.push_back(vhit[ii]->WireID().Wire);
	    int_0.push_back(vhit[ii]->Integral());
	    amp_0.push_back(vhit[ii]->PeakAmplitude());
	    channel_0.push_back(vhit[ii]->Channel());
	    hitx_0.push_back(xpos);
	    hity_0.push_back(ypos);
	    hitz_0.push_back(zpos);		
            hittruex_0.push_back(truexpos);
            hittruey_0.push_back(trueypos);
            hittruez_0.push_back(truezpos);
	  }//planenum 0
	  if(planenum==1){	 
	    peakT_1.push_back(vhit[ii]->PeakTime());
	    tpc_1.push_back(vhit[ii]->WireID().TPC);
	    wire_1.push_back(vhit[ii]->WireID().Wire);
	    int_1.push_back(vhit[ii]->Integral());
	    amp_1.push_back(vhit[ii]->PeakAmplitude());
	    channel_1.push_back(vhit[ii]->Channel());
	    hitx_1.push_back(xpos);
	    hity_1.push_back(ypos);
	    hitz_1.push_back(zpos);	
            hittruex_1.push_back(truexpos);
            hittruey_1.push_back(trueypos);
            hittruez_1.push_back(truezpos);
	  }//planenum 1
	  if(planenum==2){
	    peakT_2.push_back(vhit[ii]->PeakTime());
	    tpc_2.push_back(vhit[ii]->WireID().TPC);
	    wire_2.push_back(vhit[ii]->WireID().Wire);
	    int_2.push_back(vhit[ii]->Integral());
	    amp_2.push_back(vhit[ii]->PeakAmplitude());
	    channel_2.push_back(vhit[ii]->Channel());
	    hitx_2.push_back(xpos);
	    hity_2.push_back(ypos);
	    hitz_2.push_back(zpos);
            hittruex_2.push_back(truexpos);
            hittruey_2.push_back(trueypos);
            hittruez_2.push_back(truezpos);
	    double xyzStart[3];
	    double xyzEnd[3];
	    unsigned int wireno=vhit[ii]->WireID().Wire;
            fWireReadoutGeom.WireEndPoints(geo::WireID(0,vhit[ii]->WireID().TPC,2,wireno), xyzStart, xyzEnd);
	    hitz_wire2.push_back(xyzStart[2]);
            //std::cout<<xpos<<" "<<ypos<<" "<<zpos<<" "<<truexpos<<" "<<trueypos<<" "<<truezpos<<std::endl;
	  }//planenum 2
	}//loop over vhit
      }//fmthm valid
      //hits and calorimetry loop
      if(peakT_2.size()<20) continue;
      max_value=*std::max_element(peakT_2.begin(),peakT_2.end());
      min_value=*std::min_element(peakT_2.begin(),peakT_2.end());
      if(max_value-min_value<4300) continue;
      std::cout<<max_value<<"  "<<min_value<<std::endl;
      trks++;
      hit_peakT0.push_back(peakT_0);
      hit_tpc0.push_back(tpc_0);
      hit_wire0.push_back(wire_0);
      hit_channel0.push_back(channel_0);
      trkhitx0.push_back(hitx_0);
      trkhity0.push_back(hity_0);
      trkhitz0.push_back(hitz_0);
      trkhittruex0.push_back(hittruex_0);
      trkhittruey0.push_back(hittruey_0);
      trkhittruez0.push_back(hittruez_0);
      trkdq_amp0.push_back(amp_0);
      trkdq_int0.push_back(int_0);

      hit_peakT1.push_back(peakT_1);
      hit_tpc1.push_back(tpc_1);
      hit_wire1.push_back(wire_1);
      hit_channel1.push_back(channel_1);
      trkhitx1.push_back(hitx_1);
      trkhity1.push_back(hity_1);
      trkhitz1.push_back(hitz_1);
      trkhittruex1.push_back(hittruex_1);
      trkhittruey1.push_back(hittruey_1);
      trkhittruez1.push_back(hittruez_1);
      trkdq_amp1.push_back(amp_1);
      trkdq_int1.push_back(int_1);

      hit_peakT2.push_back(peakT_2);
      hit_tpc2.push_back(tpc_2);
      hit_wire2.push_back(wire_2);
      hit_channel2.push_back(channel_2);
      trkhitx2.push_back(hitx_2);
      trkhity2.push_back(hity_2);
      trkhitz2.push_back(hitz_2);
      trkhittruex2.push_back(hittruex_2);
      trkhittruey2.push_back(hittruey_2);
      trkhittruez2.push_back(hittruez_2);
      trkhitz_wire2.push_back(hitz_wire2);
      trkdq_amp2.push_back(amp_2);
      trkdq_int2.push_back(int_2);

      xprojectedlen.push_back(TMath::Abs(end.X()-pos.X()));
      trackthetaxz.push_back(theta_xz);
      trackthetayz.push_back(theta_yz);
      trkstartx.push_back(pos.X());
      trkstarty.push_back(pos.Y());
      trkstartz.push_back(pos.Z());
      startcosxyz.push_back(dir_start.X());
      startcosxyz.push_back(dir_start.Y());
      startcosxyz.push_back(dir_start.Z());
      endcosxyz.push_back(dir_end.X());
      endcosxyz.push_back(dir_end.Y());
      endcosxyz.push_back(dir_end.Z());

      trkstartcosxyz.push_back(startcosxyz);
      trkendcosxyz.push_back(endcosxyz);
      trkendx.push_back(end.X());
      trkendy.push_back(end.Y());
      trkendz.push_back(end.Z());
      trklen.push_back(track.Length());
      TrkID.push_back(track.ID());
      T0_values.push_back(t_zero);
    } // loop over trks...
    tot_trks.push_back(trks);
    fEventTree->Fill();
  } // end of analyze function
	   
  /////////////////// Defintion of reset function ///////////
  void velocity::reset(){
    run = -9999;
    subrun = -9999;
    event = -9999;
    evttime = -9999;
    //all_trks = -9999;
    for(int k=0; k < 6; k++)
      fNactivefembs[k] = -9999;
    trackthetaxz.clear();
    trackthetayz.clear();
    trkstartx.clear();
    trkstarty.clear();
    trkstartz.clear();
    trkstartcosxyz.clear();
    trkendcosxyz.clear();
    trkendx.clear();
    trkendy.clear();
    trkendz.clear();
    trklen.clear();
    TrkID.clear();
    tot_trks.clear();
    T0_values.clear();
    xprojectedlen.clear();
    hit_peakT0.clear();
    hit_tpc0.clear();
    hit_wire0.clear();
    trkhitx0.clear();
    trkhity0.clear();
    trkhitz0.clear();
    trkdq_int0.clear();
    trkdq_amp0.clear();
    hit_channel0.clear();
    hit_peakT1.clear();
    hit_tpc1.clear();
    hit_wire1.clear();
    trkhitx1.clear();
    trkhity1.clear();
    trkhitz1.clear();
    trkdq_int1.clear();
    trkdq_amp1.clear();
    hit_channel1.clear();
    hit_peakT2.clear();
    hit_tpc2.clear();
    hit_wire2.clear();
    trkhitx2.clear();
    trkhity2.clear();
    trkhitz2.clear();
    trkhitz_wire2.clear();
    trkdq_int2.clear();
    trkdq_amp2.clear();
    hit_channel2.clear();
    trkhittruex0.clear();
    trkhittruey0.clear();
    trkhittruez0.clear();
    trkhittruex1.clear();
    trkhittruey1.clear();
    trkhittruez1.clear();
    trkhittruex2.clear();
    trkhittruey2.clear();
    trkhittruez2.clear();

  }
  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(velocity)
}
