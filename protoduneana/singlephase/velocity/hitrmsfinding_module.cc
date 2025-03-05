//////////////////////////////////////////////////////////////////////////////
// Class:       hitrmsfinding                                                   ///
// File:        hitrmsfinding_module.cc                                         ///   
//Description:                                                             /// 
//drift hitrmsfinding calculation module                                        ///
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
  class hitrmsfinding : public art::EDAnalyzer {
  public:
    explicit hitrmsfinding(fhicl::ParameterSet const& pset);
    virtual ~hitrmsfinding();

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

    std::vector<float> trkstartx_crt2;
    std::vector<float> trkendx_crt2;
    std::vector<float> crtreco_x0;
    std::vector<float> crtreco_x1;
    std::vector<float> crtreco_y0;
    std::vector<float> crtreco_y1;
    std::vector<float> crtreco_z0;
    std::vector<float> crtreco_z1;



    std::vector<float> trklen;
    std::vector<int> TrkID; 
    std::vector<float>  xprojectedlen;
    std::vector<double> T0_values;
    //  std::vector<double> t0crt1;
    std::vector<double> t0crt2;
    std::vector<double> crt2tickoffset;
    std::vector<int> tot_trks;
    std::vector< std::vector<float> > hit_peakT0;
    std::vector< std::vector<int> > hit_tpc0;
    std::vector< std::vector<int> > hit_wire0;
    std::vector< std::vector<float> > hit_rms0;
    std::vector< std::vector<float> > hit_deltaT;
    std::vector< std::vector<float> > trkhitx0;
    std::vector< std::vector<float> > trkhity0;
    std::vector< std::vector<float> > trkhitz0;
    std::vector< std::vector<float> > trkdq_amp0;
    std::vector< std::vector<float> >trkdq_int0;
    std::vector< std::vector<float> > hit_peakT1;
    std::vector< std::vector<int> > hit_tpc1;
    std::vector< std::vector<int> > hit_wire1;
    std::vector< std::vector<float> > hit_rms1;
    std::vector< std::vector<float> > trkhitx1;
    std::vector< std::vector<float> > trkhity1;
    std::vector< std::vector<float> > trkhitz1;
    std::vector< std::vector<float> > trkdq_amp1;
    std::vector< std::vector<float> > trkdq_int1;
    std::vector< std::vector<float> > hit_peakT2;
    std::vector< std::vector<int> > hit_tpc2;
    std::vector< std::vector<int> > hit_wire2;
    std::vector< std::vector<float> > hit_rms2;
    std::vector< std::vector<float> > hit_rms_true;
    std::vector< std::vector<float> > trkhitx2;
    std::vector< std::vector<float> > trkhity2;
    std::vector< std::vector<float> > trkhitz2;
    std::vector< std::vector<float> > trkhitz_wire2;
    std::vector< std::vector<float> > trkdq_amp2;
    std::vector< std::vector<float> > trkdq_int2;
    std::vector< std::vector<float> > multiplicity2;
    std::vector< std::vector<float> > goodnessoffit2;
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveTrackInfo;
    bool  fSaveCaloInfo;
  };

  //========================================================================
  hitrmsfinding::hitrmsfinding(fhicl::ParameterSet const& pset) :
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
  hitrmsfinding::~hitrmsfinding(){
  }
  //========================================================================

  //========================================================================
  void hitrmsfinding::beginJob(){
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
    fEventTree->Branch("trkendx_crt2",&trkendx_crt2);
    fEventTree->Branch("trkstartx_crt2",&trkstartx_crt2);
    fEventTree->Branch("crtreco_x0",&crtreco_x0);
    fEventTree->Branch("crtreco_x1",&crtreco_x1);
    fEventTree->Branch("crtreco_y0",&crtreco_y0);
    fEventTree->Branch("crtreco_y1",&crtreco_y1);
    fEventTree->Branch("crtreco_z0",&crtreco_z0);
    fEventTree->Branch("crtreco_z1",&crtreco_z1);
    fEventTree->Branch("crt2tickoffset",&crt2tickoffset);
    fEventTree->Branch("trklen",&trklen);
    fEventTree->Branch("TrkID",&TrkID);
    fEventTree->Branch("tot_trks",&tot_trks);
    fEventTree->Branch("hit_peakT0",&hit_peakT0);
    fEventTree->Branch("hit_tpc0",&hit_tpc0);
    fEventTree->Branch("hit_wire0",&hit_wire0);
    fEventTree->Branch("hit_rms0",&hit_rms0);
    fEventTree->Branch("hit_deltaT",&hit_deltaT);
    fEventTree->Branch("trkhitx0",&trkhitx0);
    fEventTree->Branch("trkhity0",&trkhity0);
    fEventTree->Branch("trkhitz0",&trkhitz0);
    fEventTree->Branch("trkdq_int0",&trkdq_int0);
    fEventTree->Branch("trkdq_amp0",&trkdq_amp0);
    fEventTree->Branch("hit_peakT1",&hit_peakT1);
    fEventTree->Branch("hit_tpc1",&hit_tpc1);
    fEventTree->Branch("hit_wire1",&hit_wire1);
    fEventTree->Branch("hit_rms1",&hit_rms1);
    fEventTree->Branch("trkhitx1",&trkhitx1);
    fEventTree->Branch("trkhity1",&trkhity1);
    fEventTree->Branch("trkhitz1",&trkhitz1);
    fEventTree->Branch("trkdq_int1",&trkdq_int1);
    fEventTree->Branch("trkdq_amp1",&trkdq_amp1);
    fEventTree->Branch("hit_peakT2",&hit_peakT2);
    fEventTree->Branch("hit_tpc2",&hit_tpc2);
    fEventTree->Branch("hit_wire2",&hit_wire2);
    fEventTree->Branch("hit_rms2",&hit_rms2);
    fEventTree->Branch("hit_rms_true",&hit_rms_true);
    fEventTree->Branch("trkhitx2",&trkhitx2);
    fEventTree->Branch("trkhity2",&trkhity2);
    fEventTree->Branch("trkhitz2",&trkhitz2);
    fEventTree->Branch("trkhitz_wire2",&trkhitz_wire2);
    fEventTree->Branch("trkdq_int2",&trkdq_int2);
    fEventTree->Branch("trkdq_amp2",&trkdq_amp2);
    fEventTree->Branch("trkstartcosxyz",&trkstartcosxyz);
    fEventTree->Branch("trkendcosxyz",&trkendcosxyz);
    fEventTree->Branch("multiplicity2",&multiplicity2);
    fEventTree->Branch("goodnessoffit2",&goodnessoffit2);
    // fEventTree->Branch("t0crt1",&t0crt1);
    fEventTree->Branch("t0crt2",&t0crt2);
  }

  //========================================================================
  void hitrmsfinding::endJob(){     

  }

  //========================================================================
  void hitrmsfinding::beginRun(const art::Run&){
    mf::LogInfo("hitrmsfinding")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void hitrmsfinding::analyze( const art::Event& evt){//analyze
    reset();  

    //Detector properties service
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt);
   
    std::vector<art::Ptr<recob::Track> > tracklist;
    auto trackListHandle = evt.getHandle< std::vector<recob::Track> >("pandoraTrack");
    if(trackListHandle) {
      art::fill_ptr_vector(tracklist, trackListHandle);
    }
    else return;

    std::vector<art::Ptr<recob::PFParticle> > pfplist;
    auto PFPListHandle = evt.getHandle< std::vector<recob::PFParticle> >("pandora");
    if(PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);
   
    std::vector<art::Ptr<recob::Hit>> hitlist;
    auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel);
    if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
    art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,"pandoraTrack");
    art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");
    std::vector<const sim::SimChannel*> fSimChannels;
    try{
      evt.getView("largeant", fSimChannels);
    }catch (art::Exception const&e){
    }

    //Get 2-CRT T0
    art::FindManyP<anab::T0> fmt0crt2(trackListHandle, evt, "crtreco");
    art::FindManyP<anab::CosmicTag> fmctcrt2(trackListHandle, evt, "crtreco");
  
    //Get 1-CRT T0
    // art::FindManyP<anab::T0> fmt0crt1(trackListHandle, evt, "crttag");
    // art::FindManyP<anab::CosmicTag> fmctcrt1(trackListHandle, evt, "crttag");

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
    std::vector<float> rms_0;   std::vector<float> rms_1; std::vector<float> rms_2;
    std::vector<float> hitx_0; std::vector<float> hitx_1; std::vector<float> hitx_2;
    std::vector<float> hity_0; std::vector<float> hity_1; std::vector<float> hity_2;
    std::vector<float> hitz_0; std::vector<float> hitz_1; std::vector<float> hitz_2;
    std::vector<float> hitz_wire2; std::vector<float> startcosxyz; std::vector<float> endcosxyz;
    std::vector<float> dT_buffer; std::vector<float> gof, multi,truerms;
   
    float max_value;
    float min_value;
    int trks=0;
    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      double xoffset=0.0;
      //clearing the 1D temporary storage vectors
      peakT_0.clear(); peakT_1.clear();  peakT_2.clear(); 
      tpc_0.clear();  tpc_1.clear();  tpc_2.clear(); 
      wire_0.clear();  wire_1.clear();  wire_2.clear(); 
      int_0.clear(); int_1.clear() ;   int_2.clear(); 
      amp_0.clear(); amp_1.clear() ;   amp_2.clear(); 
      rms_0.clear();   rms_1.clear(); rms_2.clear();
      hitx_0.clear();  hitx_1.clear();  hitx_2.clear();
      hity_0.clear();  hity_1.clear();  hity_2.clear();
      hitz_0.clear();  hitz_1.clear();  hitz_2.clear();
      hitz_wire2.clear(); startcosxyz.clear();endcosxyz.clear();
      dT_buffer.clear(); gof.clear(); multi.clear();
      truerms.clear();


      art::Ptr<recob::Track> ptrack(trackListHandle, i);

      ///this block just saves the t0 values while I include entries with no T0s as well also this saves T0 coming from pandroaTrack alg only
      double t_zero=-999999;
      max_value=0.0;
      min_value=0.0;
      /*      if(fTrackModuleLabel=="pandoraTrack"){ 
	      std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(i);
	      if(!pfps.size()) continue;
	      std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	      // if(!t0s.size()) continue;
	      //auto t0 = t0s.at(0);
	      // double t_zero=t0->Time();
	      if(t0s.size()==0) continue;
	      if(t0s.size()){ 
	      auto t0=t0s.at(0);
	      t_zero=t0->Time();
	      }
	      }
     
	      if(fTrackModuleLabel=="pmtrack"){
	      std::vector<art::Ptr<anab::T0>> T0s=fmT0.at(i);
	      if(T0s.size()==0)
	      continue;
	      }
      */
      // all_trks++;
      const recob::Track& track = *ptrack;


      double this_t0crt2=-DBL_MAX;
      //double this_t0crt1=-DBL_MAX;
      double this_crt2x0 = -DBL_MAX;
      double this_crt2x1 = -DBL_MAX;
      double this_crt2y0 = -DBL_MAX;
      double this_crt2y1 = -DBL_MAX;
      double this_crt2z0 = -DBL_MAX;
      double this_crt2z1 = -DBL_MAX;

      bool test1=true;
      // bool test2=true;
      if(fmt0crt2.isValid()){
	auto const& vt0crt2 = fmt0crt2.at(i);
	if (!vt0crt2.empty()){
	  this_t0crt2 = vt0crt2[0]->Time();
	  test1=false;
	}
	
      }
      if(test1) continue;
      if (fmctcrt2.isValid()){
	auto const& vctcrt2 = fmctcrt2.at(i);
	if (!vctcrt2.empty()){
	  this_crt2x0 = vctcrt2[0]->EndPoint1()[0];
	  this_crt2x1 = vctcrt2[0]->EndPoint2()[0];
	  this_crt2y0 = vctcrt2[0]->EndPoint1()[1];
	  this_crt2y1 = vctcrt2[0]->EndPoint2()[1];
	  this_crt2z0 = vctcrt2[0]->EndPoint1()[2];
	  this_crt2z1 = vctcrt2[0]->EndPoint2()[2];

      }
    }




      /*  if (fmt0crt1.isValid()){
	auto const& vt0crt1 = fmt0crt1.at(i);
	if (!vt0crt1.empty()){
	  this_t0crt1 = vt0crt1[0]->Time();
	  test2=false;
	}
	}*/
      // if(test1 && test2) continue;




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
      auto allHits=fmthm.at(i);
      double ticksoffset=0;

      // if (this_t0crt2 > -DBL_MAX) ticksoffset = this_t0crt2/500.+detProp.GetXTicksOffset(allHits[0]->WireID());
      if (this_t0crt2 > -DBL_MAX) ticksoffset = this_t0crt2/500.+detProp.GetXTicksOffset (allHits[0]->WireID().Plane, allHits[0]->WireID().TPC, allHits[0]->WireID().Cryostat);
      //  else if (this_t0crt1 > -DBL_MAX) ticksoffset = this_t0crt1/500.+detProp.GetXTicksOffset(allHits[0]->WireID());
      xoffset = detProp.ConvertTicksToX(ticksoffset,allHits[0]->WireID());
      std::cout<<"tickoffset , x offset "<<ticksoffset<<"  "<<xoffset<<" default term "<<detProp.GetXTicksOffset (allHits[10]->WireID().Plane, allHits[10]->WireID().TPC, allHits[10]->WireID().Cryostat)<<std::endl;



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
	  planenum=vhit[ii]->WireID().Plane;
	  if(planenum==0){	 
	    peakT_0.push_back(vhit[ii]->PeakTime());
	    tpc_0.push_back(vhit[ii]->WireID().TPC);
	    wire_0.push_back(vhit[ii]->WireID().Wire);
	    int_0.push_back(vhit[ii]->Integral());
	    amp_0.push_back(vhit[ii]->PeakAmplitude());
	    rms_0.push_back(vhit[ii]->RMS());
	    hitx_0.push_back(xpos);
	    hity_0.push_back(ypos);
	    hitz_0.push_back(zpos);		
	  }//planenum 0
	  if(planenum==1){	 
	    peakT_1.push_back(vhit[ii]->PeakTime());
	    tpc_1.push_back(vhit[ii]->WireID().TPC);
	    wire_1.push_back(vhit[ii]->WireID().Wire);
	    int_1.push_back(vhit[ii]->Integral());
	    amp_1.push_back(vhit[ii]->PeakAmplitude());
	    rms_1.push_back(vhit[ii]->RMS());
	    hitx_1.push_back(xpos);
	    hity_1.push_back(ypos);
	    hitz_1.push_back(zpos);	
	  }//planenum 1
	  if(planenum==2){
	    peakT_2.push_back(vhit[ii]->PeakTime()-this_t0crt2/500.0);
	    tpc_2.push_back(vhit[ii]->WireID().TPC);
	    wire_2.push_back(vhit[ii]->WireID().Wire);
	    int_2.push_back(vhit[ii]->Integral());
	    amp_2.push_back(vhit[ii]->PeakAmplitude());
	    rms_2.push_back(vhit[ii]->RMS());
	    dT_buffer.push_back(vhit[ii]->EndTick()-vhit[ii]->StartTick());
	    hitx_2.push_back(xpos);
	    hity_2.push_back(ypos);
	    hitz_2.push_back(zpos);
	    gof.push_back(vhit[ii]->GoodnessOfFit());
	    multi.push_back(vhit[ii]->Multiplicity());
	    double xyzStart[3];
	    double xyzEnd[3];
	    unsigned int wireno=vhit[ii]->WireID().Wire;
            fWireReadoutGeom.WireEndPoints(geo::WireID(0,vhit[ii]->WireID().TPC,2,wireno), xyzStart, xyzEnd);
	    hitz_wire2.push_back(xyzStart[2]);
	    double truermsb=-1;

	    //section for using tuth information
	    unsigned int t0=vhit[ii]->PeakTime()-3*(vhit[ii]->RMS());
	    if(t0<0) t0=0;
	    unsigned int t1=vhit[ii]->PeakTime()+3*(vhit[ii]->RMS());
	    if(t1>6000) t1=6000-1;
	    
	    if(!evt.isRealData()){
	      const sim::SimChannel* chan=0;
	      for(size_t sc=0;sc<fSimChannels.size();++sc){
		if(fSimChannels[sc]->Channel()==vhit[ii]->Channel()) chan=fSimChannels[sc];
	      }
	      if(chan){
		const auto &tdcidemap = chan->TDCIDEMap();
		double hit_ch=0;
		double mean_t=0;
		double mean_t2=0;
		double hit_rmsvalue=0;
		for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
		  if(!((*mapitr).first>t0 && (*mapitr).first<=t1)) continue;
		  // loop over the vector of IDE objects.
		  int hit_nelec=0;
		  const std::vector<sim::IDE> idevec = (*mapitr).second;
		  for(size_t iv = 0; iv < idevec.size(); ++iv){ 
		    hit_nelec+=idevec[iv].numElectrons;
		  }//iv loop
		  hit_ch+=hit_nelec;
		  // std::cout<<"hit_ch "<<hit_ch<<std::endl; 
		  int j=(*mapitr).first;
		  mean_t+=double(j)*double(hit_nelec);
		  double jterm=double(j)*double(j)*double(hit_nelec);
		  mean_t2=mean_t2+jterm;
		}//mapitr loop
		mean_t/=hit_ch;
		mean_t2/=hit_ch;
		hit_rmsvalue=sqrt(mean_t2-mean_t*mean_t);
		//	std::cout<<"hit rms "<<hit_rmsvalue<<" reco rms "<<vhit[ii]->RMS()<<std::endl;
		truermsb=hit_rmsvalue;
	      }//if(chan)
	    }//!evt.IsRealData
	    truerms.push_back(truermsb);

	   
	  }//planenum 2
	}//loop over vhit
      }//fmthm valid
      //hits and calorimetry loop
      if(peakT_2.size()<10) continue;
      max_value=*std::max_element(peakT_2.begin(),peakT_2.end());
      min_value=*std::min_element(peakT_2.begin(),peakT_2.end());
      // if(max_value-min_value<4300) continue;
      std::cout<<max_value<<"  "<<min_value<<std::endl;
      trks++;
      hit_peakT0.push_back(peakT_0);
      hit_tpc0.push_back(tpc_0);
      hit_wire0.push_back(wire_0);
      hit_rms0.push_back(rms_0);
      trkhitx0.push_back(hitx_0);
      trkhity0.push_back(hity_0);
      trkhitz0.push_back(hitz_0);
      trkdq_amp0.push_back(amp_0);
      trkdq_int0.push_back(int_0);

      hit_peakT1.push_back(peakT_1);
      hit_tpc1.push_back(tpc_1);
      hit_wire1.push_back(wire_1);
      hit_rms1.push_back(rms_1);
      trkhitx1.push_back(hitx_1);
      trkhity1.push_back(hity_1);
      trkhitz1.push_back(hitz_1);
      trkdq_amp1.push_back(amp_1);
      trkdq_int1.push_back(int_1);

      hit_peakT2.push_back(peakT_2);
      hit_tpc2.push_back(tpc_2);
      hit_wire2.push_back(wire_2);
      hit_rms2.push_back(rms_2);
      hit_rms_true.push_back(truerms);
      trkhitx2.push_back(hitx_2);
      trkhity2.push_back(hity_2);
      trkhitz2.push_back(hitz_2);
      trkhitz_wire2.push_back(hitz_wire2);
      trkdq_amp2.push_back(amp_2);
      trkdq_int2.push_back(int_2);
      hit_deltaT.push_back(dT_buffer);
      multiplicity2.push_back(multi);
      goodnessoffit2.push_back(gof);

     

      xprojectedlen.push_back(TMath::Abs(end.X()-pos.X()));
      trackthetaxz.push_back(theta_xz);
      trackthetayz.push_back(theta_yz);
      trkstartx.push_back(pos.X());
      trkstartx_crt2.push_back(pos.X()-xoffset);
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
      trkendx_crt2.push_back(end.X()-xoffset);
      crtreco_x0.push_back(this_crt2x0);
      crtreco_x1.push_back(this_crt2x1);
      crtreco_y0.push_back(this_crt2y0);
      crtreco_y1.push_back(this_crt2y1);
      crtreco_z0.push_back(this_crt2z0);
      crtreco_z1.push_back(this_crt2z1);

      trkendy.push_back(end.Y());
      trkendz.push_back(end.Z());
      trklen.push_back(track.Length());
      TrkID.push_back(track.ID());
      T0_values.push_back(t_zero);
      crt2tickoffset.push_back(ticksoffset);
      t0crt2.push_back(this_t0crt2/500.0);
    } // loop over trks...
    tot_trks.push_back(trks);
    fEventTree->Fill();
  } // end of analyze function
	   
  /////////////////// Defintion of reset function ///////////
  void hitrmsfinding::reset(){
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
    trkstartx_crt2.clear();
    trkendx_crt2.clear();
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
    hit_rms0.clear();
    hit_peakT1.clear();
    hit_tpc1.clear();
    hit_wire1.clear();
    trkhitx1.clear();
    trkhity1.clear();
    trkhitz1.clear();
    trkdq_int1.clear();
    trkdq_amp1.clear();
    hit_rms1.clear();
    hit_peakT2.clear();
    hit_tpc2.clear();
    hit_wire2.clear();
    trkhitx2.clear();
    trkhity2.clear();
    trkhitz2.clear();
    trkhitz_wire2.clear();
    trkdq_int2.clear();
    trkdq_amp2.clear();
    hit_rms2.clear();
    hit_deltaT.clear();
    multiplicity2.clear();
    goodnessoffit2.clear();
    hit_rms_true.clear();
    crtreco_x0.clear();
    crtreco_x1.clear();
    crtreco_y0.clear();
    crtreco_y1.clear();
    crtreco_z0.clear();
    crtreco_z1.clear();
    crt2tickoffset.clear();

    //  t0crt1.clear();
    t0crt2.clear();

  }
  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(hitrmsfinding)
}
