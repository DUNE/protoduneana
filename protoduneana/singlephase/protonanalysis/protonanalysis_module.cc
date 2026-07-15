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

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
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
#include "TH3.h"
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

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

const int kMaxTracks  = 200;
//const int nplots =14;

using namespace std;

namespace dune{

  class protonanalysis : public art::EDAnalyzer {
  public:

    explicit protonanalysis(fhicl::ParameterSet const& pset);
    virtual ~protonanalysis();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
  private:
    TH1F* pdg_beamparticles;
    TH1F* calibration_constants;
    TH2F*  electric_field_YZ;
    TH2F*  electric_field_YZ_3;
    TH1F* EFieldX;
    TH3F* EField3DX;
    TTree* fEventTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
   
    Double_t Total_energy[kMaxTracks];
    Double_t EndE[kMaxTracks];
    Int_t pdg[kMaxTracks];
    Double_t evttime; 
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    beam_proton;
    Int_t tot_trks;
    Int_t total_trks;
    Int_t pdg_all_beam[kMaxTracks];
    Float_t  xprojectedlen[kMaxTracks];
    Float_t  trackthetaxz[kMaxTracks];
    Float_t  trackthetayz[kMaxTracks];
    Float_t trkstartx[kMaxTracks];
    Float_t trkstarty[kMaxTracks];
    Float_t trkstartz[kMaxTracks];
    Float_t trkendx[kMaxTracks];
    Float_t trkendy[kMaxTracks];
    Float_t trkendz[kMaxTracks];
    Float_t trklen[kMaxTracks];
    Int_t    TrkID[kMaxTracks]; 
    Int_t    TrackId_geant[kMaxTracks]; 
    Int_t origin[kMaxTracks];
    Float_t EndPointx[kMaxTracks];
    Float_t EndPointy[kMaxTracks];
    Float_t EndPointz[kMaxTracks];
    Float_t EndPointx_corrected[kMaxTracks];
    Float_t EndPointy_corrected[kMaxTracks];
    Float_t EndPointz_corrected[kMaxTracks];
    Float_t StartPointx[kMaxTracks];
    Float_t StartPointy[kMaxTracks];
    Float_t StartPointz[kMaxTracks];
    Float_t  trkstartcosxyz[kMaxTracks][3];
    Float_t  trkendcosxyz[kMaxTracks][3];
    Int_t    ntrkhits[kMaxTracks][3];
    Float_t  trkdqdx[kMaxTracks][3][3000];
    Float_t  trkdedx[kMaxTracks][3][3000];
    Float_t  trkresrange[kMaxTracks][3][3000];
    Float_t  trkhitx[kMaxTracks][3][3000];
    Float_t  trkhity[kMaxTracks][3][3000];
    Float_t  trkhitz[kMaxTracks][3][3000];
    Float_t  trkpitch[kMaxTracks][3][3000];
    Float_t  Efield[kMaxTracks][3][3000];
    Float_t peakT_max[kMaxTracks];
    Float_t peakT_min[kMaxTracks];
    Double_t momentum[kMaxTracks];

    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
  }; 

  //========================================================================
  protonanalysis::protonanalysis(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false))
  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
 
  //========================================================================
  protonanalysis::~protonanalysis(){
  }
  //========================================================================

  //========================================================================
  void protonanalysis::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    pdg_beamparticles=tfs->make<TH1F>("beam_pdg","beam particle pdgs;pdg;no of beam particles",6000,-2999.5,3000.5);
    calibration_constants=tfs->make<TH1F>("calibration_constants","calibration constants;calibration constants;no of hits",6000,0.003000,0.009000);
    electric_field_YZ=tfs->make<TH2F>("EField_YZ_x_0","Electric field for YZ plane at  x= -300cm;Z coordinate;Y coordinate",139,0,695,120,0,600);
    electric_field_YZ_3=tfs->make<TH2F>("EField_YZ_x_3","Electric field for YZ plane at  x= -25cm;Z coordinate;Y coordinate",139,0,695,120,0,600);
    EFieldX=tfs->make<TH1F>("EFieldX","Electric Field for Y=300cm && Z=350cm;X coordinate;Efield in kV/cm",120,-360,360);
    EField3DX=tfs->make<TH3F>("EField3DX","Electric FieldX 3D",144,-360,360,120,0,600,139,0,695);

    pdg_beamparticles->Sumw2();
    calibration_constants->Sumw2();
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event,"event/I");
  

    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
    fEventTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    fEventTree->Branch("beam_proton",&beam_proton,"beam_proton/I");//total tracks pasing broken tracks filter
    fEventTree->Branch("tot_trks",&tot_trks,"tot_trks/I");
   
    fEventTree->Branch("total_trks",&total_trks,"total_trks/I");//total proton beam tracks
    fEventTree->Branch("xprojectedlen",xprojectedlen,"xprojectedlen[beam_proton]/F");
    fEventTree->Branch("trackthetaxz",trackthetaxz,"trackthetaxz[beam_proton]/F");
    fEventTree->Branch("trackthetayz",trackthetayz,"trackthetayz[beam_proton]/F");
    fEventTree->Branch("trkstartx",trkstartx,"trkstartx[beam_proton]/F");
    fEventTree->Branch("trkstarty",trkstarty,"trkstarty[beam_proton]/F");
    fEventTree->Branch("trkstartz",trkstartz,"trkstartz[beam_proton]/F");
    fEventTree->Branch("trkendx",trkendx,"trkendx[beam_proton]/F");
    fEventTree->Branch("trkendy",trkendy,"trkendy[beam_proton]/F");
    fEventTree->Branch("trkendz",trkendz,"trkendz[beam_proton]/F");
    fEventTree->Branch("trklen",trklen,"trklen[beam_proton]/F");
    fEventTree->Branch("Total_energy",Total_energy,"Total_energy[beam_proton]/D");
    fEventTree->Branch("EndE",EndE,"EndE[beam_proton]/D");
    fEventTree->Branch("momentum",momentum,"momentum[beam_proton]/D");
    fEventTree->Branch("peakT_max",peakT_max,"peakT_max[beam_proton]/F");
    fEventTree->Branch("peakT_min",peakT_min,"peakT_min[beam_proton]/F");
    fEventTree->Branch("pdg",pdg,"pdg[beam_proton]/I");
    fEventTree->Branch("pdg_all_beam",pdg_all_beam,"pdg_all_beam[tot_trks]/I");
    fEventTree->Branch("origin",origin,"origin[beam_proton]/I");
    fEventTree->Branch("TrkID",TrkID,"TrkID[beam_proton]/I");
    fEventTree->Branch("TrackId_geant",TrackId_geant,"TrackId_geant[beam_proton]/I");
    fEventTree->Branch("EndPointx",EndPointx,"EndPointx[beam_proton]/F");
    fEventTree->Branch("EndPointy",EndPointy,"EndPointy[beam_proton]/F");
    fEventTree->Branch("EndPointz",EndPointz,"EndPointz[beam_proton]/F");
    fEventTree->Branch("EndPointx_corrected",EndPointx_corrected,"EndPointx_corrected[beam_proton]/F");
    fEventTree->Branch("EndPointy_corrected",EndPointy_corrected,"EndPointy_corrected[beam_proton]/F");
    fEventTree->Branch("EndPointz_corrected",EndPointz_corrected,"EndPointz_corrected[beam_proton]/F");
    fEventTree->Branch("StartPointx",StartPointx,"StartPointx[beam_proton]/F");
    fEventTree->Branch("StartPointy",StartPointy,"StartPointy[beam_proton]/F");
    fEventTree->Branch("StartPointz",StartPointz,"StartPointz[beam_proton]/F");
    fEventTree->Branch("trkstartcosxyz",trkstartcosxyz,"trkstartcosxyz[beam_proton][3]/F");
    fEventTree->Branch("trkendcosxyz",trkendcosxyz,"trkendcosxyz[beam_proton][3]/F");
    fEventTree->Branch("ntrkhits",ntrkhits,"ntrkhits[beam_proton][3]/I");
    fEventTree->Branch("trkdqdx",trkdqdx,"trkdqdx[beam_proton][3][3000]/F");
    fEventTree->Branch("trkdedx",trkdedx,"trkdedx[beam_proton][3][3000]/F");
    fEventTree->Branch("trkresrange",trkresrange,"trkresrange[beam_proton][3][3000]/F");
    fEventTree->Branch("trkhitx",trkhitx,"trkhitx[beam_proton][3][3000]/F");
    fEventTree->Branch("trkhity",trkhity,"trkhity[beam_proton][3][3000]/F");
    fEventTree->Branch("trkhitz",trkhitz,"trkhitz[beam_proton][3][3000]/F");
    fEventTree->Branch("trkpitch",trkpitch,"trkpitch[beam_proton][3][3000]/F");
    fEventTree->Branch("Efield",Efield,"Efield[beam_proton][3][3000]/F");
  }

  //========================================================================
  void protonanalysis::endJob(){     

  }

  //========================================================================
  void protonanalysis::beginRun(const art::Run&){
    mf::LogInfo("protonanalysis")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void protonanalysis::analyze( const art::Event& evt){
    reset();  

  
    
    std::vector<art::Ptr<recob::Track> > tracklist;
    auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel); 
    if (trackListHandle) art::fill_ptr_vector(tracklist, trackListHandle);

    std::vector<art::Ptr<recob::Hit>> hitlist;
    auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel); // to get information about the hits
    if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);
    art::FindManyP<recob::Hit> fmtht(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
   
    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime=tts.AsDouble();
     
    UInt_t year=0;
    UInt_t month=0;
    UInt_t day=0;
     
    year_month_date=tts.GetDate(kTRUE,0,&year,&month,&day);
     
    UInt_t hour=0;
    UInt_t min=0;
    UInt_t sec=0;
     
    hour_min_sec=tts.GetTime(kTRUE,0,&hour,&min,&sec);
  
    beam_proton=0;
    total_trks=0;
    tot_trks=0;



    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    std::cout<<detProp.Temperature()<<endl;
    std::cout<<"drift velocity "<<detProp.DriftVelocity(0.50,87.0)<<std::endl;
    double efield=detProp.PerPlaneEfield();
    double ycord=300;
    double zcord=350;
    for(int i=0;i<120;i++){
      double xcord=-360+i*6;
      geo::Vector_t fEfieldOffsets = SCE->GetEfieldOffsets(geo::Point_t{xcord,ycord,zcord});
      double Ef=efield+efield*fEfieldOffsets.X();
      EFieldX->SetBinContent(i+1, Ef);
 }

    for(int i=0;i<144;i++){
      for(int j=0;j<120;j++){
	for(int k=0;k<139;k++){
	  double x1=i*5.0-357.5;
	  double y1=j*5.0+2.5;
	  double z1=k*5.0+2.5;
	  geo::Vector_t fEfieldOffsets=SCE->GetEfieldOffsets(geo::Point_t{x1,y1,z1});
	  double EfX=efield+efield*fEfieldOffsets.X();
	  EField3DX->SetBinContent(i+1,j+1,k+1,EfX);
	}
      }
    }





    double EField1=efield;
    Float_t x1=-300;
    for(int i=0;i<139;i++){
      Float_t z=5*i+2.5;
      for(int j=0;j<120;j++){
	Float_t y=5*j+2.5;
	geo::Vector_t fEfieldOffsets = SCE->GetEfieldOffsets(geo::Point_t{x1,y,z});
	EField1 = std::sqrt((efield + efield*fEfieldOffsets.X())*(efield + efield*fEfieldOffsets.X())+(efield*fEfieldOffsets.Y())*(efield*fEfieldOffsets.Y())+(efield*fEfieldOffsets.Z())*(efield*fEfieldOffsets.Z()));
	electric_field_YZ->SetBinContent(i+1,j+1,EField1);
      }
    }
    Float_t x4=-25;
    for(int i=0;i<139;i++){
      Float_t z=5*i+2.5;
      for(int j=0;j<120;j++){
	Float_t y=5*j+2.5;
	geo::Vector_t fEfieldOffsets = SCE->GetEfieldOffsets(geo::Point_t{x4,y,z});
	EField1 = std::sqrt((efield + efield*fEfieldOffsets.X())*(efield + efield*fEfieldOffsets.X())+(efield*fEfieldOffsets.Y()*efield*fEfieldOffsets.Y())+(efield*fEfieldOffsets.Z()*efield*fEfieldOffsets.Z()));
	electric_field_YZ_3->SetBinContent(i+1,j+1,EField1);
      }
    }


    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);
      const recob::Track& track = *ptrack;
      // TVector3 pos, dir_start, dir_end, end;
      auto pos = track.Vertex();
      auto dir_start = track.VertexDirection();
      auto dir_end   = track.EndDirection();
      auto end = track.End();
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());     
     

     
      /*********************storing hits for each tracks*******************************/

      std::vector<art::Ptr<recob::Hit>> allHits=fmtht.at(i); //storing hits for ith track
      art::ServiceHandle<cheat::BackTrackerService> bt_serv;
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;    
      int trackid=-1;
      std::map<int,double> trkide;
      std::map<int,double> trknumelec;
      std::vector<float> hitpeakT;
      //selection based on hitpeakT
      for(size_t h1=0;h1<allHits.size();h1++){
	hitpeakT.push_back(allHits[h1]->PeakTime());
      }
      float max=-999999;
      float min=-999999;
      min=*min_element(hitpeakT.begin(),hitpeakT.end());
      max=*max_element(hitpeakT.begin(),hitpeakT.end());
      hitpeakT.clear();

      for(size_t h=0; h<allHits.size();h++){
	art::Ptr<recob::Hit> hit=allHits[h];
        std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(clockData, hit);
	for(size_t e=0;e<eveIDs.size(); ++e){
	  trkide[eveIDs[e].trackID] += eveIDs[e].energy;
	}
      }
      double  maxe = -1;
      // double tote = 0; // unused
      for(std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	//	cout<<" trkid "<<ii->first<<"  energy deposited = "<<ii->second<<endl;
	// tote += ii->second; // unused
	if((ii->second)>maxe){
	  maxe = ii->second;
	  trackid = ii->first;
	 
	}
      }

      //*******************new section************************//
      float total_energy=0.0;
      for(size_t h=0; h<allHits.size();h++){
	art::Ptr<recob::Hit> hit=allHits[h];
	if (hit->WireID().Plane!=2) continue;
        std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(clockData, hit);
	for(size_t e=0;e<eveIDs.size(); ++e){
	  if (eveIDs[e].trackID == trackid) total_energy +=  eveIDs[e].energy;
	}
      }
    
      int pdg1=-99999;
     
      const art::Ptr<simb::MCTruth> mc=pi_serv->TrackIdToMCTruth_P(trackid);
      const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trackid);
      if(!particle) continue;
      pdg1=particle->PdgCode();
      int origin_type=mc->Origin();
      // int no_daughters=particle->NumberDaughters();
      if((particle->Process()=="primary") && origin_type==4) pdg_beamparticles->Fill(pdg1);
     
      if(particle->Process()=="primary" && origin_type==4){
	tot_trks++;
	pdg_all_beam[tot_trks-1]=pdg1;
      }

      if(!(pdg1==2212 && (particle->Process()=="primary") && origin_type==4/* && no_daughters==0*/)) continue;//non interacting beam protons
     
      for(size_t h=0; h<allHits.size();h++){
	art::Ptr<recob::Hit> hit=allHits[h];
        std::vector<const sim::IDE*> eventIDs = bt_serv->HitToSimIDEs_Ps(clockData, hit);
	float charge=0;
	int number_electrons=0;
	if(hit->WireID().Plane!=2) continue;
	charge = hit->Integral();
	for(size_t e=0;e<eventIDs.size(); ++e){
	  //  trknumelec[eventIDs[e].trackID]+=eventIDs[e].numElectrons;
	  //std::cout<<eventIDs[e]->numElectrons<<" "<<eventIDs[e]->energy<<std::endl;
	  number_electrons+=eventIDs[e]->numElectrons;
	}
	calibration_constants->Fill(charge/float(number_electrons));
      }
     

    
     

      total_trks++;

      beam_proton++;
    




      /*  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();


      ****************************************************************************************************************/

      EndPointx_corrected[beam_proton-1]=particle->EndPosition()[0]-SCE->GetPosOffsets(geo::Point_t(particle->EndPosition()[0],particle->EndPosition()[1],particle->EndPosition()[2])).X();
      EndPointy_corrected[beam_proton-1]=particle->EndPosition()[1]+SCE->GetPosOffsets(geo::Point_t(particle->EndPosition()[0],particle->EndPosition()[1],particle->EndPosition()[2])).Y();
      EndPointz_corrected[beam_proton-1]=particle->EndPosition()[2]+SCE->GetPosOffsets(geo::Point_t(particle->EndPosition()[0],particle->EndPosition()[1],particle->EndPosition()[2])).Z();
      origin[beam_proton-1]=mc->Origin();
      EndE[beam_proton-1]=particle->EndE();
      momentum[beam_proton-1]=particle->Momentum().Vect().Mag();
      pdg[beam_proton-1]=particle->PdgCode();
      TrackId_geant[beam_proton-1]=particle->TrackId();
      EndPointx[beam_proton-1]=particle->EndPosition()[0];
      EndPointy[beam_proton-1]=particle->EndPosition()[1];
      EndPointz[beam_proton-1]=particle->EndPosition()[2];
      StartPointx[beam_proton-1]=particle->Vx();
      StartPointy[beam_proton-1]=particle->Vy();
      StartPointz[beam_proton-1]=particle->Vz();
      Total_energy[beam_proton-1]=total_energy;
      trackthetaxz[beam_proton-1]=theta_xz;
      trackthetayz[beam_proton-1]=theta_yz;
      trkstartx[beam_proton-1]=pos.X();
      trkstarty[beam_proton-1]=pos.Y();
      trkstartz[beam_proton-1]=pos.Z();
      trkendx[beam_proton-1]=end.X();
      trkendy[beam_proton-1]=end.Y();
      trkendz[beam_proton-1]=end.Z();
      trklen[beam_proton-1]=track.Length();
      peakT_max[beam_proton-1]=max;
      peakT_min[beam_proton-1]=min;
      TrkID[beam_proton-1]=track.ID();
      trkstartcosxyz[beam_proton-1][0]=dir_start.X();
      trkstartcosxyz[beam_proton-1][1]=dir_start.Y();
      trkstartcosxyz[beam_proton-1][2]=dir_start.Z();
      trkendcosxyz[beam_proton-1][0]=dir_end.X();
      trkendcosxyz[beam_proton-1][1]=dir_end.Y();
      trkendcosxyz[beam_proton-1][2]=dir_end.Z();
      for(size_t ical = 0; ical<calos.size(); ++ical){
	if(!calos[ical]) continue;
	if(!calos[ical]->PlaneID().isValid) continue;
	int planenum = calos[ical]->PlaneID().Plane;
	if(planenum<0||planenum>2) continue;
	const size_t NHits = calos[ical] -> dEdx().size();
	ntrkhits[beam_proton-1][planenum]=int(NHits);
	for(size_t iHit = 0; iHit < NHits; ++iHit){
	  const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	  trkdqdx[beam_proton-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
	  trkdedx[beam_proton-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
	  trkresrange[beam_proton-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	  trkhitx[beam_proton-1][planenum][iHit]=TrkPos.X();
	  trkhity[beam_proton-1][planenum][iHit]=TrkPos.Y();
	  trkhitz[beam_proton-1][planenum][iHit]=TrkPos.Z();
	  trkpitch[beam_proton-1][planenum][iHit]=(calos[ical]->TrkPitchVec())[iHit];

	  geo::Vector_t fEfieldOffsets = SCE->GetEfieldOffsets(geo::Point_t{TrkPos.X(),TrkPos.Y(),TrkPos.Z()});
	  Efield[beam_proton-1][planenum][iHit]=sqrt((efield + efield*fEfieldOffsets.X())*(efield + efield*fEfieldOffsets.X())+(efield*fEfieldOffsets.Y())*(efield*fEfieldOffsets.Y())+(efield*fEfieldOffsets.Z())*(efield*fEfieldOffsets.Z()));
	} // loop over iHit..
      } // loop over ical 2nd time...
    } // loop over trks...
    // stopping_muons->Fill(stopping_trks);
    fEventTree->Fill();
  } // end of analyze function
	   
  /////////////////// Defintion of reset function ///////////
  void protonanalysis::reset(){
    run = -99999;
    subrun = -99999;
    event = -99999;
    evttime = -99999;
    beam_proton = -99999;
    total_trks=-99999;
    year_month_date=-99999;
    hour_min_sec=-99999;
    for(int i=0; i<kMaxTracks; i++){
      trackthetaxz[i]=-99999;
      trackthetayz[i]=-99999;
      trkstartx[i]=-99999;
      trkstarty[i]=-99999;
      trkstartz[i]=-99999;
      trkendx[i]=-99999;
      trkendy[i]=-99999;
      trkendz[i]=-99999;
      trklen[i]=-99999;
      TrkID[i]=-99999;
      pdg_all_beam[i]=-99999;
      TrackId_geant[i]=-99999;
      origin[i]=-99999;
      EndE[i]=-99999;
      pdg[i]=-99999;
      Total_energy[i]=-999999;
      EndPointx[i]=-99999;
      EndPointy[i]=-99999;
      EndPointz[i]=-99999;
      EndPointx_corrected[i]=-99999;
      EndPointy_corrected[i]=-99999;
      EndPointz_corrected[i]=-99999;
      StartPointx[i]=-99999;
      StartPointy[i]=-99999;
      StartPointz[i]=-99999;
      pdg[i]=-99999;
      EndE[i]=-99999;
      momentum[i]=-99999;
      peakT_max[i]=-99999;
      peakT_min[i]=-99999;
      xprojectedlen[i]=-99999;
      trkstartcosxyz[i][0]=-99999;
      trkstartcosxyz[i][1]=-99999;
      trkstartcosxyz[i][2]=-99999; 
      trkendcosxyz[i][0]=-99999;
      trkendcosxyz[i][1]=-99999;
      trkendcosxyz[i][2]=-99999;
      ntrkhits[i][0] = -99999;
      ntrkhits[i][1] = -99999;
      ntrkhits[i][2] = -99999;
      for(int j=0; j<3; j++){
	for(int k=0; k<3000; k++){
	  trkdqdx[i][j][k]=-99999;
	  trkdedx[i][j][k]=-99999;
	  trkresrange[i][j][k]=-99999;
	  trkhitx[i][j][k]=-99999;
	  trkhity[i][j][k]=-99999;
	  trkhitz[i][j][k]=-99999;
	  trkpitch[i][j][k]=-99999;
	  Efield[i][j][k]=-99999;
	}
      }
    }
  }
  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(protonanalysis)
}
