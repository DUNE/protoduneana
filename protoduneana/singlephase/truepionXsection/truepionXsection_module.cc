#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/WireReadout.h"
//#include "larcore/Geometry/CryostatGeo.h"
//#include "larcore/Geometry/TPCGeo.h"
//#include "larcore/Geometry/PlaneGeo.h"
//#include "larcore/Geometry/WireGeo.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
//#include "larcoreobj/SimpleTypeAndConstants/RawTypes.h"
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
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
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
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

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

#include "TComplex.h"
#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include <map>

namespace protoana{

  class truepionXsection : public art::EDAnalyzer {
  public:

    explicit truepionXsection(fhicl::ParameterSet const& pset);
    virtual ~truepionXsection();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    //  double  distance(double, double, double, double, double, double );
  private:
    ///below are some addition to be able to use protodune utitilites
    // Helper utility functions
    // ProtoDUNEDataUtils dataUtil;
    // ProtoDUNEPFParticleUtils pfpUtil;
    // ProtoDUNETrackUtils trackUtil;
   

    ProtoDUNETruthUtils truthUtil;

    //fcl parameters
    const art::InputTag fBeamModuleLabel;
    std::string fCalorimetryTag;
    std::string fTrackerTag;
    std::string fShowerTag;
    std::string fPFParticleTag;
    std::string fGeneratorTag;
    ///////////////////////helper utility ends/////////


    TTree* fEventTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    // TH1D*   h_DE   ;
    // TH1D*   h_DX   ;
    //  TH1D*   h_DEDX ;
    TH1D*   h_DEUniform   ;
    TH1D*   h_DXUniform   ;
    TH1D*   h_DEDXUniform ;
    TH1D*   h_DeltaE ;
    // TH1D*   h_SimIDEDist ;
    TH1D*   h_UniformDistances ;
    TH1D *hInteractingKE; 
    TH1D *hInteractingKEEl; 
    TH1D *hInteractingKEElDep; 
    TH1D *hInteractingKEInel; 
    TH1D *hIncidentKE; 
    TH1D *hCrossSection;
    TH1D *hCrossSectionEl;
    TH1D *hCrossSectionInel;
    TH1D *hKEAtTPCFF; 
    TH1D *hInitialKE; 
    TH1D *hInitialPz; 
    TH2D *hXZ; 
    TH2D *hYZ; 
    TH2D *hXZPre; 
    TH2D *hYZPre; 
    TH2D *hdEVsdX; 
    TH2D *hdEVsKE; 
    TH1D *processtot;//number of processes for each mc particle
    TH1D *processEl;
    TH1D *processInel;
    //bool    debug = false;
    double trueVtxX ;
    double trueVtxY ;
    double trueVtxZ ;
    double trueEndX ;
    double trueEndY ;
    double trueEndZ ;
    double finalKE ;
    std::vector<std::string> G4Process; 
 
    double minX =  -360.0;
    double maxX = 360.0;
    double minY =0.0;
    double maxY = 600.0;
    double minZ =  0.0; // G10 at z=1.2
    double maxZ = 695.0;
    int evtsPreTPC  = 0;
    int evtsInTheMiddle = 0;
    int evtsPostTPC = 0;
    int throughgoing = 0;
    int interactingInTPC = 0;
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
  }; 

  //========================================================================
  truepionXsection::truepionXsection(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    ///new lines 
    //dataUtil(pset.get<fhicl::ParameterSet>("DataUtils")),
    fBeamModuleLabel(pset.get< art::InputTag >("BeamModuleLabel")),
    //fTrackerTag(pset.get<std::string>("TrackerTag")),
    //fShowerTag(pset.get<std::string>("ShowerTag")),
    //fPFParticleTag(pset.get<std::string>("PFParticleTag")),
    fGeneratorTag(pset.get<std::string>("GeneratorTag")),
    ///new lines end here

    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false))
  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
 
  //========================================================================
  truepionXsection::~truepionXsection(){
  }
  //========================================================================

  //========================================================================

  double distance(double x1, double y1, double z1, double x2, double y2, double z2)
  {
    double d = TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    return d;
  }




  void truepionXsection::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");

    //  h_DE   = tfs->make<TH1D>("h_DE","h_DE; Energy Deposited [MeV]",200, 0,100);   
    // h_DX   = tfs->make<TH1D>("h_DX","h_DX; Distance between points  [cm]",400, 0,20);   
    // h_DEDX = tfs->make<TH1D>("h_DEDEX","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   
    h_DEUniform   = tfs->make<TH1D>("h_DEUniform","h_DE; Energy Deposited [MeV]",200, 0,100);   
    h_DXUniform   = tfs->make<TH1D>("h_DXUniform","h_DX; Distance between points  [cm]",400, 0,20);   
    h_DEDXUniform = tfs->make<TH1D>("h_DEDEXUniform","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   
    h_DeltaE = tfs->make<TH1D>("h_DeltaE","h_DeltaE; dEDep - TrjDE [MeV/cm]",500, -1000,1000);   
    // h_SimIDEDist= tfs->make<TH1D>("h_SimIDEDist","h_SimIDEDist; h_SimIDEDist [cm]",1000, 0,10);   
  
    h_UniformDistances = tfs->make<TH1D>("h_UniformDistances","h_UniformDistances; Distance between uniform points  [cm]",500, 0,5);   
    hInitialPz     = tfs->make<TH1D>("hInitialPz"    , "Initial Pz [MeV/c]"    , 42, -100, 2000);
    hInitialKE     = tfs->make<TH1D>("hInitialKE"    , "Initial Kinetic Energy [MeV]"    , 42, -100, 2000);
    hKEAtTPCFF     = tfs->make<TH1D>("hKEAtTPCFF"    , "Kinetic Energy @ TPC FF [MeV]"   , 42, -100, 2000);
    hIncidentKE        = tfs->make<TH1D>("hIncidentKE"   , "Incident Kinetic Energy [MeV]"   , 42, -100, 2000); 

    processtot=tfs->make<TH1D>("processtot","Total number of interaction for each particle;number of interaction;number of entries",50,0,50);
    processEl=tfs->make<TH1D>("processEl","Total number of Elastic interaction for each particle;number of interaction;number of entries",50,0,50);
    processInel=tfs->make<TH1D>("processInel","Total number of InElastic interaction for each particle;number of interaction;number of entries",50,0,50);
    hInteractingKE     = tfs->make<TH1D>("hInteractingKE", "Interacting Kinetic Energy [MeV]", 42, -100, 2000);
    hInteractingKEEl   = tfs->make<TH1D>("hInteractingKEEl", "Elastic Interacting Kinetic Energy [MeV]", 42, -100, 2000); 
    hInteractingKEElDep   = tfs->make<TH1D>("hInteractingKEElDep", "Dep Elastic Interacting Kinetic Energy [MeV]", 42, -100, 2000); 
    hInteractingKEInel = tfs->make<TH1D>("hInteractingKEInel", "Inelastic Interacting Kinetic Energy [MeV]", 42, -100, 2000);
    hCrossSection     = tfs->make<TH1D>("hCrossSection"     , "Cross-Section [barn]"             , 42, -100, 2000);
    hCrossSectionEl   = tfs->make<TH1D>("hCrossSectionEl"   , "Elastic Cross-Section [barn]"     , 42, -100, 2000);
    hCrossSectionInel = tfs->make<TH1D>("hCrossSectionInel" , "Inelastic Cross-Section [barn]"   , 42, -100, 2000);
    hXZ     = tfs->make<TH2D>("hXZ"     , "hXZ"  , 895, -200, 695  , 720, -360, 360);  
    hYZ     = tfs->make<TH2D>("hYZ"     , "hYZ"  , 895, -200, 695  , 600, 0, 600); 
    hXZPre  = tfs->make<TH2D>("hXZPre"  , "hXZPre", 200, -100, 100 , 200, -100, 100); 
    hYZPre  = tfs->make<TH2D>("hYZPre"  , "hYZPre", 200, -100, 100 , 200, 300, 500); 
    hdEVsdX  = tfs->make<TH2D>("hdEVsdX"  , "hdEVsdX" , 504, -1, 50, 1100, -10, 100); 
    hdEVsKE  = tfs->make<TH2D>("hdEVsKE"  , "hdEVsKE" , 504, -1, 50, 220,  -10, 1000); 
  
    fEventTree->Branch("trueVtxX" ,&trueVtxX ,"trueVtxX/D");
    fEventTree->Branch("trueVtxY" ,&trueVtxY ,"trueVtxY/D");
    fEventTree->Branch("trueVtxZ" ,&trueVtxZ ,"trueVtxZ/D");
    fEventTree->Branch("trueEndX" ,&trueEndX ,"trueEndX/D");
    fEventTree->Branch("trueEndY" ,&trueEndY ,"trueEndY/D");
    fEventTree->Branch("trueEndZ" ,&trueEndZ ,"trueEndZ/D");
    fEventTree->Branch("finalKE"  ,&finalKE  ,"finalKE/D" );
    fEventTree->Branch("G4Process",&G4Process);


  }

 
  //========================================================================
  void truepionXsection::beginRun(const art::Run&){
    mf::LogInfo("truepionXsection")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void truepionXsection::analyze( const art::Event& evt){
    reset();  
    bool verbose=false;
    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist=pi_serv->ParticleList();
    // sim::ParticleList::const_iterator itPartp=plist.begin();
    bool keepInteraction=false;

    //for(size_t iPart=0;(iPart<plist.size()) && (plist.begin()!=plist.end());++iPart){  //particle loop begins here
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    const simb::MCParticle* particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);

    //const simb::MCParticle* particle=(itPartp++)->second;
    if(particle!=0x0){
      std::cout<<"particle pdg code is "<<particle->PdgCode()<<std::endl;
      if(particle->PdgCode()==211 && particle->Process()=="primary"){
	std::cout<<"found a pi+"<<std::endl;
	simb::MCTrajectory truetraj=particle->Trajectory();
	auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
        geo::View_t view = wireReadout.Plane({0, 0, 2}).View();
	auto simIDE_prim=bt_serv->TrackIdToSimIDEs_Ps(particle->TrackId(),view);
	std::map<double, sim::IDE> orderedSimIDE;
	for (auto& ide : simIDE_prim) orderedSimIDE[ide->z]= *ide;
	std::string interactionLabel="";
	double mass=particle->Mass();
	if(verbose) std::cout<<"mass "<<mass<<"\n";
	//Store the kinetic energy and momentum on z at WC4. Just for cross check 
	auto inTPCPoint  = truetraj.begin(); 
	auto Momentum0   = inTPCPoint->second;
	double KE0 = 1000*(TMath::Sqrt(Momentum0.X()*Momentum0.X() + Momentum0.Y()*Momentum0.Y() + Momentum0.Z()*Momentum0.Z() + mass*mass ) - mass); //I want this in MeV
	hInitialKE->Fill(KE0);
	hInitialPz->Fill(1000*Momentum0.Z());
	if(verbose) std::cout<<KE0;

	//--------------------------------------------------------
	// Identify the first trajectory point inside the TPC
	// Loop From First TrajPoint --> First Point in TPC 
	// Stop when you get into the TPC
	for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++)
	  {
	    auto pos = t->first;
	    if (pos.Z() < minZ) continue;
	    else if (pos.X() < minX || pos.X() > maxX ) continue;
	    else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	    else {
	      inTPCPoint = t;
	      break;
	    }
	  }// End search for first point in TPC
	evtsPreTPC++;
	hXZPre->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).X());
	hYZPre->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).Y());
	if (inTPCPoint !=truetraj.begin()){
	evtsPostTPC++;     
	hXZ   ->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).X());
	hYZ   ->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).Y());

	// Identify the last interesting trajectory point in TPC
	auto finTPCPoint = std::prev(truetraj.end()); 
	// The last point is a bit more complicated:
	// if there's no interaction, then it is simply the last point in the TPC
	// if there's one or more interaction points, it's the first interaction point deemed interesting (no coulomb)
	// Take the interaction Map... check if there's something there
	auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
    
	int el=0;
	int inel=0;
	int totint=0;
	if (thisTracjectoryProcessMap.size())
	  {
	    for(auto const& couple: thisTracjectoryProcessMap) 
	      { 
		totint++;
		// I'm not interested in the CoulombScattering, discard this case
		if(!verbose) std::cout<<(truetraj.KeyToProcess(couple.second))<<" Position "<<((truetraj.at(couple.first)).first).Z()<<"  Momentum "<<((truetraj.at(couple.first)).second).Z()<<std::endl;
		if ((truetraj.KeyToProcess(couple.second)).find("CoulombScat")!= std::string::npos) continue;
              
		// Let's check if the interaction is in the the TPC
		auto     interactionPos4D =  (truetraj.at(couple.first)).first ;        
		if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
		else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
		else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;
		if ((truetraj.KeyToProcess(couple.second)).find("hadElastic")!= std::string::npos) el++;
		if ((truetraj.KeyToProcess(couple.second)).find("pi+InElastic")!= std::string::npos) inel++;

		// If we made it here, then this is the first interesting interaction in the TPC
		// Our job is done!!! Great! Store the interaction label and the iterator for the final point
		interactionLabel = truetraj.KeyToProcess(couple.second);
		std::cout<<"interaction Label "<<interactionLabel<<" EndProcess "<<particle->EndProcess()<<std::endl;
		finTPCPoint = truetraj.begin() + couple.first; 
		keepInteraction=true;
		interactingInTPC++;
		break; //commented for now
	      }// Loop on interaction points
	  } // If there are G4 interactions 

	std::cout<<"this trjectory size "<<thisTracjectoryProcessMap.size()<<std::endl;
	processtot->Fill(totint);
	processEl->Fill(el);
	processInel->Fill(inel);
  

	sim::ParticleList::const_iterator itPart=plist.begin();
  
	// If I didn't find anything interesting in the intereaction map, let's loop back!
	if (!keepInteraction)
	  {
	    //Loop on the daughters 
	    for(size_t iPart1=0;(iPart1<plist.size()) && (plist.begin()!=plist.end());++iPart1)
	      {
		const simb::MCParticle* pPart=(itPart++)->second;
		//   art::Ptr<simb::MCTruth> const& mcDaught=pi_serv->ParticleToMCTruth_P(pPart)
		//We keep only the dauthers of the primary not coming from elastic or inelastic scattering
		if (pPart->Mother()  != 1 ) continue;
		//  std::cout<<"pPart mother "<<pPart->Mother()<<std::endl;
		if ((pPart->Process()).find("astic")!= std::string::npos) continue;
		if ((pPart->Process()).find("CoulombScat")!= std::string::npos) continue;
		//Is the daughter born inside the TPC? If yes, store the process which created it 
		simb::MCTrajectory trueDaugthTraj = pPart->Trajectory();              
		if (trueDaugthTraj.begin()->first.Z() < minZ || trueDaugthTraj.begin()->first.Z() > maxZ) continue;
		else if (trueDaugthTraj.begin()->first.X() <   minX || trueDaugthTraj.begin()->first.X() > maxX ) continue;
		else if (trueDaugthTraj.begin()->first.Y() <   minY || trueDaugthTraj.begin()->first.Y() > maxY ) continue;
		else {
		  interactionLabel = pPart->Process();
	    
		  break;
		}
	      }        
	    for ( auto t = std::prev(truetraj.end()); t!= truetraj.begin(); t--)
	      {
		auto pos = t->first;
	    
		if (pos.Z() > maxZ) continue;
		else if (pos.X() <   minX || pos.X() > maxX ) continue;
		else if (pos.Y() <   minY || pos.Y() > maxY ) continue;
		else {
		  finTPCPoint = t;
		  break;
		}
	      }
	  }
	//}//new bracket added here

	if (finTPCPoint != inTPCPoint){
	  auto posFin = finTPCPoint->first;
	  auto posIni = inTPCPoint->first;
	  //Let's record what the initial and final points are.
	  trueVtxX = posIni.X();
	  trueVtxY = posIni.Y();
	  trueVtxZ = posIni.Z();
	  trueEndX = posFin.X();
	  trueEndY = posFin.Y();
	  trueEndZ = posFin.Z();
	  // std::cout<<"initial x, y and z && final x, y and z "<<trueVtxX<<" "<<trueVtxY<<" "<<trueVtxZ<<" "<<trueEndX<<" "<<trueEndY<<" "<<trueEndZ<<std::endl;
	  auto totLength = distance(posFin.X(), posFin.Y(), posFin.Z(),posIni.X(), posIni.Y(), posIni.Z() );
	  // std::cout<<"totalLength "<<totLength<<std::endl;
      
	  // We want to chop up the points between the first and last uniformly
	  // and ordered by Z
	  // Order them in order of increasing Z
	  std::map<double, TVector3> orderedUniformTrjPts;
	  // We want the first and uniform point to coincide with the 
	  // the first and last points we just found 
	  auto positionVector0 = (inTPCPoint ->first).Vect(); 
	  auto positionVector1 = (finTPCPoint->first).Vect(); 
	  orderedUniformTrjPts[positionVector0.Z()] = positionVector0;
	  orderedUniformTrjPts[positionVector1.Z()] = positionVector1;
	  // const double trackPitch = 0.47;
	  const double trackPitch = 0.4792;
	  // I do have space for at least one extra point, so let's put it there!
	  // Calculate how many extra points I need to put between the new first point and the second TrajPoint
	  int    nPts            = (int) (totLength/trackPitch);
	  for (int iPt = 1; iPt <= nPts; iPt++ )
	    {
	      auto newPoint = positionVector0 + iPt*(trackPitch/totLength) * (positionVector1 - positionVector0);
	      orderedUniformTrjPts[newPoint.Z()] = newPoint;
	    }
	  // If the distance between the last point and the second to last is less then 0.235
	  // eliminate the second to last point
	  auto lastPt         = (orderedUniformTrjPts.rbegin())->second;
	  auto secondtoLastPt = (std::next(orderedUniformTrjPts.rbegin()))->second;
	  double lastDist = distance(lastPt.X(),lastPt.Y(),lastPt.Z(),secondtoLastPt.X(),secondtoLastPt.Y(),secondtoLastPt.Z());
	  if (lastDist < 0.240)
	    {
	      orderedUniformTrjPts.erase((std::next(orderedUniformTrjPts.rbegin()))->first );
	    } 
    
	  // Calculate the initial kinetic energy
	  auto initialMom =     inTPCPoint->second;
	  double initialKE = 1000*(TMath::Sqrt(initialMom.X()*initialMom.X() + initialMom.Y()*initialMom.Y() + initialMom.Z()*initialMom.Z() + mass*mass ) - mass); 
	  hKEAtTPCFF->Fill(initialKE);
	  double kineticEnergy = initialKE;
	  auto old_it = orderedUniformTrjPts.begin();
	  for (auto it = std::next(orderedUniformTrjPts.begin()); it != orderedUniformTrjPts.end(); it++, old_it++ )
	    {
          
	      if (verbose)  std::cout << it->first<<" : " << (it->second).Z() << std::endl ;
	      auto oldPos        = old_it->second;
	      auto currentPos    =     it->second;
        
	      double uniformDist =  (currentPos - oldPos).Mag();
	      h_UniformDistances->Fill(uniformDist);
          
	      //Calculate the energy deposited in this slice          
	      auto old_iter = orderedSimIDE.begin();
	      double currentDepEnergy = 0.;
	      for ( auto iter= orderedSimIDE.begin(); iter!= orderedSimIDE.end(); iter++,old_iter++)
		{
		  auto currentIde = iter->second;
		  // std::cout<<"Z position of the trajectory hits "<<currentIde.z<<std::endl;
		  if ( currentIde.z < oldPos.Z()) continue;
		  if ( currentIde.z > currentPos.Z()) continue;
		  currentDepEnergy += currentIde.energy;
		}// Determing which simIDE is within the current slice
	      // avoid overfilling super tiny energy depositions
	      if (currentDepEnergy/uniformDist < 0.1 ) continue;
	      //Calculate the current kinetic energy
	      kineticEnergy -= currentDepEnergy;
	      hdEVsdX->Fill(currentDepEnergy,(currentPos.Z()-oldPos.Z()) );
	      hdEVsKE->Fill(currentDepEnergy,kineticEnergy);
	      hIncidentKE->Fill(kineticEnergy);
	      h_DEUniform->Fill(currentDepEnergy);
	      h_DXUniform->Fill(uniformDist);
	      h_DEDXUniform->Fill(currentDepEnergy/uniformDist);
          
	    }// Loop on OrderedPoints

	  if(interactionLabel.find("Inelastic")!= std::string::npos )//added pi+Inelastic instead of Inelastic
	    //if(interactionLabel=="pi+Inelastic")
	    // if(particle->EndProcess()=="pi+Inelastic")
	    {
	      // std::cout<<"Interaction Label: "<<interactionLabel<<"  EndProcess "<<particle->EndProcess(); 
	      hInteractingKE->Fill(kineticEnergy);
	      hInteractingKEInel->Fill(kineticEnergy);
	      //std::cout<<"inside pi+Inelastic loop and ke "<<kineticEnergy<<std::endl;
	    }
      
	  //Fill the Elastic and Total Interacting with the last point
	  if ( interactionLabel.find("Elastic")!= std::string::npos ) //added hadElastic instead of Elastic
	    //   if(interactionLabel=="hadElastic")
	    // if(particle->EndProcess()=="hadElastic")
	    {
	      //  std::cout<<"Interaction Label: "<<interactionLabel<<"  EndProcess "<<particle->EndProcess(); 
	      h_DeltaE ->Fill(kineticEnergy -  1000*((finTPCPoint->second).E() - mass) );
	      hInteractingKEElDep->Fill(kineticEnergy);
	      hInteractingKE->Fill(kineticEnergy);
	      //	auto MomentumF = finTPCPoint->second;
	      //	double KEF = 1000*(TMath::Sqrt(MomentumF.X()*MomentumF.X() + MomentumF.Y()*MomentumF.Y() + MomentumF.Z()*MomentumF.Z() + mass*mass ) - mass); //I want this in MeV
	      //	hInteractingKEEl->Fill(KEF);
	      hInteractingKEEl->Fill(kineticEnergy);
	      //	std::cout<<"inside hadElastic Loop  and KEF "<<KEF<<std::endl;
	    }
	  finalKE = kineticEnergy;
	  keepInteraction = false;
	  if (!interactionLabel.size())
	    {
	      throughgoing++;
	      G4Process.push_back("throughgoing");
	    }else
	    {
	      G4Process.push_back(interactionLabel);
	    }



	  ////***********************new segment of the code ends here******************************************************************************//
	}
	}//if finTPC!=inTPC
      }//if primary
    } //if particle 
    fEventTree->Fill();
  } // end of analyze function

    //========================================================================
  void truepionXsection::endJob(){     

    std::cout<<"-------------------------------------------"<<std::endl;
    std::cout<<"True Events pre-TPC .............. "<<evtsPreTPC<<std::endl;
    std::cout<<"True Events pre-TPC .............. "<<evtsInTheMiddle<<std::endl;
    std::cout<<"True Events post-TPC ............. "<<evtsPostTPC<<std::endl;
    std::cout<<"True Throughgoing    ............. "<<throughgoing<<std::endl;
    std::cout<<"True interactingInTPC ............ "<<interactingInTPC<<std::endl;
    std::cout<<"-------------------------------------------"<<std::endl;
    float rho = 1396; //kg/m^3
    float molar_mass = 39.95; //g/mol
    float g_per_kg = 1000; 
    float avogadro = 6.022e+23; //number/mol
    float number_density = rho*g_per_kg/molar_mass*avogadro;
    //float slab_width = 0.0047;//in m
    float slab_width = 0.004792;//in m
    // Calculate the Cross Section
    // ###################################################################
    // #### Looping over the exiting bins to extract the cross-section ###
    // ###################################################################
    for( int iBin = 1; iBin <= hInteractingKE->GetNbinsX(); ++iBin )
      {
	// ### If an incident bin is equal to zero then skip that bin ###
	if( hIncidentKE->GetBinContent(iBin) == 0 )continue; //Temporary fix to ensure that no Infinities are propagated to pad
   
	// ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
	float TempCrossSection = (hInteractingKE->GetBinContent(iBin)/hIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width);
   
	float elCrossSection   = ((hInteractingKEEl->GetBinContent(iBin)/hIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width) ) * (1/1e-28);
	float inelCrossSection = ((hInteractingKEInel->GetBinContent(iBin)/hIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width) ) * (1/1e-28);
	// ### Convert this into Barns ###
	float crossSection = TempCrossSection * (1/1e-28); 
        
	// ### Putting the value on the plot
	hCrossSection    ->SetBinContent(iBin,crossSection);
	hCrossSectionEl  ->SetBinContent(iBin,elCrossSection);
	hCrossSectionInel->SetBinContent(iBin,inelCrossSection);
	// ###########################################################
	// ### Calculating the error on the numerator of the ratio ###
	// ###########################################################
	float denomError = pow(hIncidentKE->GetBinContent(iBin),0.5);
	float denom = hIncidentKE->GetBinContent(iBin);
	if(denom == 0) continue; 
	float term2 = denomError/denom;
      
	float numError = pow(hInteractingKE->GetBinContent(iBin),0.5);
	float num = hInteractingKE->GetBinContent(iBin);
	float numErrorEl = pow(hInteractingKEEl->GetBinContent(iBin),0.5);
	float numEl = hInteractingKEEl->GetBinContent(iBin);
	float numErrorInel = pow(hInteractingKEInel->GetBinContent(iBin),0.5);
	float numInel = hInteractingKEInel->GetBinContent(iBin);
	// ### Putting in a protection against dividing by zero ###   
	if(num != 0){
	  float term1 = numError/num;
	  float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);
	  hCrossSection->SetBinError(iBin,totalError);
	}
      
	if(numEl != 0){
	  float term1 = numErrorEl/numEl;
	  float totalError = (elCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width)  *(1e26);
	  hCrossSectionEl->SetBinError(iBin,totalError);
	}
      
	if(numInel != 0){
	  float term1 = numErrorInel/numInel;
	  float totalError = (inelCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width)  *(1e26);
	  hCrossSectionInel->SetBinError(iBin,totalError);
	}
      }//<---End iBin Loop


  }



	   
  /////////////////// Defintion of reset function ///////////
  void truepionXsection::reset(){
    run = -99999;
    subrun = -99999;
    event = -99999;
 
    trueVtxX = -999.;
    trueVtxY = -999.;
    trueVtxZ = -999.;
    trueEndX = -999.;
    trueEndY = -999.;
    trueEndZ = -999.;
    trueEndZ = -999.;
    finalKE  = -999.;
    G4Process.clear(); 


  }








  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(truepionXsection)
}
