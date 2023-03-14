////////////////////////////////////////////////////////////////////////
// Class:       protonmc
// File:        protonmc_module.cc
//
// Extract protoDUNE useful information, do a firs tpre-selection and save output to a flat tree
// 
// Some parts are copied from the beam module example
//
// Georgios Christodoulou - georgios.christodoulou at cern.ch
// Heng-Ye Liao modified it for his protonanalysis - liao@phys.ksu.edu
//
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/Exceptions.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
//#include "art/Framework/Services/Optional/TFileService.h"
#include "art_root_io/TFileService.h"                      
//#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcore/Geometry/WireReadout.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
//#include "dunesim/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
//#include "dunesim/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"

#include "lardata/ArtDataHelper/MVAReader.h"

//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNETrackUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEShowerUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNETruthUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEPFParticleUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

#include "larsim/MCCheater/ParticleInventoryService.h"

//#g4reweight
#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"
#include "geant4reweight/src/ReweightBase/G4ReweighterFactory.hh"
#include "geant4reweight/src/ReweightBase/G4Reweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/PropBase/G4ReweightParameterMaker.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "protoduneana/Utilities/G4ReweightUtils.h"



// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TTimeStamp.h"

// C++ Includes
#include <stdio.h>
#include <stdlib.h> 
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <tuple>

// Maximum number of beam particles to save
const int NMAXDAUGTHERS = 15;

//const int kMaxTracks  = 1000;
//const int kMaxHits = 10000;
//double prim_energy=0.0;

namespace protoana {
	class protonmc;
}

using protoana::G4ReweightUtils::CreateRWTraj;



class protoana::protonmc : public art::EDAnalyzer {
	public:

		explicit protonmc(fhicl::ParameterSet const & p);

		protonmc(protonmc const &) = delete;
		protonmc(protonmc &&) = delete;
		protonmc & operator = (protonmc const &) = delete;
		protonmc & operator = (protonmc &&) = delete;

		virtual void beginJob() override;
		virtual void endJob() override;

		// Required functions.
		void analyze(art::Event const & evt) override;

	private:

		// Helper utility functions
		protoana::ProtoDUNEDataUtils dataUtil;
		protoana::ProtoDUNEPFParticleUtils pfpUtil;
		protoana::ProtoDUNETrackUtils trackUtil;
		protoana::ProtoDUNETruthUtils truthUtil;

		// Initialise tree variables
		void Initialise();

		// Fill cosmics tree
		void FillCosmicsTree(art::Event const & evt, std::string pfParticleTag);

		// fcl parameters
		//const art::InputTag fNNetModuleLabel; //label of the module used for CNN tagging  
		const art::InputTag fBeamModuleLabel;
		std::string fCalorimetryTag;
		std::string fTrackerTag;
		std::string fShowerTag;
		std::string fPFParticleTag;
		std::string fGeneratorTag;
		//std::string fMuMCSInputTag;
		//trkf::TrajectoryMCSFitter fMCSFitter; /// objet holding the configuration of the ProtoDUNE-SP MCS fitting alg

		//geometry
		geo::GeometryCore const * fGeometry;
                geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

		double MCTruthT0, TickT0;
		int nT0s;   

		TTree *fPandoraBeam;
		//TTree *fPandoraCosmics;

		// Tree variables
		int fRun;
		int fSubRun;
		int fevent;
		double fTimeStamp;
		int fNactivefembs[6];

		// Beam track
		int fbeamtrigger;
		double ftof;
		int fcerenkov1;
		int fcerenkov2;
		double fbeamtrackMomentum;
		double fbeamtrackP[3]; //Px/Py/Pz
		double fbeamtrackEnergy;
		double fbeamtrackPos[3];
		double fbeamtrackDir[3];
		double fbeamtrackTime;
		int fbeamtrackPdg;
		int fbeamtrackID;
		//int NumberBeamTrajectoryPoints;
		std::vector<double> beamtrk_x;
		std::vector<double> beamtrk_y;
		std::vector<double> beamtrk_z;
		std::vector<double> beamtrk_Px;
		std::vector<double> beamtrk_Py;
		std::vector<double> beamtrk_Pz;
		std::vector<double> beamtrk_Eng;


		//beam info from spectrometer
		std::vector<double> beamMomentum_spec;
  		std::vector<double> beamPosx_spec;
  		std::vector<double> beamPosy_spec;
  		std::vector<double> beamPosz_spec;
  
  		std::vector<double> beamDirx_spec;
  		std::vector<double> beamDiry_spec;
  		std::vector<double> beamDirz_spec;


		//Front face KE & Positions
		double ke_ff;
		std::vector<double> pos_ff;

		//KE_true and KE_reco
		std::vector<double> primtrk_ke_true;
		std::vector<double> primtrk_ke_reco;
		std::string primtrk_end_g4process;
		std::vector<double> primtrk_range_true;

		//MCS
		//std::vector<float> primtrk_mcs_angles_reco;
		//std::vector<float> primtrk_mcs_radlengths_reco;
		//std::vector<float> primtrk_mcs_angles_true;
		//std::vector<float> primtrk_mcs_radlengths_true;

		//true info
		double ke_ff_true;
		//std::vector<unsigned char> primtrk_hit_processkey;
		std::vector< std::vector<double> > primtrk_hitx_true;
		std::vector< std::vector<double> > primtrk_hity_true;
		std::vector< std::vector<double> > primtrk_hitz_true;
		std::vector< std::vector<double> > primtrk_trkid_true;
		std::vector< std::vector<double> > primtrk_edept_true;

		//save all the true info  
		std::vector< std::vector<double> > primtrk_true_x;
		std::vector< std::vector<double> > primtrk_true_y;
		std::vector< std::vector<double> > primtrk_true_z;
		std::vector< std::vector<double> > primtrk_true_trkid;
		std::vector< std::vector<double> > primtrk_true_edept;

		std::vector< std::vector<int> > primtrk_true_wid;
		//std::vector< std::vector<int> > primtrk_true_tpc;
		//std::vector< std::vector<int> > primtrk_true_plane;

		//mctraj info
		std::vector< std::vector<double> > primtrj_true_x;
		std::vector< std::vector<double> > primtrj_true_y;
		std::vector< std::vector<double> > primtrj_true_z;
		std::vector< std::vector<double> > primtrj_true_edept;

		//hits from pfparticles
		std::vector<double> pfphit_peaktime_c;
		std::vector<double> pfphit_peaktime_u;
		std::vector<double> pfphit_peaktime_v;
		std::vector<double> pfphit_wireid_c;
		std::vector<double> pfphit_wireid_u;
		std::vector<double> pfphit_wireid_v;


		// Reconstructed tracks/showers
		double fvertex[3];
		double fsecvertex[3];

		int fisprimarytrack;
		int fisprimaryshower;
		double fprimaryBDTScore;
		int fprimaryNHits;
		double fprimaryTheta;
		double fprimaryPhi;
		double fprimaryLength;
		double fprimaryMomentum;
		double fprimaryEndMomentum;
		double fprimaryEndPosition[3];
		double fprimaryStartPosition[3];
		double fprimaryEndDirection[3];
		double fprimaryStartDirection[3];
		double fprimaryOpeningAngle;
		int fprimaryShowerBestPlane;
		double fprimaryShowerEnergy;
		double fprimaryShowerMIPEnergy;
		double fprimaryShowerdEdx;
		double fprimaryMomentumByRangeProton;
		double fprimaryMomentumByRangeMuon;
		double fprimaryKineticEnergy[3];
		double fprimaryRange[3];
		int fprimaryID;
		double fprimaryT0;

		//Truth info
		int ftruthpdg;  
		int fprimary_truth_TrackId;
		int fprimary_truth_Pdg;
		double fprimary_truth_StartPosition[4];
		double fprimary_truth_StartPosition_MC[4];
		double fprimary_truth_EndPosition[4];
		double fprimary_truth_EndPosition_MC[4];
		double fprimary_truth_P;
		double fprimary_truth_Momentum[4];
		double fprimary_truth_EndMomentum[4];
		//std::string G4Process_Primary_End;
		double fprimary_truth_Pt;
		double fprimary_truth_Mass;
		double fprimary_truth_Theta;
		double fprimary_truth_Phi;
		//int fprimary_truth_Process;
		//std::string fprimary_truth_Process;
		//std::string fprimary_truth_EndProcess;
		//std::string fprimary_backtrker_truth_Process;
		std::string fprimary_truth_EndProcess;
		int fprimary_truth_Isbeammatched;
		int fprimary_truth_NDaughters;

		//carlo info
		std::vector< std::vector<double> > primtrk_dqdx;
		std::vector< std::vector<double> > primtrk_resrange;
		std::vector< std::vector<double> > primtrk_dedx;
		std::vector<double> primtrk_range;
		std::vector< std::vector<double> > primtrk_hitx;
		std::vector< std::vector<double> > primtrk_hity;
		std::vector< std::vector<double> > primtrk_hitz;
		std::vector< std::vector<double> > primtrk_pitch;
		std::vector< std::vector<int> > primtrk_wid;

		//hit and CNN score

		//cnn score
		//std::vector<float> inelscore_c;
		//std::vector<float> elscore_c;
		//std::vector<float> nonescore_c;
		//associtated (x,y,z) 
		std::vector<float> x_c;
		std::vector<float> y_c;
		std::vector<float> z_c;


		//(x,y,z) of highest CNN score
		//double xyz_inelscore_c[3];
		//double xyz_elscore_c[3];


		//traj hit(x,y,z) to hit(wid,pt)
		std::vector<double> wid_c;
		std::vector<double> tt_c;
		std::vector<int> ch_c;
		std::vector<double> wirez_c;
		std::vector<double> pt_c;
		std::vector<double> q_c;
		std::vector<double> a_c;
		std::vector<int> wireno_c;

		std::vector<double> wid_v;
		std::vector<double> tt_v;
		std::vector<int> ch_v;
		std::vector<double> wirez_v;
		std::vector<double> pt_v;
		std::vector<double> q_v;
		std::vector<double> a_v;
		std::vector<int> wireno_v;

		std::vector<double> wid_u;
		std::vector<double> tt_u;
		std::vector<int> ch_u;
		std::vector<double> wirez_u;
		std::vector<double> pt_u;
		std::vector<double> q_u;
		std::vector<double> a_u;
		std::vector<int> wireno_u;


		//info of interactions  std::vector<double> interactionX;
		std::vector<double> interactionX;
		std::vector<double> interactionY;
		std::vector<double> interactionZ;
		std::vector<double> interactionE;
		std::vector<std::string> interactionProcesslist;
		std::vector<double> interactionAngles;
		std::vector<double> Zintersection;
		std::vector<double> Yintersection;
		std::vector<double> Xintersection;
		std::vector<double> timeintersection;

		std::vector<double> interaction_wid_c;
		std::vector<double> interaction_tt_c;
		std::vector<double> interaction_wid_v;
		std::vector<double> interaction_tt_v;
		std::vector<double> interaction_wid_u;
		std::vector<double> interaction_tt_u;

		// Daughters from primary
		int fNDAUGHTERS;
		int fisdaughtertrack[NMAXDAUGTHERS];
		int fisdaughtershower[NMAXDAUGTHERS];
		int fdaughterNHits[NMAXDAUGTHERS];
		double fdaughterTheta[NMAXDAUGTHERS];
		double fdaughterPhi[NMAXDAUGTHERS];
		double fdaughterLength[NMAXDAUGTHERS];
		double fdaughterMomentum[NMAXDAUGTHERS];
		double fdaughterEndMomentum[NMAXDAUGTHERS];
		double fdaughterEndPosition[NMAXDAUGTHERS][3];
		double fdaughterStartPosition[NMAXDAUGTHERS][3];
		double fdaughterEndDirection[NMAXDAUGTHERS][3];
		double fdaughterStartDirection[NMAXDAUGTHERS][3];
		double fdaughterOpeningAngle[NMAXDAUGTHERS];
		double fdaughterShowerEnergy[NMAXDAUGTHERS];
		double fdaughterShowerMIPEnergy[NMAXDAUGTHERS];
		double fdaughterShowerdEdx[NMAXDAUGTHERS];
		int fdaughterShowerBestPlane[NMAXDAUGTHERS];
		double fdaughterMomentumByRangeProton[NMAXDAUGTHERS];
		double fdaughterMomentumByRangeMuon[NMAXDAUGTHERS];
		double fdaughterKineticEnergy[NMAXDAUGTHERS][3];
		double fdaughterRange[NMAXDAUGTHERS][3];
		int fdaughterID[NMAXDAUGTHERS];
		//double fdaughterT0[NMAXDAUGTHERS];

		int fdaughter_truth_TrackId[NMAXDAUGTHERS];
		int fdaughter_truth_Pdg[NMAXDAUGTHERS];
		double fdaughter_truth_StartPosition[NMAXDAUGTHERS][4];
		double fdaughter_truth_EndPosition[NMAXDAUGTHERS][4];
		double fdaughter_truth_P[NMAXDAUGTHERS];
		double fdaughter_truth_Momentum[NMAXDAUGTHERS][4];
		double fdaughter_truth_EndMomentum[NMAXDAUGTHERS][4];
		double fdaughter_truth_Pt[NMAXDAUGTHERS];
		double fdaughter_truth_Mass[NMAXDAUGTHERS];
		double fdaughter_truth_Theta[NMAXDAUGTHERS];
		double fdaughter_truth_Phi[NMAXDAUGTHERS];
		int fdaughter_truth_Process[NMAXDAUGTHERS];

		double minX =  -360.0;
		double maxX = 360.0;
		double minY =0.0;
		double maxY = 600.0;
		double minZ =  0.0;
		double maxZ = 695.0;

		//g4reweight
		int true_beam_PDG;
		int true_beam_ID;
		//double true_beam_len;
		//int true_beam_nElasticScatters;
		std::vector<double> g4rw_primary_plus_sigma_weight;
		std::vector<double> g4rw_primary_minus_sigma_weight;
		std::vector<double> g4rw_primary_weights;
		std::vector<std::string> g4rw_primary_var;
		double g4rw_primary_singular_weight;

		std::vector<double> g4rw_p1; 
		std::vector<double> g4rw_p2; 
		std::vector<double> g4rw_set_weights;
		//std::string fGeneratorTag;


		//Geant4Reweight stuff
		int RW_PDG;
		TFile FracsFile, XSecFile;
		std::vector<fhicl::ParameterSet> ParSet;
		G4ReweightParameterMaker ParMaker;
                G4ReweightManager RWManager;
		G4MultiReweighter MultiRW;
		G4ReweighterFactory RWFactory;
		G4Reweighter * theRW;

		//bool fVerbose;

};


protoana::protonmc::protonmc(fhicl::ParameterSet const & p)
	:
		EDAnalyzer(p),
		dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
		//fNNetModuleLabel(p.get<art::InputTag>("NNetModuleLabel")),
		fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
		fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
		fTrackerTag(p.get<std::string>("TrackerTag")),
		fShowerTag(p.get<std::string>("ShowerTag")),
		fPFParticleTag(p.get<std::string>("PFParticleTag")),
		fGeneratorTag(p.get<std::string>("GeneratorTag")),
		//fMuMCSInputTag(p.get<std::string>("MuMCSInputTag")),
		//fMCSFitter(p.get<fhicl::ParameterSet>("TrajMCSFitter")),


		RW_PDG(p.get<int>("RW_PDG")),
		FracsFile( (p.get< std::string >( "FracsFile" )).c_str(), "OPEN" ),
		XSecFile( (p.get< std::string >( "XSecFile" )).c_str(), "OPEN"),
		ParSet(p.get<std::vector<fhicl::ParameterSet>>("ParameterSet")),
		ParMaker(ParSet, RW_PDG),
                RWManager({p.get<fhicl::ParameterSet>("Material")}),
		MultiRW(RW_PDG, FracsFile, ParSet,
                        p.get<fhicl::ParameterSet>("Material"),
                        &RWManager) {
			theRW = RWFactory.BuildReweighter(RW_PDG, &FracsFile,
					ParMaker.GetFSHists(),
                                        p.get<fhicl::ParameterSet>("Material"),
                                        &RWManager,
					ParMaker.GetElasticHist()/*, true*/ );
		}

//fVerbose(p.get<bool>("Verbose")) {
//}





void protoana::protonmc::beginJob(){

	art::ServiceHandle<art::TFileService> tfs;

	fPandoraBeam = tfs->make<TTree>("PandoraBeam", "Beam events reconstructed with Pandora");
	fPandoraBeam->Branch("run",                           &fRun,                          "run/I");
	fPandoraBeam->Branch("subrun",                        &fSubRun,                       "subrun/I");
	fPandoraBeam->Branch("event",                         &fevent,                        "event/I");
	fPandoraBeam->Branch("timestamp",                     &fTimeStamp,                    "timestamp/D");
	fPandoraBeam->Branch("Nactivefembs",                  &fNactivefembs,                 "Nactivefembs[5]/I");

	fPandoraBeam->Branch("beamtrigger",                   &fbeamtrigger,                  "beamtrigger/I");
	fPandoraBeam->Branch("tof",                           &ftof,                          "tof/D");
	fPandoraBeam->Branch("cerenkov1",                     &fcerenkov1,                    "cerenkov1/I");
	fPandoraBeam->Branch("cerenkov2",                     &fcerenkov2,                    "cerenkov2/I");
	fPandoraBeam->Branch("beamtrackMomentum",             &fbeamtrackMomentum,            "beamtrackMomentum/D");
	fPandoraBeam->Branch("beamtrackP",                    &fbeamtrackP,                   "beamtrackP[3]/D");
	fPandoraBeam->Branch("beamtrackEnergy",               &fbeamtrackEnergy,              "beamtrackEnergy/D");
	fPandoraBeam->Branch("beamtrackPos",                  &fbeamtrackPos,                 "beamtrackPos[3]/D");
	fPandoraBeam->Branch("beamtrackDir",                  &fbeamtrackDir,                 "beamtrackDir[3]/D");
	fPandoraBeam->Branch("beamtrackTime",                 &fbeamtrackTime,                "beamtrackTime/D");
	fPandoraBeam->Branch("beamtrackPdg",                  &fbeamtrackPdg,                 "beamtrackPdg/I");
	fPandoraBeam->Branch("beamtrackID",                   &fbeamtrackID,                  "beamtrackID/I");
	//fPandoraBeam->Branch("NumberBeamTrajectoryPoints",    &NumberBeamTrajectoryPoints,    "NumberBeamTrajectoryPoints/I");
	fPandoraBeam->Branch("beamtrk_x",&beamtrk_x);
	fPandoraBeam->Branch("beamtrk_y",&beamtrk_y);
	fPandoraBeam->Branch("beamtrk_z",&beamtrk_z);
	fPandoraBeam->Branch("beamtrk_Px",&beamtrk_Px);
	fPandoraBeam->Branch("beamtrk_Py",&beamtrk_Py);
	fPandoraBeam->Branch("beamtrk_Pz",&beamtrk_Pz);
	fPandoraBeam->Branch("beamtrk_Eng",&beamtrk_Eng);

  	fPandoraBeam->Branch("beamMomentum_spec",&beamMomentum_spec);
  	fPandoraBeam->Branch("beamPosx_spec",&beamPosx_spec);
  	fPandoraBeam->Branch("beamPosy_spec",&beamPosy_spec);
  	fPandoraBeam->Branch("beamPosz_spec",&beamPosz_spec);
  	fPandoraBeam->Branch("beamDirx_spec",&beamDirx_spec);
  	fPandoraBeam->Branch("beamDiry_spec",&beamDiry_spec);
  	fPandoraBeam->Branch("beamDirz_spec",&beamDirz_spec);


	fPandoraBeam->Branch("ke_ff", &ke_ff, "ke_ff/D");
	fPandoraBeam->Branch("pos_ff",&pos_ff);

	fPandoraBeam->Branch("primtrk_ke_true",&primtrk_ke_true);
	fPandoraBeam->Branch("primtrk_ke_reco",&primtrk_ke_reco);
	fPandoraBeam->Branch("primtrk_end_g4process",&primtrk_end_g4process);
	fPandoraBeam->Branch("primtrk_range_true",&primtrk_range_true);

	fPandoraBeam->Branch("ke_ff_true", &ke_ff_true, "ke_ff_true/D");
	//fPandoraBeam->Branch("primtrk_hit_processkey",&primtrk_hit_processkey);
	fPandoraBeam->Branch("primtrk_hitx_true",&primtrk_hitx_true);
	fPandoraBeam->Branch("primtrk_hity_true",&primtrk_hity_true);
	fPandoraBeam->Branch("primtrk_hitz_true",&primtrk_hitz_true);
	fPandoraBeam->Branch("primtrk_trkid_true",&primtrk_trkid_true);
	fPandoraBeam->Branch("primtrk_edept_true",&primtrk_edept_true);

	fPandoraBeam->Branch("pfphit_peaktime_c",&pfphit_peaktime_c);
	fPandoraBeam->Branch("pfphit_peaktime_u",&pfphit_peaktime_u);
	fPandoraBeam->Branch("pfphit_peaktime_v",&pfphit_peaktime_v);
	fPandoraBeam->Branch("pfphit_wireid_c",&pfphit_wireid_c);
	fPandoraBeam->Branch("pfphit_wireid_u",&pfphit_wireid_u);
	fPandoraBeam->Branch("pfphit_wireid_v",&pfphit_wireid_v);

	fPandoraBeam->Branch("primtrk_true_x",&primtrk_true_x);
	fPandoraBeam->Branch("primtrk_true_y",&primtrk_true_y);
	fPandoraBeam->Branch("primtrk_true_z",&primtrk_true_z);
	fPandoraBeam->Branch("primtrk_true_trkid",&primtrk_true_trkid);
	fPandoraBeam->Branch("primtrk_true_edept",&primtrk_true_edept);
	fPandoraBeam->Branch("primtrk_true_wid",&primtrk_true_wid);
	//fPandoraBeam->Branch("primtrk_true_tpc",&primtrk_true_tpc);
	//fPandoraBeam->Branch("primtrk_true_plane",&primtrk_true_plane);

	fPandoraBeam->Branch("primtrj_true_x",&primtrj_true_x);
	fPandoraBeam->Branch("primtrj_true_y",&primtrj_true_y);
	fPandoraBeam->Branch("primtrj_true_z",&primtrj_true_z);
	fPandoraBeam->Branch("primtrj_true_edept",&primtrj_true_edept);

	//fPandoraBeam->Branch("primtrk_mcs_angles_reco",&primtrk_mcs_angles_reco);
	//fPandoraBeam->Branch("primtrk_mcs_radlengths_reco",&primtrk_mcs_radlengths_reco);

	fPandoraBeam->Branch("vertex",                        &fvertex,                       "vertex[3]/D");
	fPandoraBeam->Branch("secvertex",                     &fsecvertex,                    "secvertex[3]/D");
	fPandoraBeam->Branch("isprimarytrack",                &fisprimarytrack,               "isprimarytrack/I");
	fPandoraBeam->Branch("isprimaryshower",               &fisprimaryshower,              "isprimaryshower/I");
	fPandoraBeam->Branch("primaryBDTScore",               &fprimaryBDTScore,              "primaryBDTScore/D");
	fPandoraBeam->Branch("primaryNHits",                  &fprimaryNHits,                 "fprimaryNHits/I");
	fPandoraBeam->Branch("primaryTheta",                  &fprimaryTheta,                 "primaryTheta/D");
	fPandoraBeam->Branch("primaryPhi",                    &fprimaryPhi,                   "primaryPhi/D");
	fPandoraBeam->Branch("primaryLength",                 &fprimaryLength,                "primaryLength/D");
	fPandoraBeam->Branch("primaryMomentum",               &fprimaryMomentum,              "primaryMomentum/D");
	fPandoraBeam->Branch("primaryEndMomentum",            &fprimaryEndMomentum,           "primaryEndMomentum/D");
	fPandoraBeam->Branch("primaryEndPosition",            &fprimaryEndPosition,           "primaryEndPosition[3]/D");
	fPandoraBeam->Branch("primaryStartPosition",          &fprimaryStartPosition,         "primaryStartPosition[3]/D");
	fPandoraBeam->Branch("primaryEndDirection",           &fprimaryEndDirection,          "primaryEndDirection[3]/D");
	fPandoraBeam->Branch("primaryStartDirection",         &fprimaryStartDirection,        "primaryStartDirection[3]/D");
	fPandoraBeam->Branch("primaryOpeningAngle",           &fprimaryOpeningAngle,          "primaryOpeningAngle/D");
	fPandoraBeam->Branch("primaryID",                     &fprimaryID,                    "primaryID/I");
	fPandoraBeam->Branch("primaryShowerBestPlane",        &fprimaryShowerBestPlane,       "primaryShowerBestPlane/I");
	fPandoraBeam->Branch("primaryShowerEnergy",           &fprimaryShowerEnergy,          "primaryShowerEnergy/I");
	fPandoraBeam->Branch("primaryShowerMIPEnergy",        &fprimaryShowerMIPEnergy,       "primaryShowerMIPEnergy/I");
	fPandoraBeam->Branch("primaryShowerdEdx",             &fprimaryShowerdEdx,            "primaryShowerdEdx/I");
	fPandoraBeam->Branch("primaryMomentumByRangeProton",  &fprimaryMomentumByRangeProton, "primaryMomentumByRangeProton/D");
	fPandoraBeam->Branch("primaryMomentumByRangeMuon",    &fprimaryMomentumByRangeMuon,   "primaryMomentumByRangeMuon/D");
	fPandoraBeam->Branch("primaryKineticEnergy",          &fprimaryKineticEnergy,         "primaryKineticEnergy[3]/D");
	fPandoraBeam->Branch("primaryRange",                  &fprimaryRange,                 "primaryRange[3]/D");
	fPandoraBeam->Branch("primaryT0",                     &fprimaryT0,                    "primaryT0/D");

	fPandoraBeam->Branch("primary_truth_TrackId",         &fprimary_truth_TrackId,         "primary_truth_TrackId/I");
	fPandoraBeam->Branch("primary_truth_Pdg",             &fprimary_truth_Pdg,             "primary_truth_Pdg/I");
	fPandoraBeam->Branch("truthpdg",                      &ftruthpdg,                      "truthpdg/I");
	fPandoraBeam->Branch("primary_truth_StartPosition",   &fprimary_truth_StartPosition,   "primary_truth_StartPosition[4]/D");
	fPandoraBeam->Branch("primary_truth_StartPosition_MC",   &fprimary_truth_StartPosition_MC,   "primary_truth_StartPosition_MC[4]/D");
	fPandoraBeam->Branch("primary_truth_EndPosition",     &fprimary_truth_EndPosition,     "primary_truth_EndPosition[4]/D");
	fPandoraBeam->Branch("primary_truth_EndPosition_MC",     &fprimary_truth_EndPosition_MC,     "primary_truth_EndPosition_MC[4]/D");
	fPandoraBeam->Branch("primary_truth_Momentum",        &fprimary_truth_Momentum,        "primary_truth_Momentum[4]/D");
	fPandoraBeam->Branch("primary_truth_EndMomentum",     &fprimary_truth_EndMomentum,     "primary_truth_EndMomentum[4]/D");
	fPandoraBeam->Branch("primary_truth_P",               &fprimary_truth_P,               "primary_truth_P/D");
	fPandoraBeam->Branch("primary_truth_Pt",              &fprimary_truth_Pt,              "primary_truth_Pt/D");
	fPandoraBeam->Branch("primary_truth_Mass",            &fprimary_truth_Mass,            "primary_truth_Mass/D");
	fPandoraBeam->Branch("primary_truth_Theta",           &fprimary_truth_Theta,           "primary_truth_Theta/D");
	fPandoraBeam->Branch("primary_truth_Phi",             &fprimary_truth_Phi,             "primary_truth_Phi/D");
	//fPandoraBeam->Branch("primary_truth_Process",         &fprimary_truth_Process,         "primary_truth_Process/I");
	//fPandoraBeam->Branch("primary_truth_Process",         &fprimary_truth_Process);
	fPandoraBeam->Branch("primary_truth_EndProcess",      &fprimary_truth_EndProcess);
	//fPandoraBeam->Branch("primary_backtrker_truth_Process",         &fprimary_backtrker_truth_Process);
	//fPandoraBeam->Branch("primary_backtrker_truth_EndProcess",      &fprimary_backtrker_truth_EndProcess);
	fPandoraBeam->Branch("primary_truth_Isbeammatched",   &fprimary_truth_Isbeammatched,   "primary_truth_Isbeammatched/I");
	fPandoraBeam->Branch("primary_truth_NDaughters",      &fprimary_truth_NDaughters,      "primary_truth_NDaughters/I");
	//fPandoraBeam->Branch("G4Process_Primary_End",      &G4Process_Primary_End);

	fPandoraBeam->Branch("interactionX",&interactionX);
	fPandoraBeam->Branch("interactionY",&interactionY);
	fPandoraBeam->Branch("interactionZ",&interactionZ);
	fPandoraBeam->Branch("interactionE",&interactionE);
	fPandoraBeam->Branch("interactionProcesslist",&interactionProcesslist);
	fPandoraBeam->Branch("interactionAngles",&interactionAngles);

	fPandoraBeam->Branch("Zintersection",&Zintersection);
	fPandoraBeam->Branch("Yintersection",&Yintersection);
	fPandoraBeam->Branch("Xintersection",&Xintersection);
	fPandoraBeam->Branch("timeintersection",&timeintersection);

	fPandoraBeam->Branch("interaction_wid_c",&interaction_wid_c);
	fPandoraBeam->Branch("interaction_tt_c",&interaction_tt_c);
	fPandoraBeam->Branch("interaction_wid_v",&interaction_wid_v);
	fPandoraBeam->Branch("interaction_tt_v",&interaction_tt_v);
	fPandoraBeam->Branch("interaction_wid_u",&interaction_wid_u);
	fPandoraBeam->Branch("interaction_tt_u",&interaction_tt_u);

	//fPandoraBeam->Branch("inelscore_c",&inelscore_c);
	//fPandoraBeam->Branch("elscore_c",&elscore_c);
	//fPandoraBeam->Branch("nonescore_c",&nonescore_c);
	fPandoraBeam->Branch("x_c",&x_c);
	fPandoraBeam->Branch("y_c",&y_c);
	fPandoraBeam->Branch("z_c",&z_c);
	//fPandoraBeam->Branch("xyz_inelscore_c",&xyz_inelscore_c,"xyz_inelscore_c[3]/D");
	//fPandoraBeam->Branch("xyz_elscore_c",&xyz_elscore_c,"xyz_elscore_c[3]/D");

	fPandoraBeam->Branch("wid_c",&wid_c);
	fPandoraBeam->Branch("tt_c",&tt_c);
	fPandoraBeam->Branch("ch_c",&ch_c);
	fPandoraBeam->Branch("wirez_c",&wirez_c);
	fPandoraBeam->Branch("pt_c",&pt_c);
	fPandoraBeam->Branch("q_c",&q_c);
	fPandoraBeam->Branch("a_c",&a_c);
	fPandoraBeam->Branch("wireno_c",&wireno_c);

	fPandoraBeam->Branch("wid_v",&wid_v);
	fPandoraBeam->Branch("tt_v",&tt_v);
	fPandoraBeam->Branch("ch_v",&ch_v);
	fPandoraBeam->Branch("wirez_v",&wirez_v);
	fPandoraBeam->Branch("pt_v",&pt_v);
	fPandoraBeam->Branch("q_v",&q_v);
	fPandoraBeam->Branch("a_v",&a_v);
	fPandoraBeam->Branch("wireno_v",&wireno_v);

	fPandoraBeam->Branch("wid_u",&wid_u);
	fPandoraBeam->Branch("tt_u",&tt_u);
	fPandoraBeam->Branch("ch_u",&ch_u);
	fPandoraBeam->Branch("wirez_u",&wirez_u);
	fPandoraBeam->Branch("pt_u",&pt_u);
	fPandoraBeam->Branch("q_u",&q_u);
	fPandoraBeam->Branch("a_u",&a_u);
	fPandoraBeam->Branch("wireno_u",&wireno_u);



	fPandoraBeam->Branch("NDAUGHTERS",                    &fNDAUGHTERS,                   "NDAUGHTERS/I");
	fPandoraBeam->Branch("isdaughtertrack",               &fisdaughtertrack,              "isdaughtertrack[NDAUGHTERS]/I");
	fPandoraBeam->Branch("isdaughtershower",              &fisdaughtershower,             "isdaughtershower[NDAUGHTERS]/I");
	fPandoraBeam->Branch("daughterNHits",                 &fdaughterNHits,                "daughterNHits[NDAUGHTERS]/I");
	fPandoraBeam->Branch("daughterTheta",                 &fdaughterTheta,                "daughterTheta[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterPhi",                   &fdaughterPhi,                  "daughterPhi[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterLength",                &fdaughterLength,               "daughterLength[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterMomentum",              &fdaughterMomentum,             "daughterMomentum[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterEndMomentum",           &fdaughterEndMomentum,          "daughterEndMomentum[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterEndPosition",           &fdaughterEndPosition,          "daughterEndPosition[NDAUGHTERS][3]/D");
	fPandoraBeam->Branch("daughterStartPosition",         &fdaughterStartPosition,        "daughterStartPosition[NDAUGHTERS][3]/D");
	fPandoraBeam->Branch("daughterStartDirection",        &fdaughterStartDirection,       "daughterStartDirection[NDAUGHTERS][3]/D");
	fPandoraBeam->Branch("daughterEndDirection",          &fdaughterEndDirection,         "daughterEndDirection[NDAUGHTERS][3]/D");
	fPandoraBeam->Branch("daughterOpeningAngle",          &fdaughterOpeningAngle,         "daughterOpeningAngle[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterShowerBestPlane",       &fdaughterShowerBestPlane,      "daughterShowerBestPlane[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterShowerEnergy",          &fdaughterShowerEnergy,         "daughterShowerEnergy[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterShowerMIPEnergy",       &fdaughterShowerMIPEnergy,      "daughterShowerMIPEnergy[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterShowerdEdx",            &fdaughterShowerdEdx,           "daughterShowerdEdx[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterMomentumByRangeProton", &fdaughterMomentumByRangeProton,"daughterMomentumByRangeProton[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterMomentumByRangeMuon",   &fdaughterMomentumByRangeMuon,  "daughterMomentumByRangeMuon[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughterKineticEnergy",         &fdaughterKineticEnergy,        "daughterKineticEnergy[NDAUGHTERS][3]/D");
	fPandoraBeam->Branch("daughterRange",                 &fdaughterRange,                "daughterRange[NDAUGHTERS][3]/D");
	fPandoraBeam->Branch("daughterID",                    &fdaughterID,                   "daughterID[NDAUGHTERS]/I");
	//fPandoraBeam->Branch("daughterT0",                    &fdaughterT0,                   "daughterT0[NDAUGHTERS]/D");

	fPandoraBeam->Branch("daughter_truth_TrackId",        &fdaughter_truth_TrackId,       "daughter_truth_TrackId[NDAUGHTERS]/I");
	fPandoraBeam->Branch("daughter_truth_Pdg",            &fdaughter_truth_Pdg,           "daughter_truth_Pdg[NDAUGHTERS]/I");
	fPandoraBeam->Branch("daughter_truth_StartPosition",  &fdaughter_truth_StartPosition, "daughter_truth_StartPosition[NDAUGHTERS][4]/D");
	fPandoraBeam->Branch("daughter_truth_EndPosition",    &fdaughter_truth_EndPosition,   "daughter_truth_EndPosition[NDAUGHTERS][4]/D");
	fPandoraBeam->Branch("daughter_truth_Momentum",       &fdaughter_truth_Momentum,      "daughter_truth_Momentum[NDAUGHTERS][4]/D");
	fPandoraBeam->Branch("daughter_truth_EndMomentum",    &fdaughter_truth_EndMomentum,   "daughter_truth_EndMomentum[NDAUGHTERS][4]/D");
	fPandoraBeam->Branch("daughter_truth_P",              &fdaughter_truth_P,             "daughter_truth_P[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughter_truth_Pt",             &fdaughter_truth_Pt,            "daughter_truth_Pt[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughter_truth_Mass",           &fdaughter_truth_Mass,          "daughter_truth_Mass[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughter_truth_Theta",          &fdaughter_truth_Theta,         "daughter_truth_Theta[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughter_truth_Phi",            &fdaughter_truth_Phi,           "daughter_truth_Phi[NDAUGHTERS]/D");
	fPandoraBeam->Branch("daughter_truth_Process",        &fdaughter_truth_Process,       "daughter_truth_Process[NDAUGHTERS]/I");

	fPandoraBeam->Branch("primtrk_dqdx",&primtrk_dqdx);
	fPandoraBeam->Branch("primtrk_dedx",&primtrk_dedx);
	fPandoraBeam->Branch("primtrk_resrange",&primtrk_resrange);
	fPandoraBeam->Branch("primtrk_range",&primtrk_range);
	fPandoraBeam->Branch("primtrk_hitx",&primtrk_hitx);
	fPandoraBeam->Branch("primtrk_hity",&primtrk_hity);
	fPandoraBeam->Branch("primtrk_hitz",&primtrk_hitz);
	fPandoraBeam->Branch("primtrk_pitch",&primtrk_pitch);
	fPandoraBeam->Branch("primtrk_wid",&primtrk_wid);

	//fPandoraCosmics = tfs->make<TTree>("PandoraCosmics", "Cosmic tracks reconstructed with Pandora");
	//fPandoraCosmics->Branch("run",                 &fRun,                "run/I");
	//fPandoraCosmics->Branch("subrun",              &fSubRun,             "subrun/I");
	//fPandoraCosmics->Branch("event",               &fevent,              "event/I");
	//fPandoraCosmics->Branch("timestamp",           &fTimeStamp,          "timestamp/D");
	//fPandoraCosmics->Branch("Nactivefembs",        &fNactivefembs,       "Nactivefembs[5]/I");
	//fPandoraCosmics->Branch("beamtrigger",         &fbeamtrigger,        "beamtrigger/I");
	//fPandoraCosmics->Branch("tof",                 &ftof,                "tof/D");
	//fPandoraCosmics->Branch("cerenkov1",           &fcerenkov1,          "cerenkov1/I");
	//fPandoraCosmics->Branch("cerenkov2",           &fcerenkov2,          "cerenkov2/I");
	//fPandoraCosmics->Branch("beamtrackMomentum",   &fbeamtrackMomentum,  "beamtrackMomentum/D");
	//fPandoraCosmics->Branch("beamtrackPos",        &fbeamtrackPos,       "beamtrackPos[3]/D");
	//fPandoraCosmics->Branch("beamtrackDir",        &fbeamtrackDir,       "beamtrackDir[3]/D");


  	fPandoraBeam->Branch("g4rw_primary_weights", &g4rw_primary_weights);
  	fPandoraBeam->Branch("g4rw_primary_singular_weight", &g4rw_primary_singular_weight);
  	fPandoraBeam->Branch("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
  	fPandoraBeam->Branch("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  	fPandoraBeam->Branch("g4rw_primary_var", &g4rw_primary_var);
  	fPandoraBeam->Branch("g4rw_p1", &g4rw_p1);
  	fPandoraBeam->Branch("g4rw_p2", &g4rw_p2);
  	fPandoraBeam->Branch("g4rw_set_weights", &g4rw_set_weights);
	

}

void protoana::protonmc::analyze(art::Event const & evt){

	// Initialise tree parameters
	Initialise();

	int beamid=-9999;
	int truthid=-999;  

	art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
	fGeometry = &*(art::ServiceHandle<geo::Geometry>());
	//anab::MVAReader<recob::Hit,3> hitResults(evt, fNNetModuleLabel);

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
	//T0
	std::vector<const anab::T0*> T0s;    
	nT0s = T0s.size();
	std::cout << "Got " << nT0s << " T0s" <<std::endl;

	if(nT0s > 0){
		std::cout << "T0s size: " << nT0s << std::endl;          
		MCTruthT0 = T0s[0]->Time();        
	}
	else{
		std::cout << "No T0s found" << std::endl;
		MCTruthT0 = 0;
	}        
        TickT0 = MCTruthT0 / sampling_rate(clockData);
	std::cout<<"TickT0:"<<TickT0<<std::endl;


	fRun = evt.run();
	fSubRun = evt.subRun();
	fevent  = evt.id().event(); 
	art::Timestamp ts = evt.time();
	if (ts.timeHigh() == 0){
		TTimeStamp ts2(ts.timeLow());
		fTimeStamp = ts2.AsDouble();
	}
	else{
		TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
		fTimeStamp = ts2.AsDouble();
	}

	// Get number of active fembs
	if(!evt.isRealData()){
		for(int k=0; k < 6; k++)
			fNactivefembs[k] = 20;
	}
	else{
		for(int k=0; k < 6; k++)
			fNactivefembs[k] = dataUtil.GetNActiveFembsForAPA(evt, k);
	}

	bool beamTriggerEvent = false;
	// If this event is MC then we can check what the true beam particle is
	auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);  
	if(!evt.isRealData()){
		//for prod. 3, new implementation to access the beam momentum from spectrometer //////////////////////////////
  		auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("generator");
  		std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  		if( beamHandle.isValid()){
    			art::fill_ptr_vector(beamVec, beamHandle);
  		}
  		const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one

  		//Access momentum
  		const std::vector< double > & momenta = beamEvent.GetRecoBeamMomenta();
  		std::cout << "Number of reconstructed beam momenta from spec: " << momenta.size() << std::endl;

  		if( momenta.size() > 0 ) std::cout << "Measured Beam Momentum from spec: " << momenta.at(0) << std::endl;
		
		//std::cout<<"beam mom size:"<<momenta.size()<<std::endl;
        	for (size_t i = 0; i<momenta.size(); ++i){
          		beamMomentum_spec.push_back(momenta[i]);
	  		//std::cout<<"beam mom["<<i<<"]:"<<momenta[i]<<" [GeV]"<<std::endl;
        	}

		auto & btracks = beamEvent.GetBeamTracks();
		std::cout<<"beam trk size:"<<btracks.size()<<std::endl;
        	for (size_t i = 0; i<btracks.size(); ++i){
	  		std::cout<<"beamPosx/beamPosy/beamPosz:"<<btracks[i].End().X()<<"/"<<btracks[i].End().Y()<<"/"<<btracks[i].End().Z()<<std::endl;
	  		std::cout<<"beamDirx/beamDiry/beamDirz:"<<btracks[i].StartDirection().X()<<"/"<<btracks[i].StartDirection().Y()<<"/"<<btracks[i].StartDirection().Z()<<std::endl;

          		beamPosx_spec.push_back(btracks[i].End().X());
          		beamPosy_spec.push_back(btracks[i].End().Y());
          		beamPosz_spec.push_back(btracks[i].End().Z());
          		beamDirx_spec.push_back(btracks[i].StartDirection().X());
          		beamDiry_spec.push_back(btracks[i].StartDirection().Y());
          		beamDirz_spec.push_back(btracks[i].StartDirection().Z());

		}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////



		// Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
		// simulation has fGeneratorTag = "generator"
		auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);

		// Also get the reconstructed beam information in the MC - TO DO
		//auto beamsim = evt.getValidHandle<std::vector<sim::ProtoDUNEbeamsim> >(fGeneratorTag);
		//const sim::ProtoDUNEbeamsim beamsimobj = (*beamsim)[0];
		//std::cout << beamsimobj.NInstruments() << std::endl;
		//sim::ProtoDUNEBeamInstrument beamsim_tof1 = beamsimobj.GetInstrument("TOF1");
		//sim::ProtoDUNEBeamInstrument beamsim_trig2 = beamsimobj.GetInstrument("TRIG2");
		//std::cout << beamsim_trig2.GetT() - beamsim_tof1.GetT() << " , " << beamsim_trig2.GetSmearedVar1() - beamsim_tof1.GetSmearedVar1() << std::endl;
		//std::cout << beamsimobj.GetInstrument("TRIG2").GetT() - beamsimobj.GetInstrument("TOF1").GetT() << " , " << beamsimobj.GetInstrument("TRIG2").GetSmearedVar1() - beamsimobj.GetInstrument("TOF1").GetSmearedVar1() << std::endl;

		// mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
		// of these, so we pass the first element into the function to get the good particle
		const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);



		if(geantGoodParticle != 0x0){
			std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() 
				<< " , track id = " << geantGoodParticle->TrackId()
				<< " , Vx/Vy/Vz = " << geantGoodParticle->Vx() << "/"<< geantGoodParticle->Vy() << "/" << geantGoodParticle->Vz() 
				<< std::endl;

			//NumberBeamTrajectoryPoints=geantGoodParticle->NumberTrajectoryPoints();
			//std::cout<<"NumberBeamTrajectoryPoints for beam particles:"<<NumberBeamTrajectoryPoints<<std::endl;

			beamTriggerEvent = true;
			fbeamtrigger       = 12;
			fbeamtrackPos[0]   = geantGoodParticle->Vx();
			fbeamtrackPos[1]   = geantGoodParticle->Vy();
			fbeamtrackPos[2]   = geantGoodParticle->Vz();
			fbeamtrackMomentum = geantGoodParticle->P();
			fbeamtrackP[0]     = geantGoodParticle->Px();
			fbeamtrackP[1]     = geantGoodParticle->Py();
			fbeamtrackP[2]     = geantGoodParticle->Pz();
			fbeamtrackEnergy   = geantGoodParticle->E();
			fbeamtrackPdg      = geantGoodParticle->PdgCode();
			fbeamtrackTime     = geantGoodParticle->T();
			fbeamtrackID       = geantGoodParticle->TrackId();
			//prim_energy=0;

			for(size_t i_s=0; i_s < geantGoodParticle->NumberTrajectoryPoints(); i_s++){ //loop over beam tracks
				beamtrk_x.push_back(geantGoodParticle->Position(i_s).X());
				beamtrk_y.push_back(geantGoodParticle->Position(i_s).Y());
				beamtrk_z.push_back(geantGoodParticle->Position(i_s).Z());

				beamtrk_Px.push_back(geantGoodParticle->Momentum(i_s).X());
				beamtrk_Py.push_back(geantGoodParticle->Momentum(i_s).Y());
				beamtrk_Pz.push_back(geantGoodParticle->Momentum(i_s).Z());

				beamtrk_Eng.push_back(geantGoodParticle->Momentum(i_s).E()-geantGoodParticle->Mass());
			} //loop over beam trks

			//Get KE at front face of TPC --------------------------------------------//
			int key_reach_tpc=-99;
			for (size_t kk=0; kk<beamtrk_z.size(); ++kk) {  //loop over all beam hits
				double zpos_beam=beamtrk_z.at(kk);
				if ((zpos_beam+0.49375)<0.01&&(zpos_beam+0.49375)>-0.01) { //find key at ff
					key_reach_tpc=kk;
				} //kind key at ff
			} //loop over all beam hits  

			//Define KE and Momentum at the entering point of TPC
			if (key_reach_tpc!=-99) { //ke and pos at front face
				ke_ff=1000.*(beamtrk_Eng.at(key_reach_tpc)); //unit: MeV
				pos_ff.push_back(beamtrk_x.at(key_reach_tpc));
				pos_ff.push_back(beamtrk_y.at(key_reach_tpc));
				pos_ff.push_back(beamtrk_z.at(key_reach_tpc));

				fprimary_truth_Momentum[0]=beamtrk_Px.at(key_reach_tpc);
				fprimary_truth_Momentum[1]=beamtrk_Py.at(key_reach_tpc);
				fprimary_truth_Momentum[2]=beamtrk_Pz.at(key_reach_tpc);
			} //ke and pos at front face
			else {
				std::cout<<"This particle doesn't enter TPC!!!"<<std::endl;
			}	
			//Get KE at front face of TPC --------------------------------------------//

			//Get Truth info
			fprimary_truth_TrackId          = geantGoodParticle->TrackId();
			fprimary_truth_Pdg              = geantGoodParticle->PdgCode();
			beamid                          = geantGoodParticle->TrackId();
			fprimary_truth_StartPosition[3] = geantGoodParticle->T();
			fprimary_truth_EndPosition[3]   = geantGoodParticle->EndT();

			fprimary_truth_StartPosition_MC[0] = geantGoodParticle->Vx();
			fprimary_truth_StartPosition_MC[1] = geantGoodParticle->Vy();
			fprimary_truth_StartPosition_MC[2] = geantGoodParticle->Vz();
			fprimary_truth_StartPosition_MC[3] = geantGoodParticle->T();

			fprimary_truth_EndPosition_MC[0]   = geantGoodParticle->EndX();
			fprimary_truth_EndPosition_MC[1]   = geantGoodParticle->EndY();
			fprimary_truth_EndPosition_MC[2]   = geantGoodParticle->EndZ();
			fprimary_truth_EndPosition_MC[3]   = geantGoodParticle->EndT();

			fprimary_truth_P                = geantGoodParticle->P();
			//fprimary_truth_Momentum[0]      = geantGoodParticle->Px();
			//fprimary_truth_Momentum[1]      = geantGoodParticle->Py();
			//fprimary_truth_Momentum[2]      = geantGoodParticle->Pz();
			fprimary_truth_Momentum[3]      = geantGoodParticle->E();
			fprimary_truth_Pt               = geantGoodParticle->Pt();
			fprimary_truth_Mass             = geantGoodParticle->Mass();
			//fprimary_truth_EndMomentum[0]   = geantGoodParticle->EndPx();
			//fprimary_truth_EndMomentum[1]   = geantGoodParticle->EndPy();
			//fprimary_truth_EndMomentum[2]   = geantGoodParticle->EndPz();
			fprimary_truth_EndMomentum[3]   = geantGoodParticle->EndE();
			fprimary_truth_Theta            = geantGoodParticle->Momentum().Theta();
			fprimary_truth_Phi              = geantGoodParticle->Momentum().Phi();
			fprimary_truth_NDaughters       = geantGoodParticle->NumberDaughters();
			//fprimary_truth_Process          = int(geantGoodParticle->Trajectory().ProcessToKey(geantGoodParticle->Process()));
			//fprimary_truth_Process           = geantGoodParticle->Process(); //HY::wired result
			//fprimary_truth_EndProcess           = geantGoodParticle->EndProcess(); //HY:wierd result

			//fprimary_backtrker_truth_Process          = int(geantGoodParticle->Trajectory().ProcessToKey(geantGoodParticle->Process()));

			//g4 reweight -------------------------------------------------------------------------------------------------//
			true_beam_PDG = geantGoodParticle->PdgCode();
			true_beam_ID = geantGoodParticle->TrackId();
			//true_beam_len = geantGoodParticle->Trajectory().TotalLength();

			std::cout << "Doing reweight" << std::endl;
			if (true_beam_PDG == RW_PDG) { //if PDG=RW_PDG
				std::cout<<"Got particle with PDG="<<RW_PDG<<std::endl;
				G4ReweightTraj theTraj(true_beam_ID, true_beam_PDG, 0, fevent, {0,0});
				bool created = CreateRWTraj(*geantGoodParticle, pi_serv->ParticleList(),
						*fGeometry, fevent, &theTraj);
				if (created && theTraj.GetNSteps()) {

					g4rw_primary_singular_weight = MultiRW.GetWeightFromNominal(theTraj);
					//the following method achieves the same result
					//g4rw_primary_singular_weight = theRW->GetWeight(&theTraj);

					std::vector<double> weights_vec = MultiRW.GetWeightFromAll1DThrows(
							theTraj);
					g4rw_primary_weights.insert(g4rw_primary_weights.end(),
							weights_vec.begin(), weights_vec.end());

					for (size_t i = 0; i < ParSet.size(); ++i) {
						std::pair<double, double> pm_weights =
							MultiRW.GetPlusMinusSigmaParWeight(theTraj, i);

						g4rw_primary_plus_sigma_weight.push_back(pm_weights.first);
						g4rw_primary_minus_sigma_weight.push_back(pm_weights.second);
						g4rw_primary_var.push_back(ParSet[i].get<std::string>("Name"));
					}

					//--------- rw parameters using multi-dimensional array -----------------------------------------------------------------------------------------------//
					if (ParSet.size()==2) { //if parset size=2, i.e. 2x2 array
						for (size_t i = 0; i < 300; ++i) {
							for (size_t j = 0; j < 200; ++j) {
								double tmp_p1=(.9 + i*.001); //reac
								double tmp_p2=(.9 + j*.001); //elast
								//std::cout<<"(p1,p2):("<<tmp_p1<<","<<tmp_p2<<")"<<std::endl; 
								std::vector<double> input_values = {tmp_p1,tmp_p2}; //HY:the index here based on the sequence of the fichl setting
								
								bool set_values = MultiRW.SetAllParameterValues(input_values);
								if (!set_values) continue;

								g4rw_p1.push_back(tmp_p1);
								g4rw_p2.push_back(tmp_p2);
								g4rw_set_weights.push_back(
									MultiRW.GetWeightFromSetParameters(theTraj)
								);

							}
						}
					} //if parset size=2, i.e. 2x2 array
					//--------- rw parameters using multi-dimensional array ----------------------------------------------------------------------------------------------//




				}

			} //if PDG=RW_PDG
			//g4 reweight -------------------------------------------------------------------------------------------------//

		}
	}
	else{
		// For data we can see if this event comes from a beam trigger
		beamTriggerEvent = dataUtil.IsBeamTrigger(evt);

		std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
		auto pdbeamHandle = evt.getHandle< std::vector<beam::ProtoDUNEBeamEvent> >(fBeamModuleLabel);
		if (pdbeamHandle)
			art::fill_ptr_vector(beaminfo, pdbeamHandle);

		for(unsigned int i = 0; i < beaminfo.size(); ++i){
			//if(!beaminfo[i]->CheckIsMatched()) continue;
			fbeamtrigger = beaminfo[i]->GetTimingTrigger();
			fbeamtrackTime = (double)beaminfo[i]->GetRDTimestamp();

			// If ToF is 0-3 there was a good match corresponding to the different pair-wise combinations of the upstream and downstream channels
			if(beaminfo[i]->GetTOFChan() >= 0)
				ftof =  beaminfo[i]->GetTOF();

			// Get Cerenkov
			if(beaminfo[i]->GetBITrigger() == 1){
				fcerenkov1 = beaminfo[i]->GetCKov0Status();
				fcerenkov2 = beaminfo[i]->GetCKov1Status();
			}

			// Beam particle could have more than one tracks - for now take the first one, need to do this properly
			auto & tracks = beaminfo[i]->GetBeamTracks();
			if(!tracks.empty()){
				fbeamtrackPos[0] = tracks[0].End().X();
				fbeamtrackPos[1] = tracks[0].End().Y();
				fbeamtrackPos[2] = tracks[0].End().Z();
				fbeamtrackDir[0] = tracks[0].StartDirection().X();
				fbeamtrackDir[1] = tracks[0].StartDirection().Y();
				fbeamtrackDir[2] = tracks[0].StartDirection().Z();
			}

			// Beam momentum
			auto & beammom = beaminfo[i]->GetRecoBeamMomenta();
			if(!beammom.empty())
				fbeamtrackMomentum = beammom[0];

			// For now only take the first beam particle - need to add some criteria if more than one are found
			break;

		}
	}

	/*
	// Now we want to access the output from Pandora. This comes in the form of particle flow objects (recob::PFParticle).
	// The primary PFParticles are those we want to consider and these PFParticles then have a hierarchy of daughters that
	// describe the whole interaction of a given primary particle
	//
	//                     / daughter track
	//                    /
	//  primary track    /   
	//  ---------------- ---- daughter track
	//                   \
	//                   /\-
	//                   /\\-- daughter shower
	//
	// The above primary PFParticle will have links to three daughter particles, two track-like and one shower-like
	*/

	// Track momentum algorithm calculates momentum based on track range
	trkf::TrackMomentumCalculator trmom;
	//trmom.SetMinLength(100);

	// Get all of the PFParticles, by default from the "pandora" product
	auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
	std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(evt,fPFParticleTag) << std::endl;

	// We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
	// to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
	// are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
	// PFParticles considered as primary
	//std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
	auto pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);

	//cluster information
	std::vector<art::Ptr<recob::Track> > tracklist;
	std::vector<art::Ptr<recob::PFParticle> > pfplist;

	auto trackListHandle = evt.getHandle< std::vector<recob::Track> >("pandoraTrack");
	if (trackListHandle) art::fill_ptr_vector(tracklist, trackListHandle);
	else return;

	auto PFPListHandle = evt.getHandle< std::vector<recob::PFParticle> > ("pandora");
	if (PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);

	std::vector<art::Ptr<recob::Cluster>> clusterlist;
	auto clusterListHandle = evt.getHandle< std::vector<recob::Cluster> >("pandora") ; // to get information about the hits
	if (clusterListHandle) art::fill_ptr_vector(clusterlist, clusterListHandle);

	art::FindManyP<recob::Cluster> fmcp(PFPListHandle,evt,"pandora");
	art::FindManyP<recob::Track> pftrack(PFPListHandle,evt,"pandoraTrack");

	std::cout<<"number of pfp_particles "<<pfplist.size()<<std::endl;
	std::cout<<" size of pfParticles "<<pfParticles.size()<<std::endl;
	art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt,"pandoraTrack"); // to associate tracks and hits

	/* for(size_t p1=0;p1<pfplist.size();p1++){
	   std::vector<art::Ptr<recob::Track>> trk=pftrack.at(p1);
	   if(trk.size()) std::cout<<" trk  key "<<trk[0].key()<<std::endl; 

	   std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(p1);
	   std::cout<<"cluster size for each particle "<<allClusters.size()<<std::endl;
	   for(size_t c1=0;c1<allClusters.size();c1++){
	   std::cout<<" cluster ID "<<allClusters[c1]->ID();
	   std::cout<<" plane number "<<allClusters[c1]->Plane().Plane;
	   std::cout<<" TPC number "<<allClusters[c1]->Plane().TPC;
	   std::cout<<" start wire "<<allClusters[c1]->StartWire();
	   std::cout<<" end wire "<<allClusters[c1]->EndWire();
	   std::cout<<" start tick "<<allClusters[c1]->StartTick();
	   std::cout<<" end tick "<<allClusters[c1]->EndTick();
	   }
	   }*/




	// We can now look at these particles
	for(const recob::PFParticle* particle : pfParticles){

		// Pandora's BDT beam-cosmic score
		fprimaryBDTScore = (double)pfpUtil.GetBeamCosmicScore(*particle,evt,fPFParticleTag);

		// NHits associated with this pfParticle
		fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();

		// Get the T0 for this pfParticle
		std::vector<anab::T0> pfT0vec = pfpUtil.GetPFParticleT0(*particle,evt,fPFParticleTag);
		if(!pfT0vec.empty())
			fprimaryT0 = pfT0vec[0].Time();

		//std::cout << "Pdg Code = " << particle->PdgCode() << std::endl;
		// "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
		// of this particle might be more helpful. These return null pointers if not track-like / shower-like
		const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
		const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);

		trkf::TrackMomentumCalculator trkm{1.};
		recob::MCSFitResult res;

		if(thisTrack != 0x0) { //this track
			// Get the true mc particle
                         const simb::MCParticle* mcparticle0 = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, evt, fTrackerTag);
			std::cout<<"inside the this track loop "<<std::endl;
			if(mcparticle0!=0x0) {
				std::cout<<"fTruth PDG: "<<mcparticle0->PdgCode()<<std::endl;
				ftruthpdg=mcparticle0->PdgCode();

				truthid=mcparticle0->TrackId();
				fprimary_truth_Isbeammatched=0;
				if(beamid==truthid) fprimary_truth_Isbeammatched=1;



			}

			fisprimarytrack               = 1;
			fisprimaryshower              = 0;

			//fprimaryID                    = thisTrack->ParticleId(); //Do NOT use ParticleID, u got nothing
			fprimaryID                    = thisTrack->ID();
			fprimaryTheta                 = thisTrack->Theta();
			fprimaryPhi                   = thisTrack->Phi();
			fprimaryLength                = thisTrack->Length();
			fprimaryMomentum              = thisTrack->StartMomentum();
			fprimaryEndMomentum           = thisTrack->EndMomentum();

			fprimaryEndPosition[0]        = thisTrack->Trajectory().End().X();
			fprimaryEndPosition[1]        = thisTrack->Trajectory().End().Y();
			fprimaryEndPosition[2]        = thisTrack->Trajectory().End().Z();
			fprimaryStartPosition[0]      = thisTrack->Trajectory().Start().X();
			fprimaryStartPosition[1]      = thisTrack->Trajectory().Start().Y();
			fprimaryStartPosition[2]      = thisTrack->Trajectory().Start().Z();
			fprimaryEndDirection[0]       = thisTrack->Trajectory().EndDirection().X();
			fprimaryEndDirection[1]       = thisTrack->Trajectory().EndDirection().Y();
			fprimaryEndDirection[2]       = thisTrack->Trajectory().EndDirection().Z();
			fprimaryStartDirection[0]     = thisTrack->Trajectory().StartDirection().X();
			fprimaryStartDirection[1]     = thisTrack->Trajectory().StartDirection().Y();
			fprimaryStartDirection[2]     = thisTrack->Trajectory().StartDirection().Z();

			fprimaryMomentumByRangeMuon   = trmom.GetTrackMomentum(thisTrack->Length(),13);
			fprimaryMomentumByRangeProton = trmom.GetTrackMomentum(thisTrack->Length(),2212);      

			//if (fRun==22628660&&fevent==876&&fprimaryID==52) std::cout<<"\n\nrun:"<<fRun<<" event:"<<fevent<<" ID:"<<fprimaryID<<" match:"<<fprimary_truth_Isbeammatched<<"\n\n"<<std::endl;

			//Get the clusters/hits associated with the PFParticle -------------------------------------------------------------------------//
			//const std::vector<const recob::Cluster*> thisCluster = pfpUtil.GetPFParticleClusters(*particle,evt,fPFParticleTag);
			//std::cout<<"size of cluster:"<<thisCluster.size()<<std::endl; //should be three (u,v,c)

			//Get the number of clusters associated to the PFParticle
			//unsigned int num_pfp_clusters = pfpUtil.GetNumberPFParticleClusters(*particle, evt, fPFParticleTag);
			//std::cout<<"num_pfp_clusters:"<<num_pfp_clusters<<std::endl;
			//std::cout<<"pfplist.size="<<pfplist.size()<<std::endl; //all the pfp particles in a pool

			//Get the number of hits
			//unsigned int num_hits=pfpUtil.GetNumberPFParticleHits(*particle, evt, fPFParticleTag);	

			//Get the hits associated to the PFParticle
			const std::vector<const recob::Hit*> pfpHits=pfpUtil.GetPFParticleHits(*particle, evt, fPFParticleTag);
			//std::cout<<"size of hit:"<<pfpHits.size()<<std::endl;
			//std::cout<<"num_hits:"<<num_hits<<std::endl;

			for(const recob::Hit* hi : pfpHits){
				//std::cout<<"(w,t,ID,TPC,Cryostat):("<<hi->WireID().Wire<<","<<hi->PeakTime()<<","<<hi->WireID().Plane<<","<<hi->WireID().TPC<<","<<hi->WireID().Cryostat<<")"<<std::endl;
				if (hi->WireID().Plane==2) {
					pfphit_peaktime_c.push_back((double)hi->PeakTime()); 
					pfphit_wireid_c.push_back((double)hi->WireID().Wire);
				}
				if (hi->WireID().Plane==1) {
					pfphit_peaktime_v.push_back((double)hi->PeakTime()); 
					pfphit_wireid_v.push_back((double)hi->WireID().Wire);
				}
				if (hi->WireID().Plane==0) {
					pfphit_peaktime_u.push_back((double)hi->PeakTime()); 
					pfphit_wireid_u.push_back((double)hi->WireID().Wire);
				}
			}
			//-------------------------------------------------------------------------------------------------------------------------------//


			//reco MCS angles ---------------------------------------------------------------------------------------------------//
			//res = fMCSFitter.fitMcs(*thisTrack);
			//std::vector<float> tmp_primtrk_mcs_angles_reco(res.scatterAngles()); //units are mrad
			//std::vector<float> tmp_primtrk_mcs_radlengths_reco(res.segmentRadLengths()); //seglength/X0
			//primtrk_mcs_angles_reco=tmp_primtrk_mcs_angles_reco; 
			//primtrk_mcs_radlengths_reco=tmp_primtrk_mcs_radlengths_reco;

			//std::vector<double> primtrk_mcs_sigma_hl();
			//float sum_primtrk_mcs_sigma_hl=0.; 
			//float count_primtrk_mcs_sigma_hl=0.;
			//MCS_Angles=res.scatterAngles();
			//double mcs_totlength=0.;
			//for (size_t i_mcs=0; i_mcs<primtrk_mcs_angles_reco.size(); ++i_mcs) {
			//mcs_totlength+=(double)tmp_primtrk_mcs_radlengths_reco[i_mcs]*14.;
			//float tmp_sigma_hl=(180./3.14)*sigma_hl(ke_ff, (float)tmp_primtrk_mcs_radlengths_reco[i_mcs]*14.); //unit[deg]
			//sum_primtrk_mcs_sigma_hl+=tmp_sigma_hl;
			//count_primtrk_mcs_sigma_hl+=1.;
			//}
			//res.scatterAngles();
			//tmp_primtrk_mcs_angles_reco.clear();
			//tmp_primtrk_mcs_radlengths_reco.clear();
			//reco MCS angles ---------------------------------------------------------------------------------------------------//



			std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTag);
			if(calovector.size() != 3)
				std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;

			//HY::Get the Calorimetry(s) from thisTrack
			std::vector<double> tmp_primtrk_dqdx;	
			std::vector<double> tmp_primtrk_resrange;	
			std::vector<double> tmp_primtrk_dedx;	
			std::vector<double> tmp_primtrk_hitx;	
			std::vector<double> tmp_primtrk_hity;	
			std::vector<double> tmp_primtrk_hitz;
			std::vector<double> tmp_primtrk_pitch;
			std::vector<int> tmp_primtrk_wid;
			//std::vector<double> tmp_primtrk_truth_Eng;
			for (auto & calo : calovector) {
				if (calo.PlaneID().Plane == 2){ //only collection plane
					primtrk_range.push_back(calo.Range());
					std::cout<<"primtrk_range:"<<calo.Range()<<std::endl;
					for (size_t ihit = 0; ihit < calo.dQdx().size(); ++ihit){ //loop over hits
						tmp_primtrk_dqdx.push_back(calo.dQdx()[ihit]);
						tmp_primtrk_resrange.push_back(calo.ResidualRange()[ihit]);
						tmp_primtrk_dedx.push_back(calo.dEdx()[ihit]);
						tmp_primtrk_pitch.push_back(calo.TrkPitchVec()[ihit]);

						const auto &primtrk_pos=(calo.XYZ())[ihit];
						tmp_primtrk_hitx.push_back(primtrk_pos.X());
						tmp_primtrk_hity.push_back(primtrk_pos.Y());
						tmp_primtrk_hitz.push_back(primtrk_pos.Z());

                                                geo::TPCID tpc = fGeometry->FindTPCAtPosition(primtrk_pos);
						if(tpc.isValid){
							int tpc_no=tpc.TPC;
							geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
							geo::WireID wireID;
							try{
                                                                wireID = fWireReadoutGeom.Plane(planeID).NearestWireID(primtrk_pos);
							}
							catch(geo::InvalidWireError const& e) {
								wireID = e.suggestedWireID(); // pick the closest valid wire
							}
							tmp_primtrk_wid.push_back(wireID.Wire);
						}
						if(!tpc.isValid){
							tmp_primtrk_wid.push_back(-9999);
						}

						//geo::WireID WireID_reco=fWireReadoutGeom.NearestWireID(pos_reco, 2); //2 for collection plane
						//if (!WireID_reco) WireID_reco = fWireReadoutGeom.Plane(2).ClosestWireID(WireID_reco);
						//tmp_primtrk_wid.push_back(WireID_reco.Wire);
						//std::cout<<"(x,y,z)_reco:("<<primtrk_pos.X()<<","<<primtrk_pos.Y()<<","<<primtrk_pos.Z()<<") wid:"<<WireID_reco.Wire<<std::endl;
						//std::cout<<"dqdx="<<calo.dQdx()[ihit]<<"; resrange="<<calo.ResidualRange()[ihit]<<std::endl;
						//std::cout<<"(X,Y,Z)="<<"("<<calo.XYZ()[ihit].X()<<","<<calo.XYZ()[ihit].Y()<<","<<calo.XYZ()[ihit].Z()<<")"<<std::endl;

					} //loop over hits
				} //only collection plane
			}

			//HY::Associtate the reco info with the truth info using backtracker (only track, no shower) -------------------------------------------------------------------------------//
			// Get the truth info
			//const simb::MCParticle* mcparticle2 = truthUtil.GetMCParticleFromRecoTrack(*thisTrack, evt, fTrackerTag);
			const simb::MCParticle* geantGoodParticle1 = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
			//if(mcparticle2 != 0x0) { //mcparticle

			//std::vector<unsigned char> tmp_primtrk_hit_processkey;	
			std::vector<double> tmp_primtrk_hitx_true;	
			std::vector<double> tmp_primtrk_hity_true;	
			std::vector<double> tmp_primtrk_hitz_true;
			std::vector<double> tmp_primtrk_trkid_true;
			std::vector<double> tmp_primtrk_edept_true;

			std::vector<double> tmp_primtrk_true_x;	
			std::vector<double> tmp_primtrk_true_y;	
			std::vector<double> tmp_primtrk_true_z;
			std::vector<double> tmp_primtrk_true_trkid;
			std::vector<double> tmp_primtrk_true_edept;
			std::vector<int> tmp_primtrk_true_wid;	
			//std::vector<int> tmp_primtrk_true_tpc;	
			//std::vector<int> tmp_primtrk_true_plane;	

			std::vector<double> tmp_primtrj_true_x;
			std::vector<double> tmp_primtrj_true_y;
			std::vector<double> tmp_primtrj_true_z;
			std::vector<double> tmp_primtrj_true_edept;

			if(geantGoodParticle1!= 0x0 && geantGoodParticle1->Process()=="primary") { //sansity check
				//const simb::MCParticle *geantGoodParticle=pi_serv->TrackIdToMotherParticle_P(mcparticle2->TrackId());
				//if(geantGoodParticle != 0x0 && geantGoodParticle->Process()=="primary") { //geatGoodParticle and primary p Loop
				double ke_reco=ke_ff; //ke_reco at ff
				double ke_true=ke_ff;	//ke_true at ff

				art::ServiceHandle<cheat::BackTrackerService> bt_serv;
				simb::MCTrajectory truetraj=geantGoodParticle1->Trajectory();
				auto thisTrajectoryProcessMap1 =  truetraj.TrajectoryProcesses();
				std::string int_label="";
				//double endz_true=-999;
				if (thisTrajectoryProcessMap1.size()){
					for(auto const& couple: thisTrajectoryProcessMap1){
						//endz_true=((truetraj.at(couple.first)).first).Z();
						int_label=truetraj.KeyToProcess(couple.second);

						fprimary_truth_EndPosition[0]=((truetraj.at(couple.first)).first).X();
						fprimary_truth_EndPosition[1]=((truetraj.at(couple.first)).first).Y();
						fprimary_truth_EndPosition[2]=((truetraj.at(couple.first)).first).Z();
						fprimary_truth_EndProcess=truetraj.KeyToProcess(couple.second);	

						fprimary_truth_EndMomentum[0]=((truetraj.at(couple.first)).second).X();
						fprimary_truth_EndMomentum[1]= ((truetraj.at(couple.first)).second).Y();
						fprimary_truth_EndMomentum[2]=((truetraj.at(couple.first)).second).Z();

						break;
					}
				}
				primtrk_end_g4process=int_label;

				//study of interaction angles
				//if (thisTrajectoryProcessMap1.size()) { //TrajectoryProcessMap1
				if (thisTrajectoryProcessMap1.size()) { //TrajectoryProcessMap1
					for(auto const& couple1: thisTrajectoryProcessMap1) { //go through this traj with all the interaction vertices
                                                auto const& four_position = truetraj.at(couple1.first).first;
                                                interactionX.push_back(four_position.X());
                                                interactionY.push_back(four_position.Y());
                                                interactionZ.push_back(four_position.Z());
                                                interactionE.push_back(four_position.E());
						interactionProcesslist.push_back(truetraj.KeyToProcess(couple1.second));

						//get the TPC num 
                                                auto const pos = geo::vect::toPoint(four_position.Vect());
						unsigned int tpcno=1;
                                                auto const [xval, zval] = std::make_tuple(pos.X(), pos.Z());
						if(xval<=0 && zval<232) tpcno=1;
						if(xval<=0 && zval>232 && zval<464) tpcno=5; 
						if(xval<=0 && zval>=464) tpcno=9;
						if(xval>0 && zval<232) tpcno=2; 
						if(xval>0 && zval>232 && zval<464) tpcno=6; 
						if(xval>0 && zval>=464) tpcno=10;

                                                geo::TPCID const tpcID{0, tpcno};
						//convert the position of the interaction vertex to (wireID, peak time)
                                                interaction_wid_c.push_back(fWireReadoutGeom.Plane(geo::PlaneID{tpcID, 2}).WireCoordinate(pos));
                                                interaction_wid_v.push_back(fWireReadoutGeom.Plane(geo::PlaneID{tpcID, 1}).WireCoordinate(pos));
                                                interaction_wid_u.push_back(fWireReadoutGeom.Plane(geo::PlaneID{tpcID, 0}).WireCoordinate(pos));

                                                interaction_tt_c.push_back(detProp.ConvertXToTicks(pos.X(), geo::PlaneID{tpcID, 2}));
                                                interaction_tt_v.push_back(detProp.ConvertXToTicks(pos.X(), geo::PlaneID{tpcID, 1}));
                                                interaction_tt_u.push_back(detProp.ConvertXToTicks(pos.X(), geo::PlaneID{tpcID, 0}));

                                                //interactionT.push_back(detProp.ConvertXToTicks(xval, 2, tpcno, 0));
						//interactionU.push_back(fWireReadoutGeom.WireCoordinate(((truetraj.at(couple1.first)).first).Y(), ((truetraj.at(couple1.first)).first).Z(),0, tpcno, 0));
						//interactionV.push_back(fWireReadoutGeom.WireCoordinate(((truetraj.at(couple1.first)).first).Y(), ((truetraj.at(couple1.first)).first).Z(),1, tpcno, 0));
						//interactionW.push_back(fWireReadoutGeom.WireCoordinate(((truetraj.at(couple1.first)).first).Y(), ((truetraj.at(couple1.first)).first).Z(),2, tpcno, 0));

						//convert(x,y,z) to (wireid ,time ticks)
						//geo::PlaneID collection_plane = geom->PlaneIDs(2);
						//std::cout<<"fprimaryT0:"<<fprimaryT0<<std::endl;
						//std::cout<<"\nint_vtx x/y/z:"<<((truetraj.at(couple1.first)).first).X()<<"/"<<((truetraj.at(couple1.first)).first).Y()<<"/"<<((truetraj.at(couple1.first)).first).Z()<<std::endl;
                                                //std::cout<<"(wid,tt)_c:"<<"("<<fWireReadoutGeom.WireCoordinate(yval, zval, 2, tpcno, 0)<<","<<detProp.ConvertXToTicks(xval, 2, tpcno, 0)<<")"<<std::endl;
                                                //std::cout<<"(wid,tt)_v:"<<"("<<fWireReadoutGeom.WireCoordinate(yval, zval, 1, tpcno, 0)<<","<<detProp.ConvertXToTicks(xval, 1, tpcno, 0)<<")"<<std::endl;
                                                //std::cout<<"(wid,tt)_u:"<<"("<<fWireReadoutGeom.WireCoordinate(yval, zval, 0, tpcno, 0)<<","<<detProp.ConvertXToTicks(xval, 0, tpcno, 0)<<")"<<std::endl;


						//not interested of CoulombScat	
						if ((truetraj.KeyToProcess(couple1.second)).find("CoulombScat")!= std::string::npos) continue;

						//check if the interaction is in the TPC
						auto     interactionPos4D =  (truetraj.at(couple1.first)).first ;        

						if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
						else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
						else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;

						///get the interaction angle here
						double interactionAngle = 999999.; // This needs to be changed
						//--------------------- Int Angle ---------------------------
						// Try to retreive the interaction angle
						auto  prevInteractionPos4D = (truetraj.at(couple1.first-1)).first ;
						auto  prevInteractionPos3D = prevInteractionPos4D.Vect() ;
						auto  interactionPos3D     = interactionPos4D.Vect() ;
						auto  distanceBtwPoint     = interactionPos3D - prevInteractionPos3D;

						//see if the next point exists
						if (truetraj.size() > couple1.first + 1) {
							// The particle doesn't die. No need to check for anything else.
							auto nextInteractionPos4D =  (truetraj.at(couple1.first+1)).first ;
							auto nextInteractionPos3D =  nextInteractionPos4D.Vect() ;
							auto distanceBtwPointNext =  nextInteractionPos3D - interactionPos3D;
							interactionAngles.push_back(TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  ));
						}
						else { // The particle has come to an end. Let's check the daugthers.
							if (geantGoodParticle1->NumberDaughters() == 0 ){
								interactionAngles.push_back(interactionAngle);
								break;
							}
							double maxAngle = 0.;
							int numberOfTroubleDaugther = 0;
							//Loop on the daughters 
							const sim::ParticleList& plist=pi_serv->ParticleList();
							sim::ParticleList::const_iterator itPart1=plist.begin();
							for(size_t dPart1=0;(dPart1<plist.size()) && (plist.begin()!=plist.end());++dPart1) {
								const simb::MCParticle* dPart=(itPart1++)->second;
								if (dPart->Mother()  != 1 ) continue;
								auto daugthTraj = dPart->Trajectory();
								//	  if (debug) std::cout<< dPart->PdgCode()
								//   <<" , length: "<< (dPart->Trajectory()).TotalLength() <<"\n";
								//<<"First Point: "<< ((daugthTraj[0].first).Vect()).X() <<","<< ((daugthTraj[0].first).Vect()).Y()<<","<<((daugthTraj[0].first).Vect()).Z() <<"\n";
								if ((dPart->NumberTrajectoryPoints () < 2 )||!(TMath::Abs(dPart->PdgCode())  == 13||TMath::Abs(dPart->PdgCode())  == 11 || TMath::Abs(dPart->PdgCode())  == 211 || TMath::Abs(dPart->PdgCode())  == 321 || TMath::Abs(dPart->PdgCode())  == 2212)||(daugthTraj.TotalLength() < 0.5 )) {
									interactionAngle=999999;
									continue;
								}

								numberOfTroubleDaugther++;        

								auto daughtFirstPt  =  ((daugthTraj[0]).first).Vect() ;
								auto daughtSecondPt =  ((daugthTraj[1]).first).Vect() ;
								auto distanceBtwPointNext = daughtSecondPt - daughtFirstPt;
								interactionAngle = TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag()));
								if ( maxAngle < interactionAngle ) maxAngle = interactionAngle;
								interactionAngle = maxAngle;
								// If the track finishes without visible daugthers, we're likely to see it. Give it a big angle!
								if (!numberOfTroubleDaugther) interactionAngle = 9999999.;//It's a huge angle: it's greater than the cut, so I don't remove this! But it also doesn't go into the angle plot
								if (interactionAngle < 0.0001) std::cout<<"-------------------------------------------------> \n";
								interactionAngles.push_back(interactionAngle);
								break;
							}
						} // The particle has come to an end. Let's check the daugthers.
					}  //go through this traj with all the interaction vertices
				} //TrajectoryProcessMap1

				//save all the mctraj info 
				std::cout<<"\n MCTrajectory of this prim. trk has:"<<truetraj.size()<<" hits!!"<<std::endl;
				for (size_t tt=0; tt<truetraj.size(); ++tt) { //loop over all the hits in this mc traj
					tmp_primtrj_true_x.push_back(truetraj.X(tt));
					tmp_primtrj_true_y.push_back(truetraj.Y(tt));
					tmp_primtrj_true_z.push_back(truetraj.Z(tt));
					tmp_primtrj_true_edept.push_back(truetraj.E(tt));
					//std::cout<<"[mctraj] x/y/z/E:"<<truetraj.X(tt)<<"/"<<truetraj.Y(tt)<<"/"<<truetraj.Z(tt)<<"/"<<truetraj.E(tt)<<std::endl;  
					//std::cout<<"[int] x/y/z/E:"<<interactionX.at(tt)<<"/"<<interactionY.at(tt)<<"/"<<interactionZ.at(tt)<<"/"<<interactionE.at(tt)<<std::endl;
				} //loop over all the hits in this mc traj

				//use backtracker to find the associtated true info
                                geo::View_t view = fWireReadoutGeom.Plane({0, 0, 2}).View();
				auto simIDE_prim=bt_serv->TrackIdToSimIDEs_Ps(geantGoodParticle1->TrackId(),view);
				std::map<double, sim::IDE> orderedSimIDE; //id & e
				for (auto& ide : simIDE_prim) orderedSimIDE[ide->z]= *ide; //order in z-direction
				auto inTPCPoint  = truetraj.begin(); 
				auto Momentum0   = inTPCPoint->second;
				//auto old_iter = orderedSimIDE.begin();
				//std::cout<<"True KE at front face: "<<ke_true<<"; Reconstructed KE at ff: "<<ke_reco<<std::endl;
				primtrk_ke_reco.push_back(ke_reco);
				primtrk_ke_true.push_back(ke_true);

				double trklen=0.0;
				bool run_me_onetime=true;
				double xi=0.0; double yi=0.0; double zi=0.0;
				int cnt=0;
				//sanity check on the start/end positions
				if(tmp_primtrk_dqdx.size()!=0) { //make sure the dqdx vector is not empty
					if(tmp_primtrk_hitz[0]>tmp_primtrk_hitz[tmp_primtrk_hitz.size()-1]) { //flip trk vectors if Pandora messes up the trk direction
						std::reverse(tmp_primtrk_hitz.begin(), tmp_primtrk_hitz.end());
						std::reverse(tmp_primtrk_hity.begin(), tmp_primtrk_hity.end());
						std::reverse(tmp_primtrk_hitx.begin(), tmp_primtrk_hitx.end());
						std::reverse(tmp_primtrk_pitch.begin(),tmp_primtrk_pitch.end());
						std::reverse(tmp_primtrk_dedx.begin(), tmp_primtrk_dedx.end());
						std::reverse(tmp_primtrk_dqdx.begin(), tmp_primtrk_dqdx.end());
						std::reverse(tmp_primtrk_resrange.begin(), tmp_primtrk_resrange.end());
					} //flip trk vectors if Pandora messes up the trk direction
					//p.s. simIDE can only get Edept, needs to calculate true KE by your own

					for(size_t idx1=0; idx1<tmp_primtrk_dqdx.size()-1; idx1++) { //energy deposition: reco hit loop
						ke_reco-=tmp_primtrk_pitch[idx1]*tmp_primtrk_dedx[idx1];
						double currentDepEng = 0.;
						for ( auto iter= orderedSimIDE.begin(); iter!= orderedSimIDE.end(); iter++) { //simIde loop ---------------------------------------------------------------------------//
							auto currentIde = iter->second;

							//calculation only inside TPC
							if(currentIde.z<minZ || currentIde.z > maxZ ) continue;
							else if (currentIde.x < minX || currentIde.x > maxX ) continue;
							else if (currentIde.y < minY || currentIde.y > maxY ) continue;

							//if(cnt==0) { //start positions
							if(cnt==0&&currentIde.trackID>=0) { //start positions
								fprimary_truth_StartPosition[0] = currentIde.x;
								fprimary_truth_StartPosition[1] = currentIde.y;
								fprimary_truth_StartPosition[2] = currentIde.z;

								//tmp_primtrk_hit_processkey.push_back(currentIde.ProcessToKey);
								tmp_primtrk_hitx_true.push_back(currentIde.x);
								tmp_primtrk_hity_true.push_back(currentIde.y);
								tmp_primtrk_hitz_true.push_back(currentIde.z);
								tmp_primtrk_trkid_true.push_back(currentIde.trackID);
								tmp_primtrk_edept_true.push_back(currentIde.energy);
							} //start positions

							//run the simIde loop only one time to sum over all trk length in each segment
							//if (currentIde.trackID!=-1) { //discard true shower info
							//if(cnt>0 && run_me_onetime) { //run at once
							if (currentIde.trackID>=0) { //discard true shower info (trackID<0)
								if(cnt>0 && run_me_onetime) { //run at once
									trklen+=TMath::Sqrt(std::pow(currentIde.x-xi,2)+std::pow(currentIde.y-yi,2)+std::pow(currentIde.z-zi,2));
									//std::cout<<"trklen:  "<<trklen<<std::endl;
									//std::cout<<"dx, dy and dz: "<<currentIde.x-xi<<" "<<currentIde.y-yi<<"  "<<currentIde.z-zi<<std::endl;
									//std::cout<<"x,y,z: "<<currentIde.x<<", "<<currentIde.y<<", "<<currentIde.z<<"; trkId: "<<currentIde.trackID<<"; energy: "<<currentIde.energy<<std::endl;

									//fprimary_truth_EndPosition[0]=currentIde.x;
									//fprimary_truth_EndPosition[1]=currentIde.y;
									//fprimary_truth_EndPosition[2]=currentIde.z;

									//tmp_primtrk_hit_processkey.push_back(currentIde.ProcessToKey);
									tmp_primtrk_hitx_true.push_back(currentIde.x);
									tmp_primtrk_hity_true.push_back(currentIde.y);
									tmp_primtrk_hitz_true.push_back(currentIde.z);
									tmp_primtrk_trkid_true.push_back(currentIde.trackID);
									tmp_primtrk_edept_true.push_back(currentIde.energy);
								} //run at once
								xi=currentIde.x; yi=currentIde.y; zi=currentIde.z;
								cnt++;
								//std::cout<<"cnt:"<<cnt<<std::endl;
							} //discard true shower info

							//true E dept within the thin slice (simIde can only get Edept, not KE)
							if ( currentIde.z <= tmp_primtrk_hitz[idx1]) continue;
							if ( currentIde.z > tmp_primtrk_hitz[idx1+1]) continue;
							currentDepEng += currentIde.energy; //total energy deposition in the current slice
						} //simIde loop -------------------------------------------------------------------------------------------------------------------------------------------------------//
						ke_true -= currentDepEng; //KE in the current slice
						run_me_onetime=false; //finished the summation of the true trk length

						if(currentDepEng>0.001) { //remove the zombie tracks
							primtrk_ke_reco.push_back(ke_reco);      
							primtrk_ke_true.push_back(ke_true);
						}	//remove the zombie tracks

						//std::cout<<"primtrk_ke_reco["<<idx1<<"]:"<<ke_reco<<" ; primtrk_ke_true:"<<ke_true<<std::endl; 

						//if(currentDepEng>0.001) hKE_truth_reco->Fill(kineticEnergy,KE_rec);
						}//energy deposition: reco hit loop


						//sime IDE loop again to save all the true info 
						for ( auto iter2=orderedSimIDE.begin(); iter2!=orderedSimIDE.end(); iter2++) { //simIde loop2 ---------------------------------------------------------------------------//
							auto currentIde2 = iter2->second;

							//calculation only inside TPC
							if(currentIde2.z<minZ || currentIde2.z > maxZ ) continue;
							else if (currentIde2.x < minX || currentIde2.x > maxX ) continue;
							else if (currentIde2.y < minY || currentIde2.y > maxY ) continue;

							//if(currentIde2.trackID>=0) { // no shower
							//tmp_primtrk_hit_processkey.push_back(currentIde2.ProcessToKey);
							tmp_primtrk_true_x.push_back(currentIde2.x);
							tmp_primtrk_true_y.push_back(currentIde2.y);
							tmp_primtrk_true_z.push_back(currentIde2.z);
							tmp_primtrk_true_trkid.push_back(currentIde2.trackID);
							tmp_primtrk_true_edept.push_back(currentIde2.energy);

                                                        geo::Point_t const pos_true{currentIde2.x, currentIde2.y, currentIde2.z};
							geo::TPCID tpc = fGeometry->FindTPCAtPosition(pos_true);
							if(tpc.isValid){
								int tpc_no=tpc.TPC;
								geo::PlaneID planeID = geo::PlaneID(0, tpc_no, 2);
								geo::WireID wireID;
								try{
                                                                        wireID = fWireReadoutGeom.Plane(planeID).NearestWireID(pos_true);
								}
								catch(geo::InvalidWireError const& e) {
									wireID = e.suggestedWireID(); // pick the closest valid wire
								}
								tmp_primtrk_true_wid.push_back(wireID.Wire);
							}
							if(!tpc.isValid){
								tmp_primtrk_true_wid.push_back(-9999);
							}


							//geo::WireID WireID_true=fWireReadoutGeom.NearestWireID(pos_true, 2); //2 for collection plane
							//if (!WireID_true) WireID_true = fWireReadoutGeom.Plane(2).ClosestWireID(WireID_true);

							//geo::PlaneID Plane_true(fWireReadoutGeom.FindTPCAtPosition(pos_true),geo::kZ);
							//WireID_true=fWireReadoutGeom.NearestWireID(pos_true, Plane_true);
							//geo::WireID wireID = fWireReadoutGeom.NearestWireID(point, planeID);
							//if (!wireID) wireID = fWireReadoutGeom.Plane(planeID).ClosestWireID(wireID);

							//tmp_primtrk_true_wid.push_back(WireID_true.Wire);
							//tmp_primtrk_true_tpc.push_back(WireID_true.TPC);
							//tmp_primtrk_true_plane.push_back(WireID_true.Plane);

							//if (fRun==22655920&&fevent==791&&fprimaryID==41) {
							//std::cout<<"run:"<<fRun<<"/"<<"evt:"<<fevent<<"/primaryID:"<<fprimaryID<<
							//"/ WireID_true:"<<WireID_true<<"| wid:"<<WireID_true.Wire<<"| tpc:"<<WireID_true.TPC<<"| plane:"<<WireID_true.Plane<<"(x,y,z):("<<currentIde2.x<<","<<currentIde2.y<<","<<currentIde2.z<<")"<<std::endl;
							//std::cout<<"run:"<<fRun<<"/"<<"evt:"<<fevent<<"/primaryID:"<<fprimaryID<<
							//"/ wid:"<<WireID_true.Wire<<"| (x,y,z):("<<currentIde2.x<<","<<currentIde2.y<<","<<currentIde2.z<<")"<<std::endl;
							//}



							//std::cout<<"[simeide2] x,y,z: "<<currentIde2.x<<", "<<currentIde2.y<<", "<<currentIde2.z<<"; trkId: "<<currentIde2.trackID<<"; energy: "<<currentIde2.energy<<std::endl;
							//} //no shower
						} //simIde loop2 -------------------------------------------------------------------------------------------------------------------------------------------------------//


						} //make sure that primtrk_tmp size g.t. 0

						//} //geatGoodParticle and primary p Loop

						//std::cout<<"primtrk_range_true:"<<trklen<<std::endl;
						primtrk_range_true.push_back(trklen);

					} //mcpartice      
					//---------------------------------------------------------------------------------------------------------------//


					primtrk_dqdx.push_back(tmp_primtrk_dqdx);
					primtrk_resrange.push_back(tmp_primtrk_resrange);
					primtrk_dedx.push_back(tmp_primtrk_dedx);
					primtrk_hitx.push_back(tmp_primtrk_hitx);
					primtrk_hity.push_back(tmp_primtrk_hity);
					primtrk_hitz.push_back(tmp_primtrk_hitz);
					primtrk_pitch.push_back(tmp_primtrk_pitch);
					primtrk_wid.push_back(tmp_primtrk_wid);

					//primtrk_hit_processkey.push_back(tmp_primtrk_hit_processkey);
					primtrk_hitx_true.push_back(tmp_primtrk_hitx_true);
					primtrk_hity_true.push_back(tmp_primtrk_hity_true);
					primtrk_hitz_true.push_back(tmp_primtrk_hitz_true);
					primtrk_trkid_true.push_back(tmp_primtrk_trkid_true);
					primtrk_edept_true.push_back(tmp_primtrk_edept_true);

					primtrk_true_x.push_back(tmp_primtrk_true_x);
					primtrk_true_y.push_back(tmp_primtrk_true_y);
					primtrk_true_z.push_back(tmp_primtrk_true_z);
					primtrk_true_trkid.push_back(tmp_primtrk_true_trkid);
					primtrk_true_edept.push_back(tmp_primtrk_true_edept);

					primtrk_true_wid.push_back(tmp_primtrk_true_wid);
					//primtrk_true_tpc.push_back(tmp_primtrk_true_tpc);
					//primtrk_true_plane.push_back(tmp_primtrk_true_plane);

					primtrj_true_x.push_back(tmp_primtrj_true_x);
					primtrj_true_y.push_back(tmp_primtrj_true_y);
					primtrj_true_z.push_back(tmp_primtrj_true_z);
					primtrj_true_edept.push_back(tmp_primtrj_true_edept);

					tmp_primtrk_dqdx.clear();
					tmp_primtrk_resrange.clear();
					tmp_primtrk_dedx.clear();
					tmp_primtrk_hitx.clear();
					tmp_primtrk_hity.clear();
					tmp_primtrk_hitz.clear();
					tmp_primtrk_pitch.clear();
					tmp_primtrk_wid.clear();

					//tmp_primtrk_hit_processkey.clear();
					tmp_primtrk_hitx_true.clear();
					tmp_primtrk_hity_true.clear();
					tmp_primtrk_hitz_true.clear();
					tmp_primtrk_trkid_true.clear();
					tmp_primtrk_edept_true.clear();

					tmp_primtrk_true_x.clear();
					tmp_primtrk_true_y.clear();
					tmp_primtrk_true_z.clear();
					tmp_primtrk_true_trkid.clear();
					tmp_primtrk_true_edept.clear();
					tmp_primtrk_true_wid.clear();
					//tmp_primtrk_true_tpc.clear();
					//tmp_primtrk_true_plane.clear();

					tmp_primtrj_true_x.clear();
					tmp_primtrj_true_y.clear();
					tmp_primtrj_true_z.clear();
					tmp_primtrj_true_edept.clear();

					for(size_t k = 0; k < calovector.size() && k<3; k++){
						fprimaryKineticEnergy[k] = calovector[k].KineticEnergy();
						fprimaryRange[k] = calovector[k].Range();
						//const std::vector< double > & dedxvec = calovector[k].dEdx();
					}

					//Get CNN score of each hit ------------------------------------------------------------------------//
					float zpos=-999;
					float ypos=-999;
					float xpos=-999;

					//float max_inel_score_c=-999.;
					//float max_el_score_c=-999.;
					if(fmthm.isValid()){ //if non-empty fmthm
						auto vhit=fmthm.at(fprimaryID);
						auto vmeta=fmthm.data(fprimaryID); //indices of meta data are the same as data 
						for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
							bool fBadhit = false;
							if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
								fBadhit = true;
								//cout<<"fBadHit"<<fBadhit<<endl;
								continue;
							}
							if (vmeta[ii]->Index()>=tracklist[fprimaryID]->NumberTrajectoryPoints()){
								throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[fprimaryID]->NumberTrajectoryPoints()<<" for track index "<<fprimaryID<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
							}
							if (!tracklist[fprimaryID]->HasValidPoint(vmeta[ii]->Index())){
								fBadhit = true;
								// cout<<"had valid point "<<fBadhit<<endl;
								continue;
							}

							//get (x,y,z) 
							auto loc = tracklist[fprimaryID]->LocationAtPoint(vmeta[ii]->Index());
							xpos=loc.X();
							ypos=loc.Y();
							zpos=loc.Z();
							//std::cout<<"x, y, z: "<<xpos<<"  "<<ypos<<"  "<<zpos<<std::endl;
							//std::cout<<"BadHit"<<fBadhit<<std::endl;

							//get the TPC num 
							unsigned int tpc_no=1;
							if(xpos<=0 && zpos<232) tpc_no=1;
							if(xpos<=0 && zpos>232 && zpos<464) tpc_no=5; 
							if(xpos<=0 && zpos>=464) tpc_no=9;
							if(xpos>0 && zpos<232) tpc_no=2; 
							if(xpos>0 && zpos>232 && zpos<464) tpc_no=6; 
							if(xpos>0 && zpos>=464) tpc_no=10;

							//skip the bad hit if any
							if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
							if (zpos<-100) continue; //hit not on track
                                                        unsigned int const planenum = vhit[ii]->WireID().Plane;

                                                        geo::PlaneID const planeID{0, tpc_no, planenum};
							if(planenum==2){
								//std::cout<<"inside loop"<<std::endl;
								//std::array<float,3> cnn_out=hitResults.getOutput(vhit[ii]);
								//peakT_2.push_back(vhit[ii]->PeakTime());
								//int_2.push_back(vhit[ii]->Integral());
								//hitz_2.push_back(zpos);
								//inelscore_c.push_back(cnn_out[hitResults.getIndex("inel")]);
								//elscore_c.push_back(cnn_out[hitResults.getIndex("el")]);
								//nonescore_c.push_back(cnn_out[hitResults.getIndex("none")]);

								//save the associtated 3d position
								x_c.push_back(xpos);
								y_c.push_back(ypos);
								z_c.push_back(zpos);
								//std::cout<<"(x_c/y_c/z_c): ("<<xpos<<","<<ypos<<","<<zpos<<")"<<std::endl;

								//convert the position of the interaction to (wireID, peak time)
                                                                wid_c.push_back(fWireReadoutGeom.Plane(planeID).WireCoordinate(loc));
                                                                tt_c.push_back(detProp.ConvertXToTicks(xpos, planeID));

								//save ch number
								ch_c.push_back(vhit[ii]->Channel());

								pt_c.push_back(vhit[ii]->PeakTime());
								q_c.push_back(vhit[ii]->Integral());
								a_c.push_back(vhit[ii]->PeakAmplitude());
								wireno_c.push_back(vhit[ii]->WireID().Wire);

								//get the z-position (start/end) of associtated wire 
								double xyzStart[3];
								double xyzEnd[3];
								unsigned int wireno=vhit[ii]->WireID().Wire;
                                                                fWireReadoutGeom.WireEndPoints(geo::WireID(0,vhit[ii]->WireID().TPC,2,wireno), xyzStart, xyzEnd);
								wirez_c.push_back(xyzStart[2]);


                                                                //std::cout<<"(w,t):("<<fWireReadoutGeom.WireCoordinate(ypos, zpos, planenum, tpc_no, 0)<<","<<detProp.ConvertXToTicks(xpos, planenum, tpc_no, 0)<<")|"<<
								//"[inel,el,non]:["<<cnn_out[hitResults.getIndex("inel")]<<","<<cnn_out[hitResults.getIndex("el")]<<","<<cnn_out[hitResults.getIndex("none")]<<"]"<<std::endl;
								// std::cout<<"peaktime "<<vhit[ii]->PeakTime()<<std::endl;	

								//if (cnn_out[hitResults.getIndex("inel")]>max_inel_score_c){ //select max. inel_score
								//max_inel_score_c=cnn_out[hitResults.getIndex("inel")];
								//xyz_inelscore_c[0]=xpos;	
								//xyz_inelscore_c[1]=ypos;	
								//xyz_inelscore_c[2]=zpos;	
								//} //select max. inel_score

								//if (cnn_out[hitResults.getIndex("el")]>max_el_score_c){ //select max. el_score
								//max_el_score_c=cnn_out[hitResults.getIndex("el")];
								//xyz_elscore_c[0]=xpos;	
								//xyz_elscore_c[1]=ypos;	
								//xyz_elscore_c[2]=zpos;	
								//} //select max. el_score
							}//planenum 2
							if(planenum==1){
								//int_1.push_back(vhit[ii]->Integral());
								//hitz_1.push_back(zpos);

								//convert the position of the interaction to (wireID, peak time)
                                                                wid_v.push_back(fWireReadoutGeom.Plane(planeID).WireCoordinate(loc));
                                                                tt_v.push_back(detProp.ConvertXToTicks(xpos, planeID));

								//save ch number
								ch_v.push_back(vhit[ii]->Channel());

								pt_v.push_back(vhit[ii]->PeakTime());
								q_v.push_back(vhit[ii]->Integral());
								a_v.push_back(vhit[ii]->PeakAmplitude());
								wireno_v.push_back(vhit[ii]->WireID().Wire);

								//get the z-position (start/end) of associtated wire 
								double xyzStart[3];
								double xyzEnd[3];
								unsigned int wireno=vhit[ii]->WireID().Wire;
                                                                fWireReadoutGeom.WireEndPoints(geo::WireID(0,vhit[ii]->WireID().TPC,1,wireno), xyzStart, xyzEnd);
								wirez_v.push_back(xyzStart[2]);
							}//planenum 1
							if(planenum==0){
								//int_0.push_back(vhit[ii]->Integral());
								//hitz_0.push_back(zpos);

								//convert the position of the interaction to (wireID, peak time)
                                                                wid_u.push_back(fWireReadoutGeom.Plane(planeID).WireCoordinate(loc));
                                                                tt_u.push_back(detProp.ConvertXToTicks(xpos, planeID));

								//save ch number
								ch_u.push_back(vhit[ii]->Channel());

								pt_u.push_back(vhit[ii]->PeakTime());
								q_u.push_back(vhit[ii]->Integral());
								a_u.push_back(vhit[ii]->PeakAmplitude());
								wireno_u.push_back(vhit[ii]->WireID().Wire);

								//get the z-position (start/end) of associtated wire 
								double xyzStart[3];
								double xyzEnd[3];
								unsigned int wireno=vhit[ii]->WireID().Wire;
                                                                fWireReadoutGeom.WireEndPoints(geo::WireID(0,vhit[ii]->WireID().TPC,0,wireno), xyzStart, xyzEnd);
								wirez_u.push_back(xyzStart[2]);
							}//planenum 0



						} //loop over all meta data hit
					} //if non-empty fmthm




					//Get CNN score of each hit ------------------------------------------------------------------------//


					//Use Ajib's intersection calculation function
					//[1]identify the beam track and tag other tracks
					std::vector<float> Stw, Endw, Stt, Endt, Stwires, Endwires, Stticks, Endticks, TPCb, TPCcl;
					Stw.clear(); Endw.clear(); Stt.clear();  Endt.clear();  Stwires.clear();  Endwires.clear(); Stticks.clear(); Endticks.clear(); TPCb.clear(); TPCcl.clear();
					float den;
					float numw, numt,wire_no,ticks_no;
					for(size_t p1=0;p1<pfplist.size();p1++){
						std::vector<art::Ptr<recob::Track>> trk=pftrack.at(p1);
						std::vector<art::Ptr<recob::Cluster>> allClusters=fmcp.at(p1);
						//std::cout<<"fprimaryID:"<<fprimaryID<<std::endl;
						for(size_t c1=0;c1<allClusters.size();c1++){
							if(allClusters[c1]->Plane().Plane!=2) continue;
							if(trk.size() && int(trk[0].key())==fprimaryID){
								Stw.push_back(allClusters[c1]->StartWire());
								Endw.push_back(allClusters[c1]->EndWire());
								Stt.push_back(allClusters[c1]->StartTick());
								Endt.push_back(allClusters[c1]->EndTick());
								TPCb.push_back(allClusters[c1]->Plane().TPC);
								//std::cout<<"\nst/end wire:"<<allClusters[c1]->StartWire()<<"/"<<allClusters[c1]->EndWire()<<std::endl;
							}
							else{
								Stwires.push_back(allClusters[c1]->StartWire());
								Endwires.push_back(allClusters[c1]->EndWire());
								Stticks.push_back(allClusters[c1]->StartTick());
								Endticks.push_back(allClusters[c1]->EndTick());
								TPCcl.push_back(allClusters[c1]->Plane().TPC);
							}
						}
					}
					//[2]find interaction points if any (assuming all tracks are straight, find the interaction points)
					for(size_t clt=0;clt<Stw.size();clt++){
						for(size_t cl1=0;cl1<Stwires.size();cl1++){
							if(TPCcl[cl1]!=TPCb[clt]) continue;
							//std::cout<<"tpc are equal "<<std::endl;
							den=(Stw[clt]-Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]-Endwires[cl1]);
							if(den==0) continue;
							numw=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stwires[cl1]-Endwires[cl1])-(Stw[clt]-Endw[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
							numt=(Stw[clt]*Endt[clt]-Stt[clt]*Endw[clt])*(Stticks[cl1]-Endticks[cl1])-(Stt[clt]-Endt[clt])*(Stwires[cl1]*Endticks[cl1]-Stticks[cl1]*Endwires[cl1]);
							wire_no=numw/den;
							ticks_no=numt/den;
							//  std::cout<<"wireno and ticks not solution "<<wire_no<<"  "<<ticks_no<<std::endl;
							if(((Stw[clt]<wire_no && Endw[clt]>wire_no)||(Stw[clt]>wire_no && Endw[clt]<wire_no))&&((Stt[clt]<ticks_no && Endt[clt]>ticks_no)||(Stt[clt]>ticks_no && Endt[clt]<ticks_no)) && ((Stwires[cl1]<wire_no && Endwires[cl1]>wire_no)||(Stwires[cl1]>wire_no && Endwires[cl1]<wire_no)) && ((Stticks[cl1]<ticks_no && Endticks[cl1]>ticks_no)||(Stticks[cl1]>ticks_no && Endticks[cl1]<ticks_no))) { 
								std::cout<<"intersection wire and ticks are "<<std::round(wire_no)<<"  "<<ticks_no<<" Stw Endw StT EndT "<<Stwires[cl1]<<" "<<Endwires[cl1]<<" "<<Stticks[cl1]<<" "<<Endticks[cl1]<<std::endl;
								double xyzStart[3];
								double xyzEnd[3];
								unsigned int wireno=std::round(wire_no);
								geo::WireID wireid(0,TPCb[clt],2,wireno);
                                                                fWireReadoutGeom.WireEndPoints(wireid, xyzStart, xyzEnd);
								std::cout<<"Z position of intersection = "<<xyzStart[2]<<" "<<xyzEnd[2]<<"  "<<wireno<<std::endl;
								Zintersection.push_back(xyzStart[2]);
								Yintersection.push_back(xyzStart[1]);
								Xintersection.push_back(xyzStart[0]);
								timeintersection.push_back(ticks_no);
							}

						}
					}





				} //this track
				else if(thisShower != 0x0){
					fisprimarytrack               = 0;
					fisprimaryshower              = 1;

					fprimaryID                    = thisShower->ID();
					fprimaryLength                = thisShower->Length();
					fprimaryShowerBestPlane       = thisShower->best_plane();
					fprimaryOpeningAngle          = thisShower->OpenAngle();
					fprimaryStartPosition[0]      = thisShower->ShowerStart().X();
					fprimaryStartPosition[1]      = thisShower->ShowerStart().Y();
					fprimaryStartPosition[2]      = thisShower->ShowerStart().Z();
					fprimaryStartDirection[0]     = thisShower->Direction().X();
					fprimaryStartDirection[1]     = thisShower->Direction().Y();
					fprimaryStartDirection[2]     = thisShower->Direction().Z();
					if( (thisShower->Energy()).size() > 0 )
						fprimaryShowerEnergy = thisShower->Energy()[0]; // thisShower->best_plane()
					if( (thisShower->MIPEnergy()).size() > 0 )
						fprimaryShowerMIPEnergy = thisShower->MIPEnergy()[0];
					if( (thisShower->dEdx()).size() > 0 )
						fprimaryShowerdEdx = thisShower->dEdx()[0];
				}
				else{
					std::cout << "INFO::Primary pfParticle is not track or shower. Skip!" << std::endl;
					continue;
				}

				// Find the particle vertex. We need the tracker tag here because we need to do a bit of
				// additional work if the PFParticle is track-like to find the vertex. 
				const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);
				fvertex[0] = vtx.X(); fvertex[1] = vtx.Y(); fvertex[2] = vtx.Z();

				// Now we can look for the interaction point of the particle if one exists, i.e where the particle
				// scatters off an argon nucleus. Shower-like objects won't have an interaction point, so we can
				// check this by making sure we get a sensible position
				const TVector3 interactionVtx = pfpUtil.GetPFParticleSecondaryVertex(*particle,evt,fPFParticleTag,fTrackerTag);
				fsecvertex[0] = interactionVtx.X(); fsecvertex[1] = interactionVtx.Y(); fsecvertex[2] = interactionVtx.Z();

				// Maximum number of daugthers to be processed
				if(particle->NumDaughters() > NMAXDAUGTHERS)
					std::cout << "INFO::Number of daughters is " << particle->NumDaughters() << ". Only the first NMAXDAUGTHERS are processed." << std::endl;

				// Let's get the daughter PFParticles... we can do this simply without the utility
				for(const int daughterID : particle->Daughters()){
					// Daughter ID is the element of the original recoParticle vector
					const recob::PFParticle *daughterParticle      = &(recoParticles->at(daughterID));
					std::cout << "Daughter " << daughterID << " has " << daughterParticle->NumDaughters() << " daughters" << std::endl;

					const recob::Track* daughterTrack              = pfpUtil.GetPFParticleTrack(*daughterParticle,evt,fPFParticleTag,fTrackerTag);
					const recob::Shower* daughterShower            = pfpUtil.GetPFParticleShower(*daughterParticle,evt,fPFParticleTag,fShowerTag);

					if(daughterTrack != 0x0){
						fisdaughtertrack[fNDAUGHTERS]                = 1;
						fisdaughtershower[fNDAUGHTERS]               = 0;
						fdaughterTheta[fNDAUGHTERS]                  = daughterTrack->Theta();
						fdaughterPhi[fNDAUGHTERS]                    = daughterTrack->Phi();
						fdaughterLength[fNDAUGHTERS]                 = daughterTrack->Length();
						fdaughterMomentum[fNDAUGHTERS]               = daughterTrack->StartMomentum();
						fdaughterEndMomentum[fNDAUGHTERS]            = daughterTrack->EndMomentum();
						fdaughterStartPosition[fNDAUGHTERS][0]       = daughterTrack->Trajectory().Start().X();
						fdaughterStartPosition[fNDAUGHTERS][1]       = daughterTrack->Trajectory().Start().Y();
						fdaughterStartPosition[fNDAUGHTERS][2]       = daughterTrack->Trajectory().Start().Z();
						fdaughterEndPosition[fNDAUGHTERS][0]         = daughterTrack->Trajectory().End().X();
						fdaughterEndPosition[fNDAUGHTERS][1]         = daughterTrack->Trajectory().End().Y();
						fdaughterEndPosition[fNDAUGHTERS][2]         = daughterTrack->Trajectory().End().Z();
						fdaughterStartDirection[fNDAUGHTERS][0]      = daughterTrack->Trajectory().StartDirection().X();
						fdaughterStartDirection[fNDAUGHTERS][1]      = daughterTrack->Trajectory().StartDirection().Y();
						fdaughterStartDirection[fNDAUGHTERS][2]      = daughterTrack->Trajectory().StartDirection().Z();
						fdaughterEndDirection[fNDAUGHTERS][0]        = daughterTrack->Trajectory().EndDirection().X();
						fdaughterEndDirection[fNDAUGHTERS][1]        = daughterTrack->Trajectory().EndDirection().Y();
						fdaughterEndDirection[fNDAUGHTERS][2]        = daughterTrack->Trajectory().EndDirection().Z();

						fdaughterMomentumByRangeMuon[fNDAUGHTERS]    = trmom.GetTrackMomentum(daughterTrack->Length(),13);
						fdaughterMomentumByRangeProton[fNDAUGHTERS]  = trmom.GetTrackMomentum(daughterTrack->Length(),2212);

						std::vector<anab::Calorimetry> daughtercalovector = trackUtil.GetRecoTrackCalorimetry(*daughterTrack, evt, fTrackerTag, fCalorimetryTag);
						if(daughtercalovector.size() != 3)
							std::cerr << "WARNING::Calorimetry vector size for daughter is = " << daughtercalovector.size() << std::endl;

						for(size_t k = 0; k < daughtercalovector.size() && k<3; k++){
							fdaughterKineticEnergy[fNDAUGHTERS][k] = daughtercalovector[k].KineticEnergy();
							fdaughterRange[fNDAUGHTERS][k] = daughtercalovector[k].Range();
						}

						// Get the true mc particle
                                                const simb::MCParticle* mcdaughterparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *daughterTrack, evt, fTrackerTag);
						if(mcdaughterparticle != 0x0){
							fdaughter_truth_TrackId[fNDAUGHTERS]          = mcdaughterparticle->TrackId();
							fdaughter_truth_Pdg[fNDAUGHTERS]              = mcdaughterparticle->PdgCode();
							fdaughter_truth_StartPosition[fNDAUGHTERS][0] = mcdaughterparticle->Vx();
							fdaughter_truth_StartPosition[fNDAUGHTERS][1] = mcdaughterparticle->Vy();
							fdaughter_truth_StartPosition[fNDAUGHTERS][2] = mcdaughterparticle->Vz();
							fdaughter_truth_StartPosition[fNDAUGHTERS][3] = mcdaughterparticle->T();
							fdaughter_truth_EndPosition[fNDAUGHTERS][0]   = mcdaughterparticle->EndX();
							fdaughter_truth_EndPosition[fNDAUGHTERS][1]   = mcdaughterparticle->EndY();
							fdaughter_truth_EndPosition[fNDAUGHTERS][2]   = mcdaughterparticle->EndZ();
							fdaughter_truth_EndPosition[fNDAUGHTERS][3]   = mcdaughterparticle->EndT();
							fdaughter_truth_P[fNDAUGHTERS]                = mcdaughterparticle->P();
							fdaughter_truth_Momentum[fNDAUGHTERS][0]      = mcdaughterparticle->Px();
							fdaughter_truth_Momentum[fNDAUGHTERS][1]      = mcdaughterparticle->Py();
							fdaughter_truth_Momentum[fNDAUGHTERS][2]      = mcdaughterparticle->Pz();
							fdaughter_truth_Momentum[fNDAUGHTERS][3]      = mcdaughterparticle->E();
							fdaughter_truth_Pt[fNDAUGHTERS]               = mcdaughterparticle->Pt();
							fdaughter_truth_Mass[fNDAUGHTERS]             = mcdaughterparticle->Mass();
							fdaughter_truth_EndMomentum[fNDAUGHTERS][0]   = mcdaughterparticle->EndPx();
							fdaughter_truth_EndMomentum[fNDAUGHTERS][1]   = mcdaughterparticle->EndPy();
							fdaughter_truth_EndMomentum[fNDAUGHTERS][2]   = mcdaughterparticle->EndPz();
							fdaughter_truth_EndMomentum[fNDAUGHTERS][3]   = mcdaughterparticle->EndE();
							fdaughter_truth_Theta[fNDAUGHTERS]            = mcdaughterparticle->Momentum().Theta();
							fdaughter_truth_Phi[fNDAUGHTERS]              = mcdaughterparticle->Momentum().Phi();
							fdaughter_truth_Process[fNDAUGHTERS]          = int(mcdaughterparticle->Trajectory().ProcessToKey(mcdaughterparticle->Process()));
							std::cout << "Daughter Process = " << (mcdaughterparticle->Process()).c_str() 
								<< " , mother = " << mcdaughterparticle->Mother() 
								<< std::endl;
						}
					}
					else if(daughterShower != 0x0){
						fisdaughtertrack[fNDAUGHTERS]                = 0;
						fisdaughtershower[fNDAUGHTERS]               = 1;
						fdaughterLength[fNDAUGHTERS]                 = daughterShower->Length();
						fdaughterShowerBestPlane[fNDAUGHTERS]        = daughterShower->best_plane();
						fdaughterOpeningAngle[fNDAUGHTERS]           = daughterShower->OpenAngle();
						fdaughterStartPosition[fNDAUGHTERS][0]       = daughterShower->ShowerStart().X();
						fdaughterStartPosition[fNDAUGHTERS][1]       = daughterShower->ShowerStart().Y();
						fdaughterStartPosition[fNDAUGHTERS][2]       = daughterShower->ShowerStart().Z();
						fdaughterStartDirection[fNDAUGHTERS][0]      = daughterShower->Direction().X();
						fdaughterStartDirection[fNDAUGHTERS][1]      = daughterShower->Direction().Y();
						fdaughterStartDirection[fNDAUGHTERS][2]      = daughterShower->Direction().Z();
						if( (daughterShower->Energy()).size() > 0 )
							fdaughterShowerEnergy[fNDAUGHTERS] = daughterShower->Energy()[0]; // thisShower->best_plane()
						if( (daughterShower->MIPEnergy()).size() > 0 )
							fdaughterShowerMIPEnergy[fNDAUGHTERS] = daughterShower->MIPEnergy()[0];
						if( (daughterShower->dEdx()).size() > 0 )
							fdaughterShowerdEdx[fNDAUGHTERS] = daughterShower->dEdx()[0];
					}
					else{
						std::cout << "INFO::Daughter pfParticle is not track or shower. Skip!" << std::endl;
						continue;
					}

					fdaughterID[fNDAUGHTERS]                       = daughterID;
					// NHits associated with this pfParticle
					fdaughterNHits[fNDAUGHTERS]                    = (pfpUtil.GetPFParticleHits(*daughterParticle,evt,fPFParticleTag)).size();
					// T0
					//std::vector<anab::T0> pfdaughterT0vec = pfpUtil.GetPFParticleT0(*daughterParticle,evt,fPFParticleTag);
					//if(!pfT0vec.empty())
					//	fdaughterT0[fNDAUGHTERS] = pfdaughterT0vec[0].Time();

					fNDAUGHTERS++;

					// Only process NMAXDAUGTHERS
					if(fNDAUGHTERS > NMAXDAUGTHERS) break;

				}

				// For actually studying the objects, it is easier to have the daughters in their track and shower forms.
				// We can use the utility to get a vector of track-like and a vector of shower-like daughters
				//const std::vector<const recob::Track*> trackDaughters = pfpUtil.GetPFParticleDaughterTracks(*particle,evt,fPFParticleTag,fTrackerTag);  
				//const std::vector<const recob::Shower*> showerDaughters = pfpUtil.GetPFParticleDaughterShowers(*particle,evt,fPFParticleTag,fShowerTag);

				// For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
				break;
			} 

			// Fill trees
			if(beamTriggerEvent)
				fPandoraBeam->Fill();

			//fPandoraCosmics->Fill();

			}

			void protoana::protonmc::endJob(){

			}

			void protoana::protonmc::FillCosmicsTree(art::Event const & evt, std::string pfParticleTag){

				// To fill

			}

			void protoana::protonmc::Initialise(){

				fRun = -999;
				fSubRun = -999;
				fevent = -999;
				fTimeStamp = -999.0;
				for(int k=0; k < 5; k++)
					fNactivefembs[k] = -999;

				for(int k=0; k < 3; k++){
					fvertex[k] = -999.0;
					fsecvertex[k] = -999.0;
					fprimaryEndPosition[k] = -999.0;
					fprimaryStartPosition[k] = -999.0;
					fprimaryEndDirection[k] = -999.0;
					fprimaryStartDirection[k] = -999.0;
					fprimaryKineticEnergy[k] = -999.0;
					fprimaryRange[k] = -999.0;

					//xyz_inelscore_c[k] =-999.;
					//xyz_elscore_c[k] =-999.;

				}

				fbeamtrigger = -999;
				ftof = -999.0;
				fcerenkov1 = -999;
				fcerenkov2 = -999;
				fbeamtrackMomentum = -999.0;
				fbeamtrackEnergy = 999.0;
				fbeamtrackPdg = -999;
				fbeamtrackTime = -999.0;
				fbeamtrackID = -999;
				for(int l=0; l < 3; l++){
					fbeamtrackP[l] = -999.0;
					fbeamtrackPos[l] = -999.0;
					fbeamtrackDir[l] = -999.0;
				}

				//NumberBeamTrajectoryPoints=0; 
				beamtrk_x.clear();
				beamtrk_y.clear();
				beamtrk_z.clear();
				beamtrk_Px.clear();
				beamtrk_Py.clear();
				beamtrk_Pz.clear();
				beamtrk_Eng.clear();


  				beamMomentum_spec.clear();
  				beamPosx_spec.clear();
  				beamPosy_spec.clear();
  				beamPosz_spec.clear();
  				beamDirx_spec.clear();
  				beamDiry_spec.clear();
  				beamDirz_spec.clear();


				x_c.clear();
				y_c.clear();
				z_c.clear();


				ke_ff=0.;
				pos_ff.clear();

				primtrk_ke_true.clear();
				primtrk_ke_reco.clear();
				primtrk_end_g4process="";
				primtrk_range_true.clear();

				//primtrk_mcs_angles_reco.clear();
				//primtrk_mcs_radlengths_reco.clear();

				fisprimarytrack = 0;
				fisprimaryshower = 0;
				fNDAUGHTERS = 0;

				fprimaryBDTScore = -999.0;
				fprimaryNHits = -999;
				fprimaryTheta = -999.0;
				fprimaryPhi = -999.0;
				fprimaryLength = -999.0;
				fprimaryMomentum = -999.0;
				fprimaryEndMomentum = -999.0;
				fprimaryOpeningAngle = -999.0;
				fprimaryShowerBestPlane = -999;
				fprimaryShowerEnergy = -999;
				fprimaryShowerMIPEnergy = -999;
				fprimaryShowerdEdx = -999;
				fprimaryID = -999;
				fprimaryMomentumByRangeProton = -999.0;
				fprimaryMomentumByRangeMuon = -999.0;
				fprimaryT0 = -999.0;

				fprimary_truth_TrackId = -999;
				fprimary_truth_Pdg = -999;
				fprimary_truth_P = -999.0;
				fprimary_truth_Pt = -999.0;
				fprimary_truth_Mass = -999.0;
				fprimary_truth_Theta = -999.0;
				fprimary_truth_Phi = -999.0;
				//fprimary_truth_Process = -999;
				fprimary_truth_Isbeammatched = -999;
				fprimary_truth_NDaughters = -999;

				//fprimary_truth_Process="";
				fprimary_truth_EndProcess="";
				//fprimary_backtrker_truth_Process="";
				//fprimary_backtrker_truth_EndProcess="";

				for(int l=0; l < 4; l++){
					fprimary_truth_StartPosition[l] = -999.0;
					fprimary_truth_StartPosition_MC[l] = -999.0;
					fprimary_truth_EndPosition[l] = -999.0;
					fprimary_truth_EndPosition_MC[l] = -999.0;
					fprimary_truth_Momentum[l] = -999.0;
					fprimary_truth_EndMomentum[l] = -999.0;
				}

				for(int k=0; k < NMAXDAUGTHERS; k++){
					fisdaughtertrack[k] = -999;
					fisdaughtershower[k] = -999;
					fdaughterNHits[k] = -999;
					fdaughterTheta[k] = -999.0;
					fdaughterPhi[k] = -999.0;
					fdaughterLength[k] = -999.0;
					fdaughterMomentum[k] = -999.0;
					fdaughterEndMomentum[k] = -999.0;
					for(int l=0; l < 3; l++){
						fdaughterEndPosition[k][l] = -999.0;
						fdaughterStartPosition[k][l] = -999.0;
						fdaughterEndDirection[k][l] = -999.0;
						fdaughterStartDirection[k][l] = -999.0;
						fdaughterKineticEnergy[k][l] = -999.0;
						fdaughterRange[k][l] = -999.0;
					}
					fdaughterOpeningAngle[k] = -999.0;
					fdaughterShowerBestPlane[k] = -999;
					fdaughterShowerEnergy[k] = -999;
					fdaughterShowerMIPEnergy[k] = -999;
					fdaughterShowerdEdx[k] = -999;
					fdaughterMomentumByRangeProton[k] = -999.0;
					fdaughterMomentumByRangeMuon[k] = -999.0;
					fdaughterID[k] = -999;
					//fdaughterT0[k] = -999;

					fdaughter_truth_TrackId[k] = -999;
					fdaughter_truth_Pdg[k] = -999;
					fdaughter_truth_P[k] = -999.0;
					fdaughter_truth_Pt[k] = -999.0;
					fdaughter_truth_Mass[k] = -999.0;
					fdaughter_truth_Theta[k] = -999.0;
					fdaughter_truth_Phi[k] = -999.0;
					fdaughter_truth_Process[k] = -999;
					for(int l=0; l < 4; l++){
						fdaughter_truth_StartPosition[k][l] = -999.0;
						fdaughter_truth_EndPosition[k][l] = -999.0;
						fdaughter_truth_Momentum[k][l] = -999.0;
						fdaughter_truth_EndMomentum[k][l] = -999.0;
					}
				}

				primtrk_dqdx.clear();
				primtrk_resrange.clear();
				primtrk_dedx.clear();
				primtrk_range.clear();
				primtrk_hitx.clear();
				primtrk_hity.clear();
				primtrk_hitz.clear();
				primtrk_pitch.clear();
				primtrk_wid.clear();

				//inelscore_c.clear();
				//elscore_c.clear();
				//nonescore_c.clear();


				ke_ff_true=0.;
				//primtrk_hit_processkey.clear();
				primtrk_hitx_true.clear();
				primtrk_hity_true.clear();
				primtrk_hitz_true.clear();
				primtrk_trkid_true.clear();
				primtrk_edept_true.clear();

				primtrk_true_x.clear();
				primtrk_true_y.clear();
				primtrk_true_z.clear();
				primtrk_true_trkid.clear();
				primtrk_true_edept.clear();
				primtrk_true_wid.clear();
				//primtrk_true_tpc.clear();
				//primtrk_true_plane.clear();

				primtrj_true_x.clear();
				primtrj_true_y.clear();
				primtrj_true_z.clear();
				primtrj_true_edept.clear();

				interactionX.clear();
				interactionY.clear();
				interactionZ.clear();
				interactionE.clear();
				interactionProcesslist.clear();
				interactionAngles.clear();

				Zintersection.clear();
				Yintersection.clear();
				Xintersection.clear();
				timeintersection.clear();

				pfphit_peaktime_c.clear();
				pfphit_peaktime_u.clear();
				pfphit_peaktime_v.clear();

				pfphit_wireid_c.clear();
				pfphit_wireid_u.clear();
				pfphit_wireid_v.clear();

				interaction_wid_c.clear();
				interaction_tt_c.clear();
				interaction_wid_v.clear();
				interaction_tt_v.clear();
				interaction_wid_u.clear();
				interaction_tt_u.clear();

				wid_c.clear();
				tt_c.clear();
				ch_c.clear();
				wirez_c.clear();
				pt_c.clear();
				q_c.clear();
				a_c.clear();
				wireno_c.clear();

				wid_v.clear();
				tt_v.clear();
				ch_v.clear();
				wirez_v.clear();
				pt_v.clear();
				q_v.clear();
				a_v.clear();
				wireno_v.clear();

				wid_u.clear();
				tt_u.clear();
				ch_u.clear();
				wirez_u.clear();
				pt_u.clear();
				q_u.clear();
				a_u.clear();
				wireno_u.clear();


				g4rw_primary_singular_weight=-1;
				g4rw_primary_weights.clear();
				g4rw_primary_plus_sigma_weight.clear();
				g4rw_primary_minus_sigma_weight.clear();
				g4rw_primary_var.clear();
				g4rw_set_weights.clear();
				g4rw_p1.clear();
				g4rw_p2.clear();


			}

			DEFINE_ART_MODULE(protoana::protonmc)
