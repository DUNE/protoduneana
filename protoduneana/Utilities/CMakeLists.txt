
art_make( LIBRARY_NAME ProtoDUNEUtilities
	  LIB_LIBRARIES
                        dunecore::ArtSupport
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larreco::RecoAlg
                        larreco::Calorimetry
                        larreco::RecoAlg_PMAlg
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        lardataalg::DetectorInfo
                        lardata::ArtDataHelper
                        lardata::Utilities
                        larsim::Simulation 
                        lardataobj::Simulation
                        duneprototypes::ProtoDUNEDataUtils
                        geant4reweight::ReweightBaseLib
                        geant4reweight::PropBaseLib
                        dunecore::DuneObj
                        duneopdet::OpticalDetector
                        dunepdlegacy::Services_ChannelMap_PdspChannelMapService_service
                        larsim::MCCheater_PhotonBackTrackerService_service
                        larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
                        nug4::ParticleNavigation 
                        nusimdata::SimulationBase
                        CLHEP::CLHEP
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common 
                        art::Persistency_Provenance
                        art::Utilities 
			canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
			cetlib_except::cetlib_except
                        ROOT::Core ROOT::Tree ROOT::Hist
	  MODULE_LIBRARIES
          ProtoDUNEUtilities
)

install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(FitUtils)
