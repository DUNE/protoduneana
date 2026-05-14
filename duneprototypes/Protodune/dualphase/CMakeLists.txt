

cet_build_plugin(
           CRTGen art::module LIBRARIES
           nusimdata::SimulationBase
           cetlib::cetlib
           cetlib_except::cetlib_except
           fhiclcpp::fhiclcpp
           CLHEP::CLHEP
           larcorealg::Geometry
           larcoreobj::SummaryData
           larcore::headers
           larsim::EventGenerator
           nurandom::RandomUtils_NuRandomService_service
           art::Framework_Services_Registry
           larsim::MCCheater_ParticleInventoryService_service
           art_root_io::tfile_support
           CLHEP::CLHEP
           messagefacility::MF_MessageLogger
  	   ROOT::Core ROOT::Hist ROOT::Tree
           ROOT::Geom
           ROOT::XMLIO
           ROOT::Gdml
)
if(FALSE)
cet_build_plugin(
			Purity art::module LIBRARIES
                        lardataobj::RecoBase
                        larreco::Calorimetry
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        lardata::ArtDataHelper
                        lardata::Utilities
                        lardataobj::Simulation
                        larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except::cetlib_except
                        ROOT::Core ROOT::Hist ROOT::Tree
                        BASENAME_ONLY
)
endif()
cet_build_plugin(
			DEdxdp art::module LIBRARIES
                        lardataobj::RecoBase
                        larreco::Calorimetry
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        lardata::ArtDataHelper
                        lardata::Utilities
                        lardataobj::Simulation
                        larsim::MCCheater_BackTrackerService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except::cetlib_except
                        ROOT::Core ROOT::Hist ROOT::Tree
		                PostgreSQL::PostgreSQL
                        BASENAME_ONLY
)

cet_build_plugin(
	  AnaRootParser art::module LIBRARIES
	  larcore::Geometry_Geometry_service
	  larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
	  larsim::MCCheater_BackTrackerService_service
          larsim::MCCheater_ParticleInventoryService_service
	  lardata::Utilities
	  larevt::Filters
	  lardataobj::RawData
	  lardataobj::RecoBase
	  lardata::RecoObjects
	  lardata::ArtDataHelper
	  lardataalg::DetectorInfo
	  larcorealg::Geometry
	  larcorealg::GeoAlgo
          larcore::headers
	  larreco::RecoAlg
	  larpandora::LArPandoraInterface
	  dunecore::SignalShapingServiceDUNE10kt_service
	  dunecore::SignalShapingServiceDUNE35t_service
          nusimdata::SimulationBase
          art::Framework_Core
	  art::Framework_Principal
	  art::Framework_Services_Registry
	  art_root_io::tfile_support
	  ROOT::Core

	  art::Persistency_Common
	  art::Persistency_Provenance
	  art::Utilities canvas::canvas
	  messagefacility::MF_MessageLogger
	  
	  fhiclcpp::fhiclcpp
	  cetlib::cetlib cetlib_except::cetlib_except
          ROOT::Geom
          ROOT::XMLIO
          ROOT::Gdml
	  ROOT::Core ROOT::Hist ROOT::Tree
)



install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(Light)
add_subdirectory(RawDecoding)
add_subdirectory(fcl)
#add_subdirectory(test)
