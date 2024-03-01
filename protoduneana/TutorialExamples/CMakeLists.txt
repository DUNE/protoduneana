cet_build_plugin(ClusterCounter art::module LIBRARIES
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
			messagefacility::MF_MessageLogger
			ROOT::Core ROOT::Tree ROOT::Hist
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter2 art::module LIBRARIES
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
			messagefacility::MF_MessageLogger
			ROOT::Core ROOT::Tree ROOT::Hist
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter3 art::module LIBRARIES
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
			messagefacility::MF_MessageLogger
			ROOT::Core ROOT::Tree ROOT::Hist
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter4 art::module LIBRARIES
                        lardataalg::DetectorInfo
                        lardata::headers
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
			nusimdata::SimulationBase
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
			messagefacility::MF_MessageLogger
			ROOT::Core ROOT::Tree ROOT::Hist
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter5 art::module LIBRARIES
                        lardataalg::DetectorInfo
                        lardata::headers
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
			nusimdata::SimulationBase
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
			messagefacility::MF_MessageLogger
			ROOT::Core ROOT::Tree ROOT::Hist
			BASENAME_ONLY
)

cet_build_plugin(ReadSpacePointAndCnn art::module LIBRARIES
			lardata::ArtDataHelper
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
			nusimdata::SimulationBase
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
			messagefacility::MF_MessageLogger
			ROOT::Core ROOT::Tree ROOT::Hist
			BASENAME_ONLY
)

cet_build_plugin(ReadBeamInfo art::module LIBRARIES
      nusimdata::SimulationBase
      larcorealg::Geometry
      larcore::Geometry_Geometry_service
      lardata::Utilities
      larsim::MCCheater_BackTrackerService_service
      larsim::MCCheater_ParticleInventoryService_service
      art::Framework_Services_Registry
      art_root_io::tfile_support
      ROOT::Core
      art_root_io::TFileService_service
      messagefacility::MF_MessageLogger
      ROOT::Core ROOT::Tree ROOT::Hist
      BASENAME_ONLY
)

add_subdirectory(MoreExamples)

install_headers()
install_fhicl()
install_source()
install_scripts()
