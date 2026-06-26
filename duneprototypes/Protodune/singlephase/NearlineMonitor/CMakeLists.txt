cet_build_plugin(SPLifetime art::module
              LIBRARIES
              lardataobj::RecoBase
              lardataalg::DetectorInfo
              lardata::headers
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art_root_io::tfile_support
	      ROOT::Core ROOT::Graf ROOT::Gpad
              art::Persistency_Provenance
              art::Persistency_Common
              art::Utilities
              canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except::cetlib_except
              ROOT::MathCore
              ROOT::GenVector
              ROOT::Geom
              BASENAME_ONLY)

cet_build_plugin(PlotTrackT0 art::module
              LIBRARIES
              lardataobj::RecoBase
              lardataalg::DetectorInfo
              lardata::headers
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              art::Persistency_Provenance
              art::Persistency_Common
              art::Utilities
              canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except::cetlib_except
              ROOT::MathCore
              ROOT::GenVector
              ROOT::Geom
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

cet_build_plugin(PlotTrackDetails art::module
              LIBRARIES
              lardataobj::RecoBase
              lardataalg::DetectorInfo
              lardata::headers
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              art::Persistency_Provenance
              art::Persistency_Common
              art::Utilities
              canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except::cetlib_except
              ROOT::MathCore
              ROOT::GenVector
              ROOT::Geom
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

cet_build_plugin(PlotEventDetails art::module
              LIBRARIES
      lardataobj::RecoBase
              lardata::headers
              lardataalg::DetectorInfo
      larcorealg::Geometry
      larcore::Geometry_Geometry_service
      larsim::MCCheater_BackTrackerService_service
      larsim::MCCheater_ParticleInventoryService_service
      nusimdata::SimulationBase
      art::Framework_Core
      art::Framework_Principal
      art_root_io::tfile_support
      ROOT::Core
      art::Persistency_Provenance
      art::Persistency_Common
      art::Utilities 
      canvas::canvas
      messagefacility::MF_MessageLogger
      cetlib::cetlib cetlib_except::cetlib_except
      ROOT::MathCore
      ROOT::GenVector
      ROOT::Geom
      ROOT::Core ROOT::Hist ROOT::Tree
      BASENAME_ONLY
)

cet_build_plugin(PlotOpticalDetails art::module
              LIBRARIES
              lardataobj::RecoBase
              lardataalg::DetectorInfo
              lardata::headers
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              art::Persistency_Provenance
              art::Persistency_Common
              art::Utilities
              canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except::cetlib_except
              ROOT::MathCore
              ROOT::GenVector
              ROOT::Geom
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

cet_build_plugin(RawEventDisplay art::module
              LIBRARIES
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              lardataalg::DetectorInfo
              lardataobj::RawData
              lardata::headers
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

cet_build_plugin(TpcMonitor art::module
              LIBRARIES
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              lardataalg::DetectorInfo
              lardataobj::RawData
              lardata::headers
              nusimdata::SimulationBase
              dunepdlegacy::Services_ChannelMap_PdspChannelMapService_service
              dunecore::DuneObj
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

cet_build_plugin(SSPMonitor art::module
              LIBRARIES
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              lardataobj::RawData
              lardataobj::RecoBase
              lardata::headers
              nusimdata::SimulationBase
              dunepdlegacy::Services_ChannelMap_PdspChannelMapService_service
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

cet_build_plugin(PDSPNearlineHeader art::module
              LIBRARIES
              lardataobj::RecoBase
              lardata::headers
              lardataalg::DetectorInfo
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              art::Persistency_Common
              art::Utilities
              canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except::cetlib_except
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

cet_build_plugin(PDSPHitMonitor art::module
              LIBRARIES
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              lardataobj::RawData
              lardata::headers
              lardataalg::DetectorInfo
              nusimdata::SimulationBase
              dunepdlegacy::Services_ChannelMap_PdspChannelMapService_service
              art::Framework_Core
              art::Framework_Principal
              art_root_io::tfile_support
              ROOT::Core
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT::Core ROOT::Hist ROOT::Tree
              BASENAME_ONLY)

install_headers()
install_fhicl()
install_source()
install_scripts()
