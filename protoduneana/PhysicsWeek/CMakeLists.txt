cet_build_plugin(MichelHitCounter art::module LIBRARIES
              larcorealg::Geometry
              larreco::Calorimetry
              larreco::RecoAlg_PMAlg
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardataobj::RecoBase
              lardata::ArtDataHelper
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art_root_io::tfile_support
              ROOT::Core
              art_root_io::TFileService_service
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT::Core ROOT::Tree ROOT::Hist
              BASENAME_ONLY
)

cet_build_plugin(VtxScore art::module LIBRARIES
              larcorealg::Geometry
              larreco::Calorimetry
              larreco::RecoAlg_PMAlg
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardataobj::RecoBase
              lardata::ArtDataHelper
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art_root_io::tfile_support
              ROOT::Core
              art_root_io::TFileService_service
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT::Core ROOT::Tree ROOT::Hist
              BASENAME_ONLY
)


cet_build_plugin(MichelEventSelection art::module LIBRARIES
              larcorealg::Geometry
              larreco::Calorimetry
              larreco::RecoAlg_PMAlg
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardataobj::RecoBase
              lardata::ArtDataHelper
              lardataalg::DetectorInfo
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art_root_io::tfile_support
	      ROOT::Gpad
              art_root_io::TFileService_service
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT::Core ROOT::Tree ROOT::Hist
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
