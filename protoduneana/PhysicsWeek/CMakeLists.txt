cet_build_plugin(MichelHitCounter art::module
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
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

cet_build_plugin(VtxScore art::module
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
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)


cet_build_plugin(MichelEventSelection art::module
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
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
