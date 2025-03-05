art_make(
  LIB_LIBRARIES
  lardataobj::RecoBase
  larana::OpticalDetector_OpHitFinder
  larsim::MCSTReco
  lardataobj::MCBase
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larana::OpticalDetector_OpDigiProperties_service
  lardata::Utilities
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common
  canvas::canvas
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  cetlib::cetlib 
  cetlib_except::cetlib_except
  fhiclcpp::fhiclcpp
  ROOT::Core ROOT::Tree ROOT::Hist
  MODULE_LIBRARIES larana::OpticalDetector_OpDigiProperties_service
  larana::OpticalDetector_OpHitFinder
  larana::OpticalDetector
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larevt::Filters
  lardata::RecoObjects
  lardataobj::RecoBase
  lardataobj::MCBase
  larreco::RecoAlg
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
  lardata::Utilities
  lardataobj::RawData
  larreco::Calorimetry
  lardataobj::AnalysisBase
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common
  canvas::canvas
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  cetlib::cetlib 
  cetlib_except::cetlib_except
  fhiclcpp::fhiclcpp
  ROOT::Core ROOT::Tree ROOT::Hist
  SERVICE_LIBRARIES
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  CLHEP
  ROOT::Core ROOT::Tree ROOT::Hist
  )

cet_remove_compiler_flags(CXX -pedantic)

install_headers()
install_fhicl()
install_source()
#install_fw( LIST toyWaveform.txt)
#SUBDIRNAME (OpticalDetector )
