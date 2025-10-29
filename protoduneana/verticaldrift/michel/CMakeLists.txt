art_make( MODULE_LIBRARIES
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core ROOT::Tree ROOT::Hist
  art_root_io::TFileService_service
  art::Persistency_Common 
  art::Persistency_Provenance
  art::Utilities 
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardata::RecoObjects
  lardata::ArtDataHelper
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::MCBase
  lardataobj::AnalysisBase
  lardataobj::Simulation
  nusimdata::SimulationBase  

  larsim::Simulation
  larsim::Utils
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_PhotonBackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  
  nug4::ParticleNavigation

  larreco::RecoAlg
  larpandora::LArPandoraInterface
  dunereco::AnaUtils
	
  protoduneana::ProtoDUNEUtilities
)

cet_remove_compiler_flags(CXX -pedantic)

install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
