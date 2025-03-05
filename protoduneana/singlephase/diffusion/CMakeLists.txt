
art_make(BASENAME_ONLY
  MODULE_LIBRARIES
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
  lardataobj::RecoBase
  lardata::Utilities
  larcorealg::Geometry
  nusimdata::SimulationBase
  larsim::MCCheater_PhotonBackTrackerService_service           
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  larreco::RecoAlg
  ROOT::Core ROOT::Tree ROOT::Hist
  duneprototypes::ProtoDUNEDataUtils
  ProtoDUNEUtilities
  )

cet_remove_compiler_flags(CXX -pedantic)

install_headers()
install_fhicl()
install_source()
