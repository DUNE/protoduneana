
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
  lardataobj::RawData
  nusimdata::SimulationBase
  larsim::MCCheater_PhotonBackTrackerService_service           
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  larreco::RecoAlg
  ROOT::Core ROOT::Tree ROOT::Hist
  )
if(DEFINED ENV{HEP_HPC_DIR})
  add_subdirectory(hdf)
endif()

install_headers()
install_fhicl()
install_source()
