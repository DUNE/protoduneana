include_directories(${HEP_HPC_INCLUDE_DIRS})
art_make(BASENAME_ONLY
  TOOL_LIBRARIES
  art::Framework_Core
  larcorealg::Geometry
  lardata::Utilities
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  hep_hpc_hdf5
  ${HDF5_LIBRARIES}
  larrecodnn::ImagePatternAlgs_Keras
  )

install_headers()
install_fhicl()
install_source()
