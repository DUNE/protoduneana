include_directories("${dunedaqdataformats_DIR}/../../../include")
include_directories("${dunedetdataformats_DIR}/../../../include")
include_directories("${nlohmann_json_DIR}/../../../include")
include_directories("${HighFive_DIR}/../../../include")


cet_build_plugin(PDHDSPSSpillFlagProducer art::module LIBRARIES
  lardataobj::RawData
  dunecore::HDF5Utils_HDF5RawFile2Service_service
  dunecore::dunedaqhdf5utils2
  HDF5::HDF5
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  art::Persistency_Provenance
	art::Utilities
  artdaq_core::artdaq-core_Data
  artdaq_core::artdaq-core_Utilities
  messagefacility::MF_MessageLogger
  ROOT::Core ROOT::Hist ROOT::Tree
  BASENAME_ONLY
)

install_headers()
install_source()
install_fhicl()
