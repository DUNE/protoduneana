art_make( LIBRARY_NAME PDHDSignalToNoiseAlg
  LIB_LIBRARIES
    lardataobj::RecoBase
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    messagefacility::MF_MessageLogger
    ROOT::Core ROOT::Hist

  MODULE_LIBRARIES
    PDHDSignalToNoiseAlg
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    art_root_io::TFileService_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    lardataobj::RecoBase
    lardataobj::RawData
    lardataobj::AnalysisBase
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    lardata::Utilities
    lardataalg::DetectorInfo
    larevt::CalibrationDBI_IOVData
    larevt::CalibrationDBI_Providers
    ROOT::Core ROOT::Tree ROOT::Hist
)

cet_remove_compiler_flags(CXX -pedantic)

install_headers()
install_fhicl()
install_source()
install_scripts()
