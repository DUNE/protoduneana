# dunetpc/dune/Protodune/singlephase/Tool/test/CMakeLists.txt

# David Adams
# June 2021

# Build test for each tool.
include(CetTest)

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} ${PROJECT_BINARY_DIR}/fcl)
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_transitive_paths(CONFIG_DATA_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_CONFIG_DATA_DIR} ${PROJECT_BINARY_DIR}/config_data)

cet_test(test_PdhdChannelRanges SOURCE test_PdhdChannelRanges.cxx
  LIBRARIES
    dunecore::ArtSupport
    art::Utilities
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    ROOT::Core ROOT::Hist ROOT::Tree
    SQLite::SQLite3
    Boost::filesystem
)

cet_test(test_PdhdChannelGroups SOURCE test_PdhdChannelGroups.cxx
  LIBRARIES
    dunecore::ArtSupport
    art::Utilities
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib 
    cetlib_except::cetlib_except
    ROOT::Core ROOT::Hist ROOT::Tree
    SQLite::SQLite3
    Boost::filesystem
)
