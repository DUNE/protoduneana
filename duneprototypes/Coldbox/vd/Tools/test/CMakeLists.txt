# dunetpc/dune/VDColdbox/Tools/test/CMakeLists.txt

# David Adams
# November 2021

# Build test for each tool.

include(CetTest)

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} ${dunedataprep_BINARY_DIR}/fcl)
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_transitive_paths(CONFIG_DATA_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_CONFIG_DATA_DIR} ${PROJECT_BINARY_DIR}/config_data ${dunedataprep_BINARY_DIR}/config_data)

cet_test(test_VDColdboxChannelRanges SOURCE test_VDColdboxChannelRanges.cxx
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

cet_test(test_VDColdboxChannelGroups SOURCE test_VDColdboxChannelGroups.cxx
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

cet_test(test_VDColdboxOnlineChannel SOURCE test_VDColdboxOnlineChannel.cxx
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

cet_test(test_CrpChannelRanges SOURCE test_CrpChannelRanges.cxx
  LIBRARIES
    dunecore::ArtSupport
    dunecore::DuneCommon_Utility
    art::Utilities
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib 
    cetlib_except::cetlib_except
    ROOT::Core ROOT::Hist ROOT::Tree
    SQLite::SQLite3
    Boost::filesystem
)

cet_test(test_CrpChannelGroups SOURCE test_CrpChannelGroups.cxx
  LIBRARIES
    dunecore::ArtSupport
    dunecore::DuneCommon_Utility
    art::Utilities
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    ROOT::Core ROOT::Hist ROOT::Tree
    SQLite::SQLite3
    Boost::filesystem
)

