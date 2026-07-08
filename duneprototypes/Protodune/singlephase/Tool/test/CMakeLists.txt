# dunetpc/dune/Protodune/singlephase/Tool/test/CMakeLists.txt

# David Adams
# June 2018

# Build test for each tool.

include(CetTest)

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR})
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_test_env_prepend(FW_SEARCH_PATH ${dunepdlegacy_BINARY_DIR}/files ${PROJECT_BINARY_PATH}/config_files)

cet_test(test_PdspOnlineChannel SOURCE test_PdspOnlineChannel.cxx
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

cet_test(test_IcebergOnlineChannel SOURCE test_IcebergOnlineChannel.cxx
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

cet_test(test_ProtoduneOnlineChannel SOURCE test_ProtoduneOnlineChannel.cxx
  LIBRARIES
    dunecore::ArtSupport
    duneprototypes_Protodune_singlephase_Utility
    art::Utilities
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    ROOT::Core ROOT::Hist ROOT::Tree
    SQLite::SQLite3
    Boost::filesystem
)

cet_test(test_ProtoDuneChannelRanges SOURCE test_ProtoDuneChannelRanges.cxx
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

cet_test(test_ApaChannelRanges SOURCE test_ApaChannelRanges.cxx
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

cet_test(test_HalfApaChannelRanges SOURCE test_HalfApaChannelRanges.cxx
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

cet_test(test_ProtoDuneChannelGroups SOURCE test_ProtoDuneChannelGroups.cxx
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

cet_test(test_ApaChannelGroups SOURCE test_ApaChannelGroups.cxx
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

