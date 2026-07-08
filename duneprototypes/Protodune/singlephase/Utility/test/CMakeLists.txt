# dunetpc/dune/DataPrep/Utility/test/CMakeLists.txt

# David Adams
# June 2018

# Build test for each utility.

include(CetTest)

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR})
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})

cet_test(test_ProtoduneChannelHelper SOURCE test_ProtoduneChannelHelper.cxx
  LIBRARIES
    duneprototypes_Protodune_singlephase_Utility
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    SQLite::SQLite3
    Boost::filesystem
    ROOT::Core ROOT::Hist ROOT::Tree
)

