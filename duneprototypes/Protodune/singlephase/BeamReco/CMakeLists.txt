
include_directories("${dunedetdataformats_DIR}/../../../include")
cet_build_plugin( BeamEvent art::module LIBRARIES
                        dunecore::ArtSupport
                        ifbeam::ifbeam
                        ifdh_art::IFBeam_service
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larreco::RecoAlg
                        larreco::Calorimetry
                        larreco::RecoAlg_PMAlg
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        lardata::ArtDataHelper
                        lardata::Utilities
                        lardataalg::UtilitiesHeaders
                        larsim::Simulation
                        lardataobj::Simulation
                        dunecore::DuneObj
                        #dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        larsim::MCCheater_PhotonBackTrackerService_service
                        larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
                        nug4::ParticleNavigation
                        nusimdata::SimulationBase
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib cetlib_except::cetlib_except
                        ROOT::Core ROOT::Hist ROOT::Tree
	  #ProtoDUNEDataUtils
)

find_library(IFBEAM ifbeam)
find_library(IFDH ifdh)
#find_library(IFDH ifdh)
cet_make_exec(
  test_ifbeam
  SOURCE test_ifbeam.c++
  LIBRARIES
    ifdh::ifdh
    ifdh_art::IFBeam_service
    #IFBEAM
    #IFDH
    #wda
    #curl
)

target_include_directories(
  test_ifbeam
  PUBLIC
    $ENV{IFDH_ART_INC}
    $ENV{IFDHC_INC}
    "$ENV{LIBWDA_FQ_DIR}/include"
    "$ENV{IFBEAM_FQ_DIR}/include"
)

target_link_libraries(
   test_ifbeam
   PRIVATE
   ifdh::ifdh
   ifbeam::ifbeam
)


install_headers()
install_fhicl()
install_source()
install_scripts()
