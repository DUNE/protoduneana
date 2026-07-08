include_directories($ENV{NLOHMANN_JSON_INC})
include_directories("${nlohmann_json_DIR}/../../../include")

install_fhicl()
install_headers()
install_source()

file(GLOB channel_map_file *ChannelMap*.txt)
install_fw( LIST ${channel_map_file} )

art_make(SERVICE_LIBRARIES 
                           art::Framework_Services_Registry
                           art::Framework_Principal
                           art::Framework_Core
                           art::Persistency_Provenance
                           messagefacility::MF_MessageLogger
                           ROOT::Core
 	                   duneprototypes::Protodune_hd_ChannelMap
         LIB_LIBRARIES
                           nlohmann_json::nlohmann_json
)
