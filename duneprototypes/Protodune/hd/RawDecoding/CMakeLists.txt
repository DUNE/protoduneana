include_directories("${dunedaqdataformats_DIR}/../../../include")
include_directories("${dunedetdataformats_DIR}/../../../include")
include_directories("${nlohmann_json_DIR}/../../../include")
include_directories("${HighFive_DIR}/../../../include")

cet_build_plugin(PDHDTimingRawDecoder art::module LIBRARIES
		 dunecore::HDF5Utils_HDF5RawFile3Service_service
		 dunecore::dunedaqhdf5utils2
                 HDF5::HDF5
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support
                 ROOT::Core
                 art_root_io::TFileService_service
                 art::Persistency_Provenance
                 messagefacility::MF_MessageLogger
		 lardataobj::RawData
                 ROOT::Core ROOT::Hist ROOT::Tree
                 BASENAME_ONLY
)

cet_build_plugin(PDHDFEMBFilter art::module LIBRARIES
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support
                 ROOT::Core
                 art::Persistency_Provenance
                 duneprototypes::Protodune_hd_ChannelMap_PD2HDChannelMapService_service
                 messagefacility::MF_MessageLogger
		 lardataobj::RawData
                 ROOT::Core ROOT::Hist ROOT::Tree
                 BASENAME_ONLY
                 dunecore::HDF5Utils_HDF5RawFile3Service_service
                 dunecore::dunedaqhdf5utils2
)


cet_build_plugin(PDHDGroundShakeFilter art::module LIBRARIES
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support
                 ROOT::Core
                 art::Persistency_Provenance
                 duneprototypes_Protodune_hd_ChannelMap_PD2HDChannelMapService_service
                 messagefacility::MF_MessageLogger
         lardataobj::RawData
                 ROOT::Core ROOT::Hist ROOT::Tree
                 BASENAME_ONLY
                 dunecore::HDF5Utils_HDF5RawFile3Service_service
                 dunecore::dunedaqhdf5utils2
)

cet_build_plugin(DAPHNEReaderPDHD art::module LIBRARIES
                 #PDHDReadoutUtils
                 lardataobj::RawData
		 dunecore::HDF5Utils_HDF5RawFile2Service_service
		 dunecore::dunedaqhdf5utils2
                 duneprototypes_Protodune_hd_ChannelMap_DAPHNEChannelMapService_service
                 HDF5::HDF5
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support
                 ROOT::Core
                 art_root_io::TFileService_service
                 art::Persistency_Provenance
                 messagefacility::MF_MessageLogger
                 ROOT::Core ROOT::Hist ROOT::Tree
                 DAPHNEUtils
                 BASENAME_ONLY
)

cet_build_plugin(PDHDTPCReader art::module LIBRARIES
                        lardataobj::RawData
			dunecore::HDF5Utils_HDF5RawFile2Service_service
			dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
			art::Utilities
                        messagefacility::MF_MessageLogger
                        ROOT::Core ROOT::Hist ROOT::Tree
                        BASENAME_ONLY
)

cet_build_plugin(PDHDDataInterfaceWIB3   art::tool LIBRARIES
                        canvas::canvas
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        lardataobj::RawData
                        dunepdlegacy::Overlays
                        artdaq_core::artdaq-core_Data
                        artdaq_core::artdaq-core_Utilities
                        duneprototypes_Protodune_hd_ChannelMap_PD2HDChannelMapService_service
			art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        ROOT::Core ROOT::Hist ROOT::Tree
			dunecore::HDF5Utils_HDF5RawFile2Service_service
			dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
             )

cet_build_plugin(PDHDDataInterfaceWIBEth   art::tool LIBRARIES
                        canvas::canvas
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        lardataobj::RawData
                        dunepdlegacy::Overlays
                        artdaq_core::artdaq-core_Data
                        artdaq_core::artdaq-core_Utilities
                        duneprototypes_Protodune_hd_ChannelMap_PD2HDChannelMapService_service
			art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        ROOT::Core ROOT::Hist ROOT::Tree
			dunecore::HDF5Utils_HDF5RawFile2Service_service
			dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
             )

cet_build_plugin(PDHDDataInterfaceWIBEth3   art::tool LIBRARIES
                        canvas::canvas
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        lardataobj::RawData
                        dunepdlegacy::Overlays
                        artdaq_core::artdaq-core_Data
                        artdaq_core::artdaq-core_Utilities
                        duneprototypes_Protodune_hd_ChannelMap_PD2HDChannelMapService_service
			art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        ROOT::Core ROOT::Hist ROOT::Tree
			dunecore::HDF5Utils_HDF5RawFile3Service_service
			dunecore::dunedaqhdf5utils3
                        HDF5::HDF5
             )


cet_build_plugin(PDHDTriggerReader art::module LIBRARIES
                        lardataobj::RawData
                        dunecore::HDF5Utils_HDF5RawFile2Service_service
                        dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        artdaq_core::artdaq-core_Data
                        artdaq_core::artdaq-core_Utilities
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
			art::Utilities
                        messagefacility::MF_MessageLogger
                        dunecore::dunedaqhdf5utils2
                        ROOT::Core ROOT::Hist ROOT::Tree
                        BASENAME_ONLY
                )

cet_build_plugin(PDHDTriggerReader3 art::module LIBRARIES
                        lardataobj::RawData
                        dunecore::HDF5Utils_HDF5RawFile3Service_service
                        dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        artdaq_core::artdaq-core_Data
                        artdaq_core::artdaq-core_Utilities
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
			art::Utilities
                        messagefacility::MF_MessageLogger
                        dunecore::dunedaqhdf5utils3
                        ROOT::Core ROOT::Hist ROOT::Tree
                        BASENAME_ONLY
                )

cet_build_plugin(PDHDCTBRawDecoder art::module LIBRARIES
		lardataobj::RawData
		dunecore::HDF5Utils_HDF5RawFile3Service_service
		dunecore::dunedaqhdf5utils2
		HDF5::HDF5
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support
		artdaq_core::artdaq-core_Data
		artdaq_core::artdaq-core_Utilities
		ROOT::Core
		art_root_io::TFileService_service
		art::Persistency_Provenance
		art::Utilities
		messagefacility::MF_MessageLogger
		dunecore::dunedaqhdf5utils3
		ROOT::Core ROOT::Hist ROOT::Tree
		BASENAME_ONLY
)

cet_build_plugin(PDHDCherenkovTriggerBitsFilter art::module LIBRARIES
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support
		ROOT::Core
		art::Persistency_Provenance
		duneprototypes_Protodune_hd_ChannelMap_PD2HDChannelMapService_service
		messagefacility::MF_MessageLogger
		lardataobj::RawData
		ROOT::Core ROOT::Hist ROOT::Tree
		BASENAME_ONLY
		dunecore::HDF5Utils_HDF5RawFile3Service_service
		dunecore::dunedaqhdf5utils2
)

cet_make_library(LIBRARY_NAME DAPHNEUtils INTERFACE
                 SOURCE DAPHNEUtils.cxx
                 LIBRARIES
                 lardataobj::RawData
                 ROOT::Core ROOT::Hist ROOT::Tree

)
cet_build_plugin(DAPHNEInterface1   art::tool LIBRARIES
                        canvas::canvas
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        lardataobj::RawData
                        dunepdlegacy::Overlays
                        #artdaq_core::artdaq-core_Data
                        #artdaq_core::artdaq-core_Utilities
                        #duneprototypes_Protodune_hd_ChannelMap_PD2HDChannelMapService_service
                        duneprototypes_Protodune_hd_ChannelMap_DAPHNEChannelMapService_service
                        DAPHNEUtils
			art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        ROOT::Core ROOT::Hist ROOT::Tree
			dunecore::HDF5Utils_HDF5RawFile2Service_service
			dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
             )

cet_build_plugin(DAPHNEInterface2   art::tool LIBRARIES
                        canvas::canvas
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        lardataobj::RawData
                        dunepdlegacy::Overlays
                        #artdaq_core::artdaq-core_Data
                        #artdaq_core::artdaq-core_Utilities
                        #duneprototypes_Protodune_hd_ChannelMap_PD2HDChannelMapService_service
                        duneprototypes_Protodune_hd_ChannelMap_DAPHNEChannelMapService_service
                        DAPHNEUtils
			art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        ROOT::Core ROOT::Hist ROOT::Tree
			dunecore::HDF5Utils_HDF5RawFile3Service_service
			dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
             )

cet_build_plugin(DAPHNEInterface3   art::tool LIBRARIES
                        canvas::canvas
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        lardataobj::RawData
                        dunepdlegacy::Overlays
                        duneprototypes_Protodune_hd_ChannelMap_DAPHNEChannelMapService_service
                        DAPHNEUtils
			art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        ROOT::Core ROOT::Hist ROOT::Tree
			dunecore::HDF5Utils_HDF5RawFile3Service_service
			dunecore::dunedaqhdf5utils2
                        HDF5::HDF5
                        lardataalg::DetectorInfo
                        lardata::DetectorClocksService
             )

#cet_make_library(PDHDReadoutUtils
#                 SOURCE PDHDReadoutUtils.cxx
#                 LIBRARIES INTERFACE
#                 dunecore::dunedaqhdf5utils2
#)


add_subdirectory(fcl)
install_headers()
install_fhicl()
install_source()
install_scripts()
