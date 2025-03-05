
cet_build_plugin(PDSPmatch art::module LIBRARIES
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larcorealg::Geometry
                        larana::OpticalDetector
                        lardata::Utilities
                        dunepdlegacy::Overlays
                        duneprototypes::Protodune_singlephase_Utility
                        artdaq_core::artdaq-core_Data
                        dunecore::Geometry
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        canvas::canvas
                        art::Persistency_Provenance
                        art::Utilities
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib 
                        cetlib_except::cetlib_except
                        ROOT::Core ROOT::Tree ROOT::Hist
                        dunepdlegacy::rce_dataaccess
                        BASENAME_ONLY
)

cet_build_plugin(PDSPmatchmc art::module LIBRARIES
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larcorealg::Geometry
                        larana::OpticalDetector
                        dunepdlegacy::Overlays
                        lardata::Utilities
                        duneprototypes::Protodune_singlephase_Utility
                        artdaq_core::artdaq-core_Data
                        dunecore::Geometry
                        nusimdata::SimulationBase
                        larsim::MCCheater_PhotonBackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        canvas::canvas
                        art::Persistency_Provenance
                        art::Utilities
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib 
                        cetlib_except::cetlib_except
                        ROOT::Core ROOT::Tree ROOT::Hist
                        dunepdlegacy::rce_dataaccess
                        BASENAME_ONLY
)


install_fhicl()
install_source()
install_scripts()
