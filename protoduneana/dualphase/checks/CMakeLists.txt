
cet_build_plugin(       HitNormCheck art::module LIBRARIES
			lardataobj::RecoBase
                        lardata::ArtDataHelper
                        lardata::Utilities
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
                        BASENAME_ONLY
)

cet_build_plugin(	CosmicsdQdx art::module LIBRARIES
			lardataobj::RecoBase
                        lardata::ArtDataHelper
                        lardata::Utilities
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
                        duneprototypes::ProtoDUNEDataUtils
			ProtoDUNEUtilities
                        BASENAME_ONLY
)

cet_build_plugin(	TrackHitInfo art::module LIBRARIES
			lardataobj::RecoBase
                        lardata::ArtDataHelper
                        lardata::Utilities
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
                        duneprototypes::ProtoDUNEDataUtils
			ProtoDUNEUtilities
                        BASENAME_ONLY
)

install_fhicl()
install_source()
