
cet_build_plugin(PDHDHitChecker art::module LIBRARIES
                 art::Framework_Core
                 art::Framework_Principal
                 art::Framework_Services_Registry
                 art_root_io::tfile_support
                 larcorealg::Geometry
                 ROOT::Core
                 art_root_io::TFileService_service
                 art::Persistency_Provenance
                 messagefacility::MF_MessageLogger
                 ROOT::Core ROOT::Hist ROOT::Tree
                 BASENAME_ONLY
)

target_include_directories(PDHDHitChecker_module PUBLIC $ENV{LARDATAOBJ_INC} $ENV{LARCOREOBJ_INC} $ENV{LARCOREALG_INC} $ENV{LARCORE_INC})

install_headers()
install_fhicl()
install_source()
