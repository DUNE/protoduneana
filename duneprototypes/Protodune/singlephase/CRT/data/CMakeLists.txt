include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( LIB_LIBRARIES
                        lardataobj::Simulation
                        nusimdata::SimulationBase
			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io_tfile_support
			ROOT::Core
 			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities canvas::canvas
			cetlib::cetlib cetlib_except::cetlib_except
			lardataobj::AnalysisBase
            DICT_LIBRARIES
                        lardataobj::Simulation
                        nusimdata::SimulationBase
                        larcorealg::Geometry

)

install_headers()
install_source()
