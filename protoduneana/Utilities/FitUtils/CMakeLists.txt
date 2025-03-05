include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories(AFTER $ENV{VALGRIND_INC})

FILE (GLOB src_files *.cxx)

art_make_library( LIBRARY_NAME protodune_fitutils
		  SOURCE ${src_files}
		  LIBRARIES ${ART_FRAMEWORK_CORE}
	       		    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_ROOT_IO_TFILE_SUPPORT} 
			    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			    ${ART_PERSISTENCY_COMMON}
			    canvas::canvas
			    art::Persistency_Provenance
			    cetlib::cetlib cetlib_except::cetlib_except
			    ${ART_UTILITIES}
			    ${MF_MESSAGELOGGER}
			    ${FHICLCPP}
			    ${CETLIB}
			    ROOT::Core ROOT::Tree ROOT::Hist
			    ROOT::RooFit
			    ${ROOT_ROOFITCORE}
			    ${ROOT_ROOFIT}
			    ROOT::RooStats
			    ROOT::HistFactory
			    ${ROOT_FOAM}
			    ${ROOT_GEOM}
			    ${ROOT_MINUIT}
			    ${ROOT_MINUIT2}
			    ${ROOT_CORE}
			    ${ROOT_MATHCORE}
                            ProtoDUNEUtilities
	    )
            #message(WARNING ${ROOT_TREEPLAYER})

cet_make_exec( NAME runProtoDUNEFit
	       SOURCE runProtoDUNEFit.cc
	       LIBRARIES protodune_fitutils
	       ${ROOT_ROOSTATS}
	    )

cet_make_exec( NAME runEventSelection
	       SOURCE runEventSelection.cc
	       LIBRARIES protodune_fitutils
                         ${ROOT_ROOTDATAFRAME}
                         ${ROOT_TREEPLAYER}
                         ${ROOT_ROOTVECOPS}
			 ROOT::Tree
			 ROOT::Core
			 ROOT::ROOTDataFrame
	    )

cet_make_exec( NAME Geant4ReweightSysts 
	       SOURCE Geant4ReweightSysts.cc
	       LIBRARIES protodune_fitutils
	    )

cet_make_exec( NAME runPDSPThinSliceFit
	       SOURCE runPDSPThinSliceFit.c++
	       LIBRARIES protodune_fitutils
	     )

install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(fcl_cfg_files)
add_subdirectory(FitObjs)
