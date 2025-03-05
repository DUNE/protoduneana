cet_make_exec( NAME dEdX_calibration 
               SOURCE protoDUNE_dEdx_calib.C
               LIBRARIES
	       ROOT::Core ROOT::Tree ROOT::Hist ROOT::Gpad
               CaloTools
	       cetlib::cetlib 
	       cetlib_except::cetlib_except
	       fhiclcpp::fhiclcpp
               messagefacility::MF_MessageLogger
               )

art_make_library( LIBRARY_NAME protodune_dedx_fit
		  SOURCE ProtoDUNEdEdXFitter.cxx 
		  LIBRARIES 
			    art::Persistency_Provenance
			    cetlib::cetlib
			    fhiclcpp::fhiclcpp
			    ROOT::Core ROOT::Tree ROOT::Hist ROOT::Gpad
                            CaloTools
                            messagefacility::MF_MessageLogger
	    )

cet_make_exec( NAME runProtoDUNEdEdXFit
	       SOURCE runProtoDUNEdEdXFit.c++
	       LIBRARIES protodune_dedx_fit 
	    )
install_scripts()
