cet_make_exec( NAME validate_calibration
               SOURCE protoDUNE_validate_calib.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
	       ROOT::Physics
               fhiclcpp::fhiclcpp
               )

cet_make_exec( NAME calculate_calibration_constants
               SOURCE calconst.C
               LIBRARIES
	       ROOT::Core ROOT::Tree ROOT::Hist ROOT::Gpad
               )

cet_make_library( LIBRARY_NAME CaloTools
                  SOURCE
                  LanGausFit.C CaloUtils.cxx
                  LIBRARIES
                  ROOT::Core ROOT::Tree ROOT::Hist
                  )

cet_make_exec( NAME fitlangaumpv
               SOURCE fitlangaumpv.C
               LIBRARIES
	       ROOT::Core ROOT::Tree ROOT::Hist ROOT::Gpad
               CaloTools
               )

cet_make_exec( NAME murecombfit
               SOURCE murecombfit.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               ROOT_LIBRARIES
	       ROOT::Minuit
	       cetlib::cetlib cetlib_except::cetlib_except
	       fhiclcpp::fhiclcpp
               messagefacility::MF_MessageLogger
               CaloTools
               )

cet_make_exec( NAME plotmurecombfitresult
               SOURCE plotmurecombfitresult.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               ROOT_LIBRARIES
               ROOT::Minuit
               cetlib::cetlib cetlib_except::cetlib_except
               fhiclcpp::fhiclcpp
               messagefacility::MF_MessageLogger
               CaloTools
               )

cet_make_exec( NAME resrangefit
               SOURCE resrangefit.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               ROOT_LIBRARIES
	       ROOT::Minuit
               cetlib::cetlib cetlib_except::cetlib_except
	       fhiclcpp::fhiclcpp
               messagefacility::MF_MessageLogger
               CaloTools
               )

cet_make_exec( NAME ConvertXCorrection
               SOURCE ConvertXCorrection.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               )

cet_make_exec( NAME ConvertYZCorrection
               SOURCE ConvertYZCorrection.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               )

install_scripts()
install_headers()
