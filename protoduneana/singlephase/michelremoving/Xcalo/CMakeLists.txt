find_library(CaloTools ../scripts)
cet_make_exec( make_x_correction
               SOURCE protoDUNE_X_calib.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               CaloTools              
               )
cet_make_exec( make_x_correctionAlt
              SOURCE protoDUNE_X_calibAlt.C
              LIBRARIES
              ROOT_BASIC_LIB_LIST
              CaloTools
              )
install_headers()
install_fhicl()
install_source()
install_scripts()

