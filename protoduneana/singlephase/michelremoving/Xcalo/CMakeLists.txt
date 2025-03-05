find_library(CaloTools ../scripts)
cet_make_exec( NAME make_x_correction
               SOURCE protoDUNE_X_calib.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               CaloTools              
               )
cet_make_exec( make_x_correctionAlt
              SOURCE protoDUNE_X_calibAlt.C
              LIBRARIES
              ROOT::Core
              CaloTools
              )
install_headers()
install_fhicl()
install_source()
install_scripts()

