find_library(CaloTools ../scripts)
cet_make_exec( NAME make_x_correction
               SOURCE protoDUNE_X_calib.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               CaloTools              
               )
install_headers()
install_fhicl()
install_source()
install_scripts()

