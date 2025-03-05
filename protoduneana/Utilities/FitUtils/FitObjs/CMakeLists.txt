#cet_make( LIBRARIES canvas::canvas
#		    messagefacility::MF_MessageLogger
#		    MF_Utilities
#		    cetlib::cetlib 
#                   cetlib_except::cetlib_except
#		    ROOT_BASIC_LIB_LIST
        #        NO_DICTIONARY )

      #art_dictionary( DICTIONARY_LIBRARIES dunedataobj_DuneBase )
#art_make(
#  LIB_LIBRARIES
#		    messagefacility::MF_MessageLogger
#		    cetlib::cetlib 
#		    cetlib_except::cetlib_except
#		    ROOT_BASIC_LIB_LIST
#                    Rint
#                    ROOTVecOps
#  DICT_LIBRARIES
#  #  EXCLUDE OpDetDivRec.h OpDetDivRec.cc
#  )

cet_make(LIBRARIES PUBLIC 
		    messagefacility::MF_MessageLogger
		    cetlib::cetlib 
		    cetlib_except::cetlib_except
		    ROOT::Core
        DICT_LIBRARIES
                    ROOT::ROOTVecOps
)


install_headers()
install_source()
