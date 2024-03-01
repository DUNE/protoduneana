# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================

cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19)

find_package(cetmodules REQUIRED)
project(protoduneana LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.83.01d00)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

# Find our own CMake modules, and tell everyone else where to find them.
cet_cmake_module_directories(Modules BINARY)

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-expansion-to-defined
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_package( Boost REQUIRED COMPONENTS system ) # Needs to be first to cover for a missing transitive
                  # dependency from artdaq_core.
find_package( Geant4 REQUIRED )
find_package(ROOT REQUIRED)

find_package( nufinder REQUIRED )
find_package( dunecore REQUIRED )
find_package( duneprototypes REQUIRED )
find_package( dunecalib REQUIRED )
find_package( nuevdb REQUIRED )
find_package( duneopdet REQUIRED )
find_package( dunesim REQUIRED )
find_package( larcoreobj REQUIRED )
find_package( larcore  REQUIRED )
find_package( larevt REQUIRED )
find_package( lardataobj REQUIRED )
find_package( lardataalg REQUIRED )
find_package( lardata  REQUIRED )
find_package( larsim REQUIRED )
find_package( larpandora  REQUIRED )
#find_package( pandora REQUIRED )
find_package( larreco REQUIRED )
find_package( larrecodnn REQUIRED )
find_package( larana REQUIRED )
find_package( nutools REQUIRED )
find_package( nug4 REQUIRED )
find_package( nugen REQUIRED )
find_package( GENIE REQUIRED )
find_package( nurandom REQUIRED )
find_package( art REQUIRED )
find_package( art_root_io REQUIRED)
find_package( CLHEP REQUIRED )
find_package( messagefacility REQUIRED )
find_package( PostgreSQL REQUIRED )
find_package( dunepdlegacy REQUIRED )
find_package( artdaq_core REQUIRED )
find_package( TRACE REQUIRED )
find_package( ifdh_art REQUIRED )
find_package( ifdhc REQUIRED )
find_package( systematicstools REQUIRED )
find_package( nusystematics REQUIRED )
find_package( geant4reweight REQUIRED )
#find_package( valgrind REQUIRED )

if(DEFINED ENV{HEP_HPC_DIR} )
  find_package( hdf5 REQUIRED )
  find_package( hep_hpc REQUIRED )
endif()

# macros for dictionary and cet_build_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(protoduneana)
add_subdirectory(fcl)

# Unit tests.
#add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
