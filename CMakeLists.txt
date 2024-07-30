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
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.91.02d01)

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
find_package( dunecore REQUIRED EXPORT )
find_package( duneprototypes REQUIRED EXPORT )
find_package( dunecalib REQUIRED EXPORT )
find_package( nuevdb REQUIRED EXPORT )
find_package( duneopdet REQUIRED EXPORT )
find_package( dunesim REQUIRED EXPORT )
find_package( larcoreobj REQUIRED EXPORT )
find_package( larcore  REQUIRED EXPORT )
find_package( larevt REQUIRED EXPORT )
find_package( lardataobj REQUIRED EXPORT )
find_package( lardataalg REQUIRED EXPORT )
find_package( lardata  REQUIRED EXPORT )
find_package( larsim REQUIRED EXPORT )
find_package( larpandora  REQUIRED EXPORT )
#find_package( pandora REQUIRED )
find_package( larreco REQUIRED EXPORT )
find_package( larrecodnn REQUIRED EXPORT )
find_package( larana REQUIRED EXPORT )
find_package( nutools REQUIRED EXPORT )
find_package( nug4 REQUIRED EXPORT )
find_package( nugen REQUIRED EXPORT )
find_package( GENIE REQUIRED EXPORT )
find_package( nurandom REQUIRED EXPORT )
find_package( art REQUIRED EXPORT )
find_package( art_root_io REQUIRED EXPORT )
find_package( CLHEP REQUIRED EXPORT )
find_package( messagefacility REQUIRED )
find_package( PostgreSQL REQUIRED )
find_package( dunepdlegacy REQUIRED EXPORT )
find_package( artdaq_core REQUIRED EXPORT )
find_package( TRACE REQUIRED )
find_package( ifdh_art REQUIRED )
find_package( ifdhc REQUIRED )
find_package( systematicstools REQUIRED )
find_package( nusystematics REQUIRED )
find_package( geant4reweight REQUIRED )
#find_package( valgrind REQUIRED )
find_package( nlohmann_json REQUIRED )

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
