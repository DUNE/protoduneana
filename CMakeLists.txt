# ======================================================================
#  duneprototypes main build file
# ======================================================================

cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19)


find_package(cetmodules REQUIRED)
project(duneprototypes LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 10.20.07d00)
set(${PROJECT_NAME}_FW_DIR "config_data")

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules)
list(APPEND CMAKE_MODULE_PATH $ENV{ART_DIR}/Modules)


include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
              -Wno-unused-result
)

cet_report_compiler_flags()

find_package( Boost REQUIRED COMPONENTS system )
find_package( art REQUIRED EXPORT )
find_package( art_root_io REQUIRED EXPORT )
find_package( canvas_root_io REQUIRED EXPORT )
find_package(ROOT REQUIRED)
find_package( larcore REQUIRED EXPORT )
find_package( larcorealg REQUIRED EXPORT )
find_package( larcoreobj REQUIRED EXPORT )
find_package( larpandora REQUIRED EXPORT )
find_package( dunecore REQUIRED EXPORT )
find_package( lardataobj REQUIRED EXPORT )
find_package( lardataalg REQUIRED EXPORT )
find_package( lardata REQUIRED EXPORT )
find_package( larana REQUIRED EXPORT )
find_package( larevt REQUIRED EXPORT )
find_package( larsim REQUIRED EXPORT )
find_package( dunecalib REQUIRED EXPORT )
find_package( nusimdata REQUIRED EXPORT )
find_package( nug4 REQUIRED EXPORT )
find_package( dunepdlegacy REQUIRED EXPORT )
find_package( artdaq_core REQUIRED EXPORT )
find_package( TRACE REQUIRED )
find_package( CLHEP REQUIRED )
find_package( dunedaqdataformats REQUIRED EXPORT )
find_package( dunedetdataformats REQUIRED EXPORT )
find_package( PostgreSQL REQUIRED )
find_package( libwda REQUIRED )
find_package( ifdh_art REQUIRED )
find_package( ifdhc REQUIRED )
find_package( ifbeam REQUIRED )
find_package( Geant4 REQUIRED )
find_package( nlohmann_json REQUIRED )
find_package( HDF5 REQUIRED EXPORT )
find_package( HighFive REQUIRED )
find_package( nuevdb REQUIRED )
# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(duneprototypes)

cet_cmake_config()
