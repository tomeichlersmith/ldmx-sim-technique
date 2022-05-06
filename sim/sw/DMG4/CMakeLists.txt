#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(DMG4)

# This below is to compile statements with nullptr (compilers 4.9, 6.3). At some point will not be needed?
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" OFF)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

list(APPEND CMAKE_MODULE_PATH /afs/cern.ch/sw/lcg/external/cmaketools/1.1/modules)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules_cmake_standalone)

#----------------------------------------------------------------------------
# Find GSL
#
find_package(GSL REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include/DMG4 ${GSL_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
file(GLOB sources_DMG4 ${PROJECT_SOURCE_DIR}/src/DMG4/*.cc)
file(GLOB sources_DarkMatter ${PROJECT_SOURCE_DIR}/src/DarkMatter/*.cc)
file(GLOB sources_UtilsDM ${PROJECT_SOURCE_DIR}/src/UtilsDM/*.cc)

add_library(UtilsDM SHARED ${sources_UtilsDM})
set_target_properties(UtilsDM PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

add_library(DMG4 SHARED ${sources_DMG4})
target_link_libraries(DMG4 ${Geant4_LIBRARIES} ${GSL_LIBRARIES} UtilsDM)
set_target_properties(DMG4 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

add_library(DarkMatter SHARED ${sources_DarkMatter})
target_link_libraries(DarkMatter ${Geant4_LIBRARIES} ${GSL_LIBRARIES} UtilsDM)
set_target_properties(DarkMatter PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

#---Install build products---------------------------------------------------
install(TARGETS DMG4 DarkMatter UtilsDM DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN *.hh PATTERN .svn EXCLUDE)
