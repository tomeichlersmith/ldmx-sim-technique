# Setup the project
#
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(NA64simulation)

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
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../modules_cmake_standalone)

#----------------------------------------------------------------------------
# Find GSL
#
find_package(GSL REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/../include ${PROJECT_SOURCE_DIR}/../include/DMG4 ${GSL_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and libraries for this project
#
file(GLOB sources ${PROJECT_BINARY_DIR}/*.cc)
link_directories (${PROJECT_SOURCE_DIR}/../lib)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and other libraries
#
add_executable(mkgeant4 ${sources})
target_link_libraries(mkgeant4 DarkMatter DMG4 UtilsDM ${Geant4_LIBRARIES} ${GSL_LIBRARIES})

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS mkgeant4 DESTINATION bin)
