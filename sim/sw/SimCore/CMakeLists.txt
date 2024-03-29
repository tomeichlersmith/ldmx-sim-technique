# Set the minimum version of CMake that's required
cmake_minimum_required(VERSION 3.12)

# Set the project name
project(
  SimCore
  VERSION 2.1.0
  DESCRIPTION "Core classes needed to run the Simulator."
  LANGUAGES CXX)

# Search and configure ROOT
find_package(ROOT CONFIG REQUIRED)

# Search for Boost
find_package(Boost REQUIRED
             COMPONENTS log)

option(BUILD_EVENT_ONLY "Build the event library." ON)
if(BUILD_EVENT_ONLY)

  register_event_object(module_path "SimCore/Event" namespace "ldmx" 
                        class "SimCalorimeterHit" type "collection")
  register_event_object(module_path "SimCore/Event" namespace "ldmx"
                        class "SimTrackerHit" type "collection")
  register_event_object(module_path "SimCore/Event" namespace "ldmx"
                        class "SimParticle" type "map" key "int")

  # Generate the files needed to build the event classes.
  setup_library(module SimCore name Event
                dependencies ROOT::Core
                register_target)

  return()

endif()

# Configure Geant4
setup_geant4_target()

# Get a list of all of the source files.
file(GLOB SRC_FILES CONFIGURE_DEPDENDS ${PROJECT_SOURCE_DIR}/src/SimCore/Persist/[a-zA-Z]*.cxx
                                       ${PROJECT_SOURCE_DIR}/src/SimCore/[a-zA-Z]*.cxx)

# Setup the library
setup_library(module SimCore
              dependencies Geant4::Interface
                           ROOT::Physics
                           Framework::Configure
                           Framework::Framework
                           G4DarkBreM
                           DMG4 UtilsDM DarkMatter
                           Boost::log
                           "${registered_targets}"
              sources ${SRC_FILES})

# Setup the biasing operators library
setup_library(module SimCore
              name BiasOperators
              dependencies SimCore::SimCore)

setup_library(module SimCore
              name Actions
              dependencies SimCore::SimCore)

# Set some target properties
set_target_properties(SimCore
                      PROPERTIES CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

setup_python(package_name ${PYTHON_PACKAGE_NAME}/SimCore)

