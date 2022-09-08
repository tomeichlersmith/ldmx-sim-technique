cmake_minimum_required(VERSION 3.12)

project(G4DarkBreM
  VERSION 1.0.0
  DESCRIPTION "Dark brem simulation integrated within Geant4"
  LANGUAGES CXX)

find_package(ROOT CONFIG REQUIRED)
setup_geant4_target()

setup_library(module G4DarkBreM
  dependencies Geant4::Interface
               DMG4 UtilsDM DarkMatter
               ROOT::Physics
               Framework::Configure
               Framework::Framework)

setup_python(package_name ${PYTHON_PACKAGE_NAME}/G4DarkBreM)

add_executable(print-dark-brem-xsec-table src/G4DarkBreM/print_dark_brem_xsec_table.cxx)
target_link_libraries(print-dark-brem-xsec-table PRIVATE G4DarkBreM::G4DarkBreM)
install(TARGETS print-dark-brem-xsec-table DESTINATION bin)
