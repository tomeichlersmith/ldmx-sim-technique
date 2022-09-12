
//----------------//
//   C++ StdLib   //
//----------------//
#include <fstream>
#include <iostream>

//----------//
//   LDMX   //
//----------//
#include "Framework/Configure/Parameters.h"
#include "G4DarkBreM/G4DarkBremsstrahlung.h"
#include "G4DarkBreM/G4APrime.h"

/**
 * @func printUsage
 *
 * Print how to use this executable to the terminal.
 */
void printUsage();

/**
 * The executable main for printing out the table.
 */
int main(int argc, char* argv[]) {
  if (argc < 2) {
    printUsage();
    return 1;
  }

  std::ofstream table_file(argv[1]);

  if (!table_file.is_open()) {
    std::cerr << "File '" << argv[1] << "' was not able to be opened."
              << std::endl;
    return 2;
  }

  double ap_mass_MeV = 100.;

  framework::config::Parameters model;
  model.addParameter<std::string>("name", "vertex_library");
  model.addParameter<std::string>("library_path", "NOTNEEDED");
  model.addParameter<std::string>("method", "forward_only");
  model.addParameter("threshold", 0.0);
  model.addParameter("epsilon", 0.01);
  model.addParameter("load_library", false);

  framework::config::Parameters process;
  process.addParameter("model", model);
  process.addParameter("ap_mass", ap_mass_MeV);
  process.addParameter("enable", true);
  process.addParameter("only_one_per_event", false);
  process.addParameter("cache_xsec", true);
  process.addParameter("global_bias", 1.);
  process.addParameter("muons", false);

  // the process accesses the A' mass from the G4 particle
  simcore::darkbrem::G4APrime::APrime(process.getParameter<double>("ap_mass"));
  // create the process to do proper initializations
  //    this calculates "common" cross sections as well
  simcore::darkbrem::G4DarkBremsstrahlung db_process(process);

  table_file << db_process.getCache();

  table_file.close();

  db_process.PrintInfo();

  return 0;
}

void printUsage() {
  std::cout << 
    "\n"
    "USAGE: print-dark-brem-xsec-table {xsec_table.csv}\n"
    "\n"
    "  ARGUMENTS:\n"
    "     xsec_table.csv  (required) file to print table to\n"
    << std::endl;
}
