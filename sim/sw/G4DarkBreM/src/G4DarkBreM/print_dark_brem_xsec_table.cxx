
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
  if (argc < 3) {
    printUsage();
    return 1;
  }

  std::ofstream table_file(argv[1]);

  if (!table_file.is_open()) {
    std::cerr << "File '" << argv[1] << "' was not able to be opened."
              << std::endl;
    return 2;
  }

  framework::config::Parameters model;
  model.addParameter<std::string>("name", "vertex_library");
  model.addParameter<std::string>("library_path", argv[2]);
  model.addParameter<std::string>("method", "forward_only");
  model.addParameter("threshold", 2.0);
  model.addParameter("epsilon", 0.01);

  framework::config::Parameters process;
  process.addParameter("model", model);
  process.addParameter("ap_mass", 10.);
  process.addParameter("enable", true);
  process.addParameter("only_one_per_event", false);
  process.addParameter("cache_xsec", true);
  process.addParameter("global_bias", 1.);

  simcore::darkbrem::G4DarkBremsstrahlung db_process(process);

  table_file << db_process.getCache();

  table_file.close();

  return 0;
}

void printUsage() {
  std::cout << 
    "\n"
    "USAGE: print-dark-brem-xsec-table {xsec_table.csv} {evetn_lib}\n"
    "\n"
    "  ARGUMENTS:\n"
    "     xsec_table.csv  (required) file to print table to\n"
    "     event_lib       (required) directory of darkbrem event library\n"
    << std::endl;
}
