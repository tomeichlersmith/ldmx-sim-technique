
//----------------//
//   C++ StdLib   //
//----------------//
#include <fstream>
#include <iostream>

#include <boost/program_options.hpp>

//----------//
//   LDMX   //
//----------//
#include "Framework/Configure/Parameters.h"
#include "G4DarkBreM/G4DarkBremsstrahlung.h"
#include "G4DarkBreM/G4APrime.h"

/**
 * The executable main for printing out the table.
 */
int main(int argc, char* argv[]) try {
  boost::program_options::options_description desc(
      "Calculate dark brem cross sections and write them out to a CSV table"
      );
  desc.add_options()
    ("help,h", 
      "produce this help and exit")
    ("output,o", 
      boost::program_options::value<std::string>()->default_value("xsec.csv"), 
      "output file to write xsec to")
    ("model,m", 
      boost::program_options::value<std::string>()->default_value("g4db"), 
      "model to use for calculating")
    ("ap-mass",
      boost::program_options::value<double>()->default_value(100.),
      "mass of A' in MeV")
  ;

  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv)
      .options(desc).run(), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc;
    return 0;
  }


  std::string model_name{vm["model"].as<std::string>()};
  if (model_name == "g4db") {
    model_name = "vertex_library";
  } else if (model_name != "dmg4") {
    std::cerr << "Model '" << model_name << "' is not 'g4db' or 'dmg4'." << std::endl;
    return -1;
  }

  std::ofstream table_file(vm["output"].as<std::string>());
  if (!table_file.is_open()) {
    std::cerr << "File '" << vm["output"].as<std::string>() << "' was not able to be opened."
              << std::endl;
    return 2;
  }

  double ap_mass_MeV = vm["ap-mass"].as<double>();

  framework::config::Parameters model;
  model.addParameter<std::string>("name", model_name);
  model.addParameter<std::string>("library_path", "NOTNEEDED");
  model.addParameter<std::string>("method", "forward_only");
  model.addParameter("threshold", 0.0);
  model.addParameter("epsilon", 1.);
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
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
