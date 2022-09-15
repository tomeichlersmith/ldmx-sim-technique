
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
    ("muons",
      "use muons as incident lepton rather than electrons")
    ("min-energy",
      boost::program_options::value<double>(),
      "minimum energy [GeV] to start xsec calculation from (defaults to 2*ap-mass)")
    ("max-energy",
      boost::program_options::value<double>(),
      "maximum energy [GeV] to calculate xsec to (defaults to 100GeV for electrons and 1000GeV for muons)")
    ("energy-step",
      boost::program_options::value<double>()->default_value(10.),
      "difference between adjacent steps in energy [MeV]")
    ("target-Z",
      boost::program_options::value<int>()->default_value(74),
      "Z of target nucleus")
    ("target-A",
      boost::program_options::value<double>()->default_value(183.84),
      "A of target nucleus")
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
  bool muons = vm.count("muons") > 0;

  framework::config::Parameters model;
  model.addParameter<std::string>("name", model_name);
  model.addParameter("epsilon", 1.);
  if (model_name == "vertex_library") {
    model.addParameter<std::string>("library_path", "NOTNEEDED");
    model.addParameter<std::string>("method", "forward_only");
    model.addParameter("threshold", 0.0);
    model.addParameter("load_library", false);
  } else {
    // no other DMG4 parameters at the moment
  }

  framework::config::Parameters process;
  process.addParameter("model", model);
  process.addParameter("ap_mass", ap_mass_MeV);
  process.addParameter("enable", true);
  process.addParameter("only_one_per_event", false);
  process.addParameter("cache_xsec", true);
  process.addParameter("global_bias", 1.);
  process.addParameter("muons", muons);
  process.addParameter("calc_common",false);

  G4double current_energy = 2 * ap_mass_MeV * MeV;
  if (vm.count("min-energy")) {
    current_energy = vm["min-energy"].as<double>() * GeV;
  }
  G4double maximum_energy = (muons ? 1000. : 100) * GeV;
  if (vm.count("max-energy")) {
    maximum_energy = vm["max-energy"].as<double>() * GeV;
  }
  G4double energy_step = vm["energy-step"].as<double>() * MeV;

  // the process accesses the A' mass from the G4 particle
  simcore::darkbrem::G4APrime::APrime(process.getParameter<double>("ap_mass"));
  // create the process to do proper initializations
  //    this calculates "common" cross sections as well
  simcore::darkbrem::G4DarkBremsstrahlung db_process(process);

  int bar_width = 80;
  int pos = 0;
  while (current_energy <= maximum_energy) {
    db_process.getCache().get(current_energy, vm["target-A"].as<double>(), vm["target-Z"].as<int>());
    current_energy += energy_step;
    int old_pos{pos};
    pos = bar_width * current_energy / maximum_energy;
    if (pos != old_pos) {
      std::cout << "[";
      for (int i{0}; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "] " << int(current_energy / maximum_energy * 100.0) << " %\r";
      std::cout.flush();
    }
  }
  std::cout << std::endl;

  table_file << db_process.getCache();

  table_file.close();

  db_process.PrintInfo();

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
