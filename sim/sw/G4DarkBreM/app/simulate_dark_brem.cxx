
#include <boost/program_options.hpp>

#include "TFile.h"
#include "TTree.h"

#include "G4Electron.hh"
#include "G4MuonMinus.hh"

#include "Framework/Configure/Parameters.h"
#include "G4DarkBreM/G4DarkBreMModel.h"
#include "G4DarkBreM/G4APrime.h"

template<typename T>
const T& get(boost::program_options::variables_map& vm, const std::string& name, const T& def) {
  if (vm.count(name)) return vm[name].as<T>();
  else return def;
}

int main(int argc, char* argv[]) try {
  boost::program_options::options_description desc(
      "Run the scaling procedure for the input beam energy and madgraph file"
      );
  desc.add_options()
    ("help,h", "produce this help and exit")
    ("output,o", boost::program_options::value<std::string>()->default_value("scaled.root"),
     "output file to write scaled outgoing kinematics to")
    ("incident-energy,E", boost::program_options::value<double>(),
     "incident beam energy in GeV (defaults to 4 for electrons and 100 for muons)")
    ("num-events,N", boost::program_options::value<int>()->default_value(100),
     "number of dark brems to simulate")
    ("db-lib,L", boost::program_options::value<std::string>(),
     "DB event library to load")
    ("ap-mass",
      boost::program_options::value<double>(),
      "mass of A' in MeV (defaults to 100 for electrons and 200 for muons)")
    ("muons",
      "use muons as incident lepton rather than electrons")
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

  if (vm.count("db-lib") == 0) {
    std::cerr << "ERROR: DB event library not provided." << std::endl;
    return 1;
  }

  int n_events = vm["num-events"].as<int>();
  bool muons = vm.count("muons") > 0;
  double ap_mass, beam_energy, lepton_mass;
  if (muons) {
    ap_mass     = get(vm, "ap-mass"    , 200. ) * MeV;
    beam_energy = get(vm, "beam-energy", 100. );
    lepton_mass = G4MuonMinus::MuonMinus()->GetPDGMass() / GeV;
  } else {
    ap_mass     = get(vm, "ap-mass"    , 100. ) * MeV;
    beam_energy = get(vm, "beam-energy", 4.   );
    lepton_mass = G4Electron::Electron()->GetPDGMass() / GeV;
  }

  framework::config::Parameters model;
  model.addParameter<std::string>("name", "vertex_library");
  model.addParameter("epsilon", 1.);
  model.addParameter("library_path", vm["db-lib"].as<std::string>());
  model.addParameter<std::string>("method", "forward_only");
  model.addParameter("threshold", 0.0);

  TFile f{vm["output"].as<std::string>().c_str(), "recreate"};
  TTree t("dbint","dbint");
  double recoil_energy, recoil_px, recoil_py, recoil_pz;
  t.Branch("recoil_energy", &recoil_energy);
  t.Branch("recoil_px", &recoil_px);
  t.Branch("recoil_py", &recoil_py);
  t.Branch("recoil_pz", &recoil_pz);

  // the process accesses the A' mass from the G4 particle
  simcore::darkbrem::G4APrime::APrime(ap_mass/MeV);
  // create the process to do proper initializations
  //    this calculates "common" cross sections as well
  simcore::darkbrem::G4DarkBreMModel db_model(model, muons);

  int bar_width = 80;
  int pos = 0;
  bool is_redirected = (isatty(STDOUT_FILENO) == 0);
  for (int i_event{0}; i_event < n_events; ++i_event) {
    G4ThreeVector recoil = db_model.scample(beam_energy, lepton_mass);

    recoil_energy = sqrt(recoil.mag2() + lepton_mass*lepton_mass);
    recoil_px = recoil.x();
    recoil_py = recoil.y();
    recoil_pz = recoil.z();

    t.Fill();

    float prog = float(i_event) / n_events;

    int old_pos{pos};
    pos = bar_width * prog;
    if (pos != old_pos) {
      std::cout << "[";
      for (int i{0}; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "] " << int(prog * 100.0) << " %";
      if (is_redirected) std::cout << "\n";
      else std::cout << "\r";
      std::cout.flush();
    }
  }

  t.Write();
  f.Close();

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
