
#include <boost/program_options.hpp>

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

  int n_events = vm["num-events"].as<int>();
  bool muons = vm.count("muons") > 0;
  G4double ap_mass, beam_energy;
  if (muons) {
    ap_mass = get(vm, "ap-mass", 200. ) * MeV;
    beam_energy = get(vm, "beam-energy", 100) * GeV;
  } else {
    ap_mass = get(vm, "ap-mass", 100. ) * MeV;
    beam_energy = get(vm, "beam-energy", 4) * GeV;
  }

  framework::config::Parameters model;
  model.addParameter<std::string>("name", "vertex_library");
  model.addParameter("epsilon", 1.);
  model.addParameter("library_path", vm["db-lib"].as<std::string>());
  model.addParameter<std::string>("method", "forward_only");
  model.addParameter("threshold", 0.0);

  // the process accesses the A' mass from the G4 particle
  simcore::darkbrem::G4APrime::APrime(ap_mass/MeV);
  // create the process to do proper initializations
  //    this calculates "common" cross sections as well
  simcore::darkbrem::G4DarkBreMModel db_model(model, muons);

  G4ParticleDefinition* lepton_def = G4Electron::Electron();
  if (muons) lepton_def = G4MuonMinus::MuonMinus();

  // G4Track cleans up its dynamic particle
  auto incident_lepton = new G4DynamicParticle(lepton_def, G4ThreeVector(0.,0.,1.), beam_energy);
  G4Track track(incident_lepton, 0., G4ThreeVector(0.,0.,0.));
  G4Step step;
  step.InitializeStep(&track);
  step.GetPostStepPoint()->SetKineticEnergy(beam_energy);

  int bar_width = 80;
  int pos = 0;
  bool is_redirected = (isatty(STDOUT_FILENO) == 0);
  for (int i_event{0}; i_event < n_events; ++i_event) {

    G4ParticleChange change;
    change.Initialize(track);

    db_model.GenerateChange(change, track, step);

    // extract info from change
    auto dphoton = change.GetSecondary(0);
    auto recoil  = change.GetSecondary(1);

    /*
    trk->GetDynamicParticle()
       ->GetMomentumDirection()
       ->GetKineticEnergy()
       */

    // ~G4ParticleChange cleans up its list of secondaries (I think)

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

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
