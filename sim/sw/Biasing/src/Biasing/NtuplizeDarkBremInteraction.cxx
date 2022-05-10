#include "Framework/EventProcessor.h"
#include "SimCore/Event/SimParticle.h"

namespace dqm {
class NtuplizeDarkBremInteraction : public framework::Analyzer {
 public:
  NtuplizeDarkBremInteraction(const std::string& n, framework::Process& p)
    : framework::Analyzer(n,p) {}
  virtual void onProcessStart() final override;
  virtual void analyze(const framework::Event& e) final override;
};

void NtuplizeDarkBremInteraction::onProcessStart() {
  getHistoDirectory();
  ntuple_.create("dbint");
  ntuple_.addVar<double>("dbint","x");
  ntuple_.addVar<double>("dbint","y");
  ntuple_.addVar<double>("dbint","z");
  ntuple_.addVar<int>("dbint","incident_pdg");
  ntuple_.addVar<int>("dbint","incident_genstatus");
  ntuple_.addVar<double>("dbint","incident_mass");
  ntuple_.addVar<double>("dbint","incident_energy");
  ntuple_.addVar<double>("dbint","incident_px");
  ntuple_.addVar<double>("dbint","incident_py");
  ntuple_.addVar<double>("dbint","incident_pz");
  ntuple_.addVar<int>("dbint","recoil_pdg");
  ntuple_.addVar<int>("dbint","recoil_genstatus");
  ntuple_.addVar<double>("dbint","recoil_mass");
  ntuple_.addVar<double>("dbint","recoil_energy");
  ntuple_.addVar<double>("dbint","recoil_px");
  ntuple_.addVar<double>("dbint","recoil_py");
  ntuple_.addVar<double>("dbint","recoil_pz");
  ntuple_.addVar<int>("dbint","aprime_pdg");
  ntuple_.addVar<int>("dbint","aprime_genstatus");
  ntuple_.addVar<double>("dbint","aprime_mass");
  ntuple_.addVar<double>("dbint","aprime_energy");
  ntuple_.addVar<double>("dbint","aprime_px");
  ntuple_.addVar<double>("dbint","aprime_py");
  ntuple_.addVar<double>("dbint","aprime_pz");
}

void NtuplizeDarkBremInteraction::analyze(const framework::Event& e) {
  const auto& particle_map{e.getMap<int,ldmx::SimParticle>("SimParticles")};
  const ldmx::SimParticle *incident{nullptr}, *recoil{nullptr}, *aprime{nullptr};
  for (const auto& [track_id, particle] : particle_map) {
    if (particle.getProcessType() == ldmx::SimParticle::ProcessType::eDarkBrem) {
      if (particle.getPdgID() == 622) {
        if (aprime != nullptr) {
          EXCEPTION_RAISE("BadEvent", "Found multiple A' in event.");
        }
        aprime = &particle;
      } else {
        recoil = &particle;
        auto parent_id = particle.getParents().at(0);
        if (particle_map.find(parent_id) != particle_map.end()) {
          incident = &(particle_map.at(parent_id));
        }
      }
    }
  }

  if (recoil == nullptr or aprime == nullptr) {
    EXCEPTION_RAISE("BadEvent","Unable to find both of the products of the dark brem.");
  }

  ntuple_.setVar<double>("x", aprime->getVertex().at(0));
  ntuple_.setVar<double>("y", aprime->getVertex().at(1));
  ntuple_.setVar<double>("z", aprime->getVertex().at(2));
  if (incident != nullptr) {
    ntuple_.setVar<int>("incident_pdg", incident->getPdgID());
    ntuple_.setVar<int>("incident_genstatus", incident->getGenStatus());
    ntuple_.setVar<double>("incident_mass", incident->getMass());
    ntuple_.setVar<double>("incident_energy", incident->getEnergy());
    ntuple_.setVar<double>("incident_px", incident->getMomentum().at(0));
    ntuple_.setVar<double>("incident_py", incident->getMomentum().at(1));
    ntuple_.setVar<double>("incident_pz", incident->getMomentum().at(2));
  } else {
    ntuple_.setVar<int>("incident_pdg", recoil->getPdgID());
    ntuple_.setVar<int>("incident_genstatus", -1);
    ntuple_.setVar<double>("incident_mass", recoil->getMass());
    ntuple_.setVar<double>("incident_energy", recoil->getEnergy()+aprime->getEnergy());
    ntuple_.setVar<double>("incident_px", recoil->getMomentum().at(0)+aprime->getMomentum().at(0));
    ntuple_.setVar<double>("incident_py", recoil->getMomentum().at(1)+aprime->getMomentum().at(1));
    ntuple_.setVar<double>("incident_pz", recoil->getMomentum().at(2)+aprime->getMomentum().at(2));
  }
  ntuple_.setVar<int>("recoil_pdg", recoil->getPdgID());
  ntuple_.setVar<int>("recoil_genstatus", recoil->getGenStatus());
  ntuple_.setVar<double>("recoil_mass", recoil->getMass());
  ntuple_.setVar<double>("recoil_energy", recoil->getEnergy());
  ntuple_.setVar<double>("recoil_px", recoil->getMomentum().at(0));
  ntuple_.setVar<double>("recoil_py", recoil->getMomentum().at(1));
  ntuple_.setVar<double>("recoil_pz", recoil->getMomentum().at(2));
  ntuple_.setVar<int>("aprime_pdg", aprime->getPdgID());
  ntuple_.setVar<int>("aprime_genstatus", aprime->getGenStatus());
  ntuple_.setVar<double>("aprime_mass", aprime->getMass());
  ntuple_.setVar<double>("aprime_energy", aprime->getEnergy());
  ntuple_.setVar<double>("aprime_px", aprime->getMomentum().at(0));
  ntuple_.setVar<double>("aprime_py", aprime->getMomentum().at(1));
  ntuple_.setVar<double>("aprime_pz", aprime->getMomentum().at(2));
}

}

DECLARE_ANALYZER_NS(dqm,NtuplizeDarkBremInteraction);
