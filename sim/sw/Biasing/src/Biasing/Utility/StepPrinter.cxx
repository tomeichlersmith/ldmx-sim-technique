
#include "Biasing/Utility/StepPrinter.h"

/*~~~~~~~~~~~~*/
/*   Geant4   */
/*~~~~~~~~~~~~*/
#include "G4Step.hh"
#include "G4Event.hh"

#include "G4VBiasingOperator.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"

namespace biasing {
namespace utility {

StepPrinter::StepPrinter(const std::string& name, framework::config::Parameters& parameters)
    : simcore::UserAction(name, parameters) {
  trackID_ = parameters.getParameter<int>("track_id");
}

StepPrinter::~StepPrinter() {}

void StepPrinter::stepping(const G4Step* step) {
  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) return;

  // Get the track associated with this step
  auto track{step->GetTrack()};

  if (auto trackID{track->GetTrackID()}; (trackID_ <= 0) or (trackID != trackID_))
    return;

  // Get the particle name.
  auto particleName{track->GetParticleDefinition()->GetParticleName()};

  // Get the energy of the particle
  auto energy{step->GetPostStepPoint()->GetTotalEnergy()};

  // Get the volume the particle is in.
  auto volume{track->GetVolume()->GetName()};

  // Get the next volume only if there IS a next volume
  G4String nextVolume{"NULL"};
  if (track->GetNextVolume()) {
    nextVolume = track->GetNextVolume()->GetName();
  }

  // Get the region
  auto region{track->GetVolume()->GetLogicalVolume()->GetRegion()->GetName()};

  G4String process{"NA"};
  if (step->GetPostStepPoint()->GetProcessDefinedStep()) {
    process = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  }

  auto bias_op_ptr{G4VBiasingOperator::GetBiasingOperator(track->GetVolume()->GetLogicalVolume())};

  bool first_step_in_vol{
    step->GetPreStepPoint()->GetStepStatus() == G4StepStatus::fGeomBoundary
      or
    track->GetCurrentStepNumber() == 1
  };

  std::cout << " Step " << track->GetCurrentStepNumber() << " {"
            << " Energy: " << energy << " Track ID: " << track->GetTrackID()
            << " Particle currently in: " << volume << " Region: " << region
            << " Next volume: " << nextVolume
            << " Weight: " << track->GetWeight() 
            << " Proc: " << process
            << " StepStatus: " << step->GetPostStepPoint()->GetStepStatus()
            << " BiasOp: " << ((bias_op_ptr)?(bias_op_ptr->GetName()):"NONE")
            << " FirstStepInVol: " << std::boolalpha << first_step_in_vol
            << " Children:";
  for (auto const& track : *(step->GetSecondaryInCurrentStep()))
    std::cout << " " << track->GetParticleDefinition()->GetPDGEncoding();

  std::cout << " }" << std::endl;

  auto print_proc_table = [](const G4ProcessVector& pvec) {
    for (std::size_t i_proc{0}; i_proc < pvec.size(); i_proc++) {
      const G4VProcess* proc{pvec[i_proc]};
      std::cout << "  " << proc->GetProcessName();
      auto bp{dynamic_cast<const G4BiasingProcessInterface*>(proc)};
      if (bp) {
        std::cout << "  biased ";
        auto bop{bp->GetCurrentBiasingOperator()};
        if (bop) {
          std::cout << "with " << bop->GetName();
        } else {
          std::cout << "no operator";
        }
        if (bp->GetIsFirstPostStepGPILInterface(false)) {
          std::cout << " is first GPIL";
        } else {
          std::cout << " is NOT first GPIL";
        }
      }
      std::cout << std::endl;
    } 
  };

  std::cout << "[ StepPrinter ] : process table" << std::endl;
  auto processes{track->GetDefinition()->GetProcessManager()->GetProcessList()};
  if (processes) {
    print_proc_table(*processes);
  } else {
    std::cout << "  EMPTY process vector" << std::endl;
  }
}

}  // namespace utility
}  // namespace biasing

DECLARE_ACTION(biasing::utility, StepPrinter)
