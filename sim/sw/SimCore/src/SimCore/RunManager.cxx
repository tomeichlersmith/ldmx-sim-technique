/**
 * @file RunManager.cxx
 * @brief Class providing a Geant4 run manager implementation.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "SimCore/RunManager.h"

//-------------//
//   ldmx-sw   //
//-------------//
#include "G4DarkBreM/G4DarkBremsstrahlung.h"  //for process name
#include "SimCore/APrimePhysics.h"
#include "SimCore/DetectorConstruction.h"
#include "SimCore/GammaPhysics.h"
#include "SimCore/PluginFactory.h"
#include "SimCore/PrimaryGeneratorAction.h"
#include "SimCore/USteppingAction.h"
#include "SimCore/UserEventAction.h"
#include "SimCore/UserRunAction.h"
#include "SimCore/UserStackingAction.h"
#include "SimCore/UserTrackingAction.h"

//------------//
//   Geant4   //
//------------//
#include "FTFP_BERT.hh"
#include "G4GDMLParser.hh"
#include "G4GenericBiasingPhysics.hh"
#include "G4ProcessTable.hh"
#include "G4VModularPhysicsList.hh"
#include "G4MuonMinus.hh"

namespace simcore {

RunManager::RunManager(framework::config::Parameters& parameters,
                       ConditionsInterface& ci)
    : conditionsIntf_(ci) {
  parameters_ = parameters;

  // Set whether the ROOT primary generator should use the persisted seed.
  auto rootPrimaryGenUseSeed{
      parameters.getParameter<bool>("rootPrimaryGenUseSeed")};

  // Validate the geometry if specified.
  setUseRootSeed(rootPrimaryGenUseSeed);
}

RunManager::~RunManager() {}

void RunManager::setupPhysics() {
  auto pList{physicsListFactory_.GetReferencePhysList("FTFP_BERT")};
  pList->RegisterPhysics(new GammaPhysics);
  pList->RegisterPhysics(new darkbrem::APrimePhysics(
      parameters_.getParameter<framework::config::Parameters>("dark_brem")));

  auto biasing_operators{
      parameters_.getParameter<std::vector<framework::config::Parameters>>(
          "biasing_operators", {})};
  if (!biasing_operators.empty()) {
    std::cout << "[ RunManager ]: Biasing enabled with "
              << biasing_operators.size() << " operator(s)." << std::endl;

    // create all the biasing operators that will be used
    for (framework::config::Parameters& bop : biasing_operators) {
      simcore::PluginFactory::getInstance().createBiasingOperator(
          bop.getParameter<std::string>("class_name"),
          bop.getParameter<std::string>("instance_name"), bop);
    }

    // Instantiate the constructor used when biasing
    G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();

    // specify which particles are going to be biased
    //  this will put a biasing interface wrapper around *all* processes
    //  associated with these particles
    for (const simcore::XsecBiasingOperator* bop :
         simcore::PluginFactory::getInstance().getBiasingOperators()) {
      std::cout << "[ RunManager ]: Biasing operator '" << bop->GetName()
                << "' set to bias '" << bop->getParticleToBias() << "'" << std::endl;
      biasingPhysics->Bias(bop->getParticleToBias());
    }

    // Register the physics constructor to the physics list:
    pList->RegisterPhysics(biasingPhysics);
  }

  this->SetUserInitialization(pList);
}

void RunManager::Initialize() {
  setupPhysics();

  // This is where the physics lists are told to construct their particles and
  // their processes
  //  They are constructed in order, so it is important to register the biasing
  //  physics *after* any other processes that need to be able to be biased
  G4RunManager::Initialize();

  // Instantiate the primary generator action
  auto primaryGeneratorAction{new PrimaryGeneratorAction(parameters_)};
  SetUserAction(primaryGeneratorAction);

  // Get instances of all G4 actions
  //      also create them in the factory
  auto actions{PluginFactory::getInstance().getActions()};

  // Create all user actions
  auto userActions{
      parameters_.getParameter<std::vector<framework::config::Parameters>>(
          "actions", {})};
  for (auto& userAction : userActions) {
    PluginFactory::getInstance().createAction(
        userAction.getParameter<std::string>("class_name"),
        userAction.getParameter<std::string>("instance_name"), userAction);
  }

  // Register all actions with the G4 engine
  for (const auto& [key, act] : actions) {
    std::visit([this](auto&& arg) { this->SetUserAction(arg); }, act);
  }

  std::cout << "[ RunManager ] Biasing Operators:";
  for (const G4VBiasingOperator* bop : G4VBiasingOperator::GetBiasingOperators()) {
    std::cout << "  " << bop->GetName();
  }
  std::cout << std::endl;

  auto print_proc_table = [](G4ProcessManager* pman) {
    const G4ProcessVector& pvec{*(pman->GetProcessList())};
    for (std::size_t i_proc{0}; i_proc < pvec.size(); i_proc++) {
      const G4VProcess* proc{pvec[i_proc]};
      std::cout << "  " << proc->GetProcessName();
      std::cout << " " << pman->GetProcessOrdering(const_cast<G4VProcess*>(proc), 
          G4ProcessVectorDoItIndex::idxAlongStep);
      std::cout << " " << pman->GetProcessOrdering(const_cast<G4VProcess*>(proc), 
          G4ProcessVectorDoItIndex::idxPostStep);
      std::cout << " " << pman->GetProcessOrdering(const_cast<G4VProcess*>(proc), 
          G4ProcessVectorDoItIndex::idxAtRest);
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

  std::cout << "[ RunManager ] : lepton process tables" << std::endl;
  std::cout << " MUONS" << std::endl;
  print_proc_table(G4MuonMinus::Definition()->GetProcessManager());
  std::cout << " ELECTRONS" << std::endl;
  print_proc_table(G4Electron::Definition()->GetProcessManager());
}

void RunManager::TerminateOneEvent() {
  // have geant4 do its own thing
  G4RunManager::TerminateOneEvent();

  auto reactivate_dark_brem = [](G4ProcessManager* pman) {
    for (std::size_t i_proc{0}; i_proc < pman->GetProcessList()->size(); i_proc++) {
      G4VProcess* p{(*(pman->GetProcessList()))[i_proc]};
      if (p->GetProcessName().contains(G4DarkBremsstrahlung::PROCESS_NAME)) {
        pman->SetProcessActivation(p, true);
        break;
      }
    }
  };

  reactivate_dark_brem(G4MuonMinus::Definition()->GetProcessManager());
  reactivate_dark_brem(G4Electron::Definition()->GetProcessManager());
  if (this->GetVerboseLevel() > 1) {
    std::cout << "[ RunManager ] : "
              << "Reset the dark brem process (if it was activated)."
              << std::endl;
  }
}

DetectorConstruction* RunManager::getDetectorConstruction() {
  return static_cast<DetectorConstruction*>(this->userDetector);
}

}  // namespace simcore
