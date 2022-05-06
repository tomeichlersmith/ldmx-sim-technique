#include "globals.hh"

#include "NA64RunAction.hh"
#include "NA64RunActionDMG4.hh"
#include "NA64EventAction.hh"
#include "NA64SteppingActionDMG4.hh"
#include "NA64DetectorConstruction.hh"

#include "QGSP_BERT.hh"
#include "FTFP_BERT.hh"

#include "NA64PrimaryGeneratorAction.hh"

//#include "MockPhysics.hh"
#include "DarkMatterPhysics.hh"
#include "DarkMatter.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
#include "G4PhysListFactory.hh"

#ifdef G4VIS_USE
#include "NA64VisManager.hh"
#include "G4VisExecutive.hh"
#endif

#include "G4ios.hh"


int main(int argc,char** argv) {

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  NA64::DetectorConstruction* mkexp = new NA64::DetectorConstruction;
  mkexp->AddParameter("SetupNumber", 2);
  runManager->SetUserInitialization(mkexp);

//    This is a standard initialization, without DMG4
//  FTFP_BERT* thePL = new FTFP_BERT;
//  //thePL->RegisterPhysics( new G4RadioactiveDecayPhysics );
//  runManager->SetUserInitialization( thePL );


  // ___ Here the "extension" part starts ___
  G4PhysListFactory factory;
  G4VModularPhysicsList * phys = factory.GetReferencePhysList("FTFP_BERT");
  // ^^^ most of the standard physics lists are available by this interface

//  G4PhysicsListHelper * phLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
//  phLHelper->DumpOrdingParameterTable();

  DarkMatterPhysics* myPhysics = new DarkMatterPhysics();
  phys->RegisterPhysics(myPhysics);
  // ^^^ Here the "extension" part ends ^^^
  runManager->SetUserInitialization(phys);  // init phys


#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
//  G4VisManager* visManager = new NA64::VisManager;
  G4VisManager* visManager = new G4VisExecutive;

  visManager->Initialize();
#endif
   
  // UserAction classes
  NA64::RunAction* runAction = new NA64::RunActionDMG4(mkexp, myPhysics->GetDarkMatterPointer(),
                                                              myPhysics->GetBiasSigmaFactor() );
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(new NA64::PrimaryGeneratorAction(mkexp));

  NA64::EventAction* eventAction = new NA64::EventAction(mkexp, runAction);
  runManager->SetUserAction(eventAction);
  //runManager->SetUserAction(new A01TrackingAction);
//  runManager->SetUserAction(new NA64::SteppingActionWithA(mkexp, eventAction));
  runManager->SetUserAction(new NA64::SteppingActionDMG4(mkexp, eventAction));

  // User interactions
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerun.g4mac");    
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
