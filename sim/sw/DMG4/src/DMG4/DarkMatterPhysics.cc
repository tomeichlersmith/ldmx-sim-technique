#include "DarkMatterPhysics.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkZ.hh"
#include "ALP.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkAxials.hh"

#include "DMProcessDMBrem.hh"
#include "DMProcessPrimakoffALP.hh"
#include "DMProcessAnnihilation.hh"

#include "DMParticleAPrime.hh"
#include "DMParticleXBoson.hh"
#include "DMParticleZPrime.hh"
#include "DMParticleALP.hh"
#include "DMParticleScalar.hh"
#include "DMParticleXScalar.hh"
#include "DMParticlePseudoScalar.hh"
#include "DMParticleXPseudoScalar.hh"
#include "DMParticleAxial.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"


DarkMatterPhysics::DarkMatterPhysics() 
: G4VPhysicsConstructor("DarkMatterPhysics")
{
  SetPhysicsType(bUnknown);
  //fMessenger = new DarkMatterPhysicsMessenger();

  if(!DarkMatterPhysicsConfigure()) {
    G4cout << "Dark Matter physics is not properly configured, exiting" << G4endl;
    exit(1);
  }
}

DarkMatterPhysics::DarkMatterPhysics(double Amass,double ratio,double alphaD,double Bias)
: G4VPhysicsConstructor("DarkMatterPhysics")
{
  SetPhysicsType(bUnknown);
  //fMessenger = new DarkMatterPhysicsMessenger();

  if(!DarkMatterPhysicsConfigureWithPars(Amass,ratio,alphaD,Bias)) {
    G4cout << "Dark Matter physics is not properly configured, exiting" << G4endl;
    exit(1);
  }
}



DarkMatterPhysics::~DarkMatterPhysics()
{
  if(myDarkMatter) delete myDarkMatter;
}


void DarkMatterPhysics::ConstructParticle()
{
  // This call to particle definition must be first or at least go before
  // Physics::ConstructProcess()
  DMParticleAPrime::Definition(myDarkMatter->GetMA()*GeV);
  DMParticleXBoson::Definition(myDarkMatter->GetMA()*GeV, myDarkMatter->Getepsil());
  DMParticleZPrime::Definition(myDarkMatter->GetMA()*GeV, myDarkMatter->Getepsil());
  DMParticleALP::Definition(myDarkMatter->GetMA()*GeV, myDarkMatter->Getepsil());
  DMParticleScalar::Definition(myDarkMatter->GetMA()*GeV);
  DMParticleXScalar::Definition(myDarkMatter->GetMA()*GeV, myDarkMatter->Getepsil());
  DMParticlePseudoScalar::Definition(myDarkMatter->GetMA()*GeV);
  DMParticleXPseudoScalar::Definition(myDarkMatter->GetMA()*GeV, myDarkMatter->Getepsil());
  DMParticleAxial::Definition(myDarkMatter->GetMA()*GeV);
}


void DarkMatterPhysics::ConstructProcess()
{
  // Which DM particle?
  G4ParticleDefinition* theDMParticlePtr = 0;
  if(myDarkMatter->GetParentPDGID() == 11) {
    if(myDarkMatter->GetDMType() == 1) {
      theDMParticlePtr = DMParticleAPrime::Definition();
      if(myDarkMatter->Decay()) theDMParticlePtr = DMParticleXBoson::Definition();
    }
    if(myDarkMatter->GetDMType() == 2) {
      theDMParticlePtr = DMParticleScalar::Definition();
      if(myDarkMatter->Decay()) theDMParticlePtr = DMParticleXScalar::Definition();
    }
    if(myDarkMatter->GetDMType() == 3) {
      theDMParticlePtr = DMParticleAxial::Definition();
      //if(myDarkMatter->Decay()) theDMParticlePtr = DMParticleXAxial::Definition();
    }
    if(myDarkMatter->GetDMType() == 4) {
      theDMParticlePtr = DMParticlePseudoScalar::Definition();
      if(myDarkMatter->Decay()) theDMParticlePtr = DMParticleXPseudoScalar::Definition();
    }
  }
  if(myDarkMatter->GetParentPDGID() == -11) { // Annihilation: only vector and scalar for the moment
    if(myDarkMatter->GetDMType() == 1) {
      theDMParticlePtr = DMParticleAPrime::Definition();
    }
    if(myDarkMatter->GetDMType() == 2) {
      theDMParticlePtr = DMParticleScalar::Definition();
    }
  }
  if(myDarkMatter->GetParentPDGID() == 13) {
    theDMParticlePtr = DMParticleAPrime::Definition(); // Fully invisible decay (or almost)
    if(myDarkMatter->Decay()) theDMParticlePtr = DMParticleZPrime::Definition(); // Decay to SM particles
  }
  if(myDarkMatter->GetParentPDGID() == 22) {
    theDMParticlePtr = DMParticleALP::Definition();
  }

  if(!theDMParticlePtr) {G4cout << "DarkMatterPhysics::ConstructProcess: did not manage to determine the DM particle type, exiting" << G4endl; exit(1);}

  myDarkMatter->SetMA(theDMParticlePtr->GetPDGMass()/GeV);
  myDarkMatter->PrepareTable();

  G4PhysicsListHelper * phLHelper = G4PhysicsListHelper::GetPhysicsListHelper();

  phLHelper->DumpOrdingParameterTable();

  // if one need to (re-)associate certain process with a particle, note
  // the following snippet
  //G4ProcessManager * pMgr = Mocktron::Definition()->GetProcessManager();
  //pmanager->RemoveProcess(idxt);
  //pmanager->AddProcess(new G4MonopoleTransportation(fMpl),-1, 0, 0);

  // ... here one can set up the model parameters from external config
  //     sources, internal attributes previously set by messengers, etc

  // ... here the processes asociated with new physics should be registered
  //     as follows
  
  if(myDarkMatter->GetParentPDGID() == 11) {
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Electron::ElectronDefinition() );
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Positron::PositronDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == -11) {
    phLHelper->RegisterProcess( new DMProcessAnnihilation(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Positron::PositronDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == 13) {
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4MuonMinus::MuonMinusDefinition() );
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4MuonPlus::MuonPlusDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == 22) {
    phLHelper->RegisterProcess( new DMProcessPrimakoffALP(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Gamma::GammaDefinition() );
  }
}
