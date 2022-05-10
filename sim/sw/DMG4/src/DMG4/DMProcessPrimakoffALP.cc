#include "DMProcessPrimakoffALP.hh"

#include "DarkMatter.hh"
#include "ALP.hh"

#include "DMParticleALP.hh"

#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"


DMProcessPrimakoffALP::DMProcessPrimakoffALP(DarkMatter* DarkMatterPointerIn, G4ParticleDefinition* theDMParticlePtrIn,
                                             G4double BiasSigmaFactorIn)
: G4VDiscreteProcess( "DMProcessPrimakoffALP", fUserDefined ),  // fElectromagnetic
  myDarkMatter(DarkMatterPointerIn),
  theDMParticlePtr(theDMParticlePtrIn),
  BiasSigmaFactor(BiasSigmaFactorIn)
{
  SetProcessSubType( 1 ); //fBremsstrahlung? // TODO: verify this
}

G4bool DMProcessPrimakoffALP::IsApplicable(const G4ParticleDefinition & pDef)
{
  return ("gamma" == pDef.GetParticleName());
}

G4double DMProcessPrimakoffALP::GetMeanFreePath( const G4Track& aTrack,
                                                 G4double, /*previousStepSize*/
                                                 G4ForceCondition* /*condition*/ )
{
  G4double DensityMat = aTrack.GetMaterial()->GetDensity()/(g/cm3);
  G4double ekin = aTrack.GetKineticEnergy()/GeV;
  if( myDarkMatter->EmissionAllowed(ekin, DensityMat) ) {

    G4double XMeanFreePath = myDarkMatter->GetMeanFreePathFactor()/myDarkMatter->GetSigmaTot(ekin);
    XMeanFreePath /= BiasSigmaFactor;

//    std::cout << "DMMeanFreePath = " << XMeanFreePath << std::endl;

    return XMeanFreePath;

  }
  return DBL_MAX;
}

G4VParticleChange* DMProcessPrimakoffALP::PostStepDoIt( const G4Track& aTrack,
                                                        const G4Step & aStep )
{
  const G4double incidentE = aTrack.GetKineticEnergy();
  //const G4double DMMass = theDMParticleAPrimePtr->GetPDGMass();
  G4ThreeVector incidentDir = aTrack.GetMomentumDirection();

  G4double XAcc, angles[2];
  XAcc = myDarkMatter->SimulateEmission(incidentE/GeV, angles);

  // Check if it failed? In this case XAcc = 0

  if(XAcc > 0.001) myDarkMatter->EmissionSimulated();

  G4double DMTheta = angles[0], DMPhi = angles[1];
  G4double DME = incidentE * XAcc;

  // Initialize DM direction vector:
  G4ThreeVector DMDirection(0., 0., .1);
  {
    DMDirection.setMag(1.);
    DMDirection.setTheta( DMTheta );
    DMDirection.setPhi( DMPhi );
    DMDirection.rotateUz(incidentDir);
  }
  
  G4DynamicParticle* movingDM = new G4DynamicParticle( theDMParticlePtr,
                                                       DMDirection,
                                                       DME );
  aParticleChange.Initialize( aTrack );

  // Set DM:
  aParticleChange.SetNumberOfSecondaries( 1 );
  aParticleChange.AddSecondary( movingDM );
  // Kill projectile:
  aParticleChange.ProposeEnergy( 0. );
  aParticleChange.ProposeTrackStatus( fStopAndKill );

  std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() 
            << " emitted by " << aTrack.GetDefinition()->GetParticleName()
            << " with energy = " << incidentE/GeV << " DM energy = " << DME/GeV << std::endl;

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}
