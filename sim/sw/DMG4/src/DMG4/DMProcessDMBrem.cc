#include "DMProcessDMBrem.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkZ.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkAxials.hh"

#include "DMParticleAPrime.hh"
#include "DMParticleXBoson.hh"
#include "DMParticleZPrime.hh"
#include "DMParticleScalar.hh"
#include "DMParticlePseudoScalar.hh"
#include "DMParticleAxial.hh"

#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"


DMProcessDMBrem::DMProcessDMBrem(DarkMatter* DarkMatterPointerIn, G4ParticleDefinition* theDMParticlePtrIn,
                                 G4double BiasSigmaFactorIn)
: G4VDiscreteProcess( "DMProcessDMBrem", fUserDefined ),  // fElectromagnetic
  myDarkMatter(DarkMatterPointerIn),
  theDMParticlePtr(theDMParticlePtrIn),
  BiasSigmaFactor(BiasSigmaFactorIn)
{
  SetProcessSubType( 1 ); //fBremsstrahlung? // TODO: verify this
}

G4bool DMProcessDMBrem::IsApplicable(const G4ParticleDefinition & pDef)
{
  if(myDarkMatter->GetParentPDGID() == 11)
    return ("e-" == pDef.GetParticleName() || "e+" == pDef.GetParticleName());
  if(myDarkMatter->GetParentPDGID() == 13)
    return ("mu-" == pDef.GetParticleName() || "mu+" == pDef.GetParticleName());
  //return fModelPtr->IsApplicableToParticle(pDef);
  return false;
}

G4double DMProcessDMBrem::GetMeanFreePath( const G4Track& aTrack,
                                           G4double, /*previousStepSize*/
                                           G4ForceCondition* /*condition*/ )
{
  G4double DensityMat = aTrack.GetMaterial()->GetDensity()/(g/cm3);
  G4double ekin = aTrack.GetKineticEnergy()/GeV;

  if( myDarkMatter->EmissionAllowed(ekin, DensityMat) ) {

    G4double XMeanFreePath = myDarkMatter->GetMeanFreePathFactor()/myDarkMatter->GetSigmaTot(ekin);
    XMeanFreePath /= BiasSigmaFactor;

    return XMeanFreePath;

  }
  return DBL_MAX;
}

G4VParticleChange* DMProcessDMBrem::PostStepDoIt( const G4Track& aTrack,
                                                  const G4Step & aStep )
{
  const G4double incidentE = aTrack.GetTotalEnergy();
  G4ThreeVector incidentDir = aTrack.GetMomentumDirection();
  const G4double incidentKinE = aTrack.GetKineticEnergy();

  G4double XAcc=0., angles[2];
  if(myDarkMatter->GetParentPDGID() == 11)
    if(myDarkMatter->Decay()) {
      XAcc = myDarkMatter->SimulateEmissionWithAngle2(incidentE/GeV, angles);
    } else {
      XAcc = myDarkMatter->SimulateEmission(incidentE/GeV, angles);
    }
  if(myDarkMatter->GetParentPDGID() == 13) {
    //XAcc = myDarkMatter->SimulateEmissionByMuon(incidentE/GeV, angles);  // 2-dim sampling, angles are for the recoil muon
    XAcc = myDarkMatter->SimulateEmissionByMuon2(incidentE/GeV, angles); // 2-step sampling, angles are for the recoil muon
  }

  // Check if it failed? In this case XAcc = 0

  if(XAcc > 0.) myDarkMatter->EmissionSimulated();

  G4double recoilE = incidentKinE - incidentE * XAcc,
           recoilTheta = 0.,
           recoilPhi = 0.;
  G4double DMTheta = angles[0], DMPhi = angles[1];
  if(myDarkMatter->GetParentPDGID() == 13) {
    recoilTheta = angles[0];
    recoilPhi = angles[1];
    DMTheta = 0.;
    DMPhi = 0.;
  }
  G4double DME = incidentE * XAcc;
  G4double DMM = myDarkMatter->GetMA()*GeV;
  G4double DMKinE = DME*DME - DMM*DMM;
  if(DMKinE < 0.) DMKinE = 0.;
  DMKinE = sqrt(DMKinE);

  // Initialize DM direction vector:
  G4ThreeVector DMDirection(0., 0., .1);
  {
    DMDirection.setMag(1.);
    DMDirection.setTheta( DMTheta );
    DMDirection.setPhi( DMPhi );
    DMDirection.rotateUz(incidentDir);
  }
  // Initialize new projectile particle direction vector:
  G4ThreeVector projDirection(0., 0., 1.);
  projDirection.setMag(1.);
  projDirection.setTheta( recoilTheta );
  projDirection.setPhi( recoilPhi );
  projDirection.rotateUz(incidentDir);

  G4DynamicParticle* movingDM = new G4DynamicParticle( theDMParticlePtr,
                                                       DMDirection,
                                                       DMKinE );
  aParticleChange.Initialize( aTrack );

  // Set DM:
  aParticleChange.SetNumberOfSecondaries( 1 );
  aParticleChange.AddSecondary( movingDM );
  // Set projectile changes:
  aParticleChange.ProposeEnergy( recoilE );
  aParticleChange.ProposeMomentumDirection( projDirection );

  std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() 
            << " emitted by " << aTrack.GetDefinition()->GetParticleName()
            << " with energy = " << incidentE/GeV << " DM kinetic energy = " << DMKinE/GeV << std::endl;

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}
