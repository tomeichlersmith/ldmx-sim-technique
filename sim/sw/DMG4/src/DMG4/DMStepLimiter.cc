#include "DMStepLimiter.hh"
#include "G4TransportationProcessType.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"


DMStepLimiter::DMStepLimiter(const G4String& aName)
  : G4StepLimiter(aName),
    DMMaxStep(DBL_MAX)
{}


DMStepLimiter::~DMStepLimiter()
{}


G4double 
  DMStepLimiter::PostStepGetPhysicalInteractionLength( 
                                       const G4Track& aTrack,
                                       G4double, // previousStepSize
                                       G4ForceCondition* condition  )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  return DMMaxStep;
}

G4VParticleChange*
  DMStepLimiter::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&  )
// Do Nothing
//
{
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}
