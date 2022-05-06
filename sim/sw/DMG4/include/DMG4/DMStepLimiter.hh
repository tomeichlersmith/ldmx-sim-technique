#ifndef DMStepLimiter_h
#define DMStepLimiter_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"
#include "G4StepLimiter.hh"

class DMStepLimiter : public G4StepLimiter
{
  public:

     DMStepLimiter(const G4String& processName ="StepLimiter");

     virtual ~DMStepLimiter();

     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

     virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            );
                            
     void SetMaxStep(G4double MaxStepIn) {DMMaxStep = MaxStepIn;}

  private:
  
  // hide assignment operator as private 
      DMStepLimiter(DMStepLimiter&);
      DMStepLimiter& operator=(const DMStepLimiter& right);

    G4double DMMaxStep;

};

#endif
