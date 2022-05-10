#pragma once

#include <G4VDiscreteProcess.hh>

class DarkMatter;
class G4ParticleDefinition;


class DMProcessAnnihilation : public G4VDiscreteProcess
{
  public:

    DMProcessAnnihilation(DarkMatter* DarkMatterPointer, G4ParticleDefinition* theDMParticlePtrIn, G4double BiasSigmaFactor);

    // Implements final state parameters when the process won.
    virtual G4VParticleChange* PostStepDoIt( const G4Track &, const G4Step & ) override;

    virtual G4double GetMeanFreePath( const G4Track & aTrack,
                                      G4double previousStepSize,
                                      G4ForceCondition * condition ) override;

    virtual G4bool IsApplicable(const G4ParticleDefinition &) override;

  private:

    DarkMatter* myDarkMatter;
    G4ParticleDefinition* theDMParticlePtr;
    G4double BiasSigmaFactor;
};
