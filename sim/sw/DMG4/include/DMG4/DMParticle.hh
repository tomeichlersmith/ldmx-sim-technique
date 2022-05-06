#pragma once

#include <G4ParticleDefinition.hh>

class DMParticle : public G4ParticleDefinition {
private:
    DMParticle();
    ~DMParticle();
public:
    inline void CalculateLifeTime() {SetPDGLifeTime(CLHEP::hbar_Planck/GetPDGWidth());}
    inline void SetLongLived() {SetPDGLifeTime(1000.);}
};
