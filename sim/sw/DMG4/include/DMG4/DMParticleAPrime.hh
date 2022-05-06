#pragma once

#include <DMParticle.hh>

class DMParticleAPrime : public DMParticle {
private:
    static DMParticleAPrime * theInstance;
    DMParticleAPrime();
    ~DMParticleAPrime();
public:
    static DMParticleAPrime * Definition(G4double MassIn, G4double epsilIn=0.0001);
    static DMParticleAPrime * Definition();
};
