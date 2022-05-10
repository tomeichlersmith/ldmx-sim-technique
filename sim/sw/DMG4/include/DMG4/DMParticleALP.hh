#pragma once

#include <DMParticle.hh>

class DMParticleALP : public DMParticle {
private:
    static DMParticleALP * theInstance;
    DMParticleALP();
    ~DMParticleALP();
public:
    static DMParticleALP * Definition(G4double MassIn, G4double epsilIn);
    static DMParticleALP * Definition();
};
