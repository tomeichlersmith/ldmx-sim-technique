#pragma once

#include <DMParticle.hh>

class DMParticleXScalar : public DMParticle {
private:
    static DMParticleXScalar * theInstance;
    DMParticleXScalar();
    ~DMParticleXScalar();
public:
    static DMParticleXScalar * Definition(G4double MassIn, G4double epsilIn);
    static DMParticleXScalar * Definition();
};
