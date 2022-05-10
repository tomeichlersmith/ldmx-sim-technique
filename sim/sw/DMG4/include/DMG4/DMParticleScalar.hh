#pragma once

#include <DMParticle.hh>

class DMParticleScalar : public DMParticle {
private:
    static DMParticleScalar * theInstance;
    DMParticleScalar();
    ~DMParticleScalar();
public:
    static DMParticleScalar * Definition(G4double MassIn, G4double epsilIn=0.0001);
    static DMParticleScalar * Definition();
};
