#pragma once

#include <DMParticle.hh>

class DMParticleAxial : public DMParticle {
private:
    static DMParticleAxial * theInstance;
    DMParticleAxial();
    ~DMParticleAxial();
public:
    static DMParticleAxial * Definition(G4double MassIn, G4double epsilIn=0.0001);
    static DMParticleAxial * Definition();
};
