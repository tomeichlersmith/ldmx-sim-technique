#pragma once

#include <DMParticle.hh>

class DMParticleXBoson : public DMParticle {
private:
    static DMParticleXBoson * theInstance;
    DMParticleXBoson();
    ~DMParticleXBoson();
public:
    static DMParticleXBoson * Definition(G4double MassIn, G4double epsilIn);
    static DMParticleXBoson * Definition();
};
