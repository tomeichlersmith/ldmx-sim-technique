#pragma once

#include <DMParticle.hh>

class DMParticlePseudoScalar : public DMParticle {
private:
    static DMParticlePseudoScalar * theInstance;
    DMParticlePseudoScalar();
    ~DMParticlePseudoScalar();
public:
    static DMParticlePseudoScalar * Definition(G4double MassIn, G4double epsilIn=0.0001);
    static DMParticlePseudoScalar * Definition();
};
