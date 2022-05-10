#pragma once

#include <DMParticle.hh>

class DMParticleXPseudoScalar : public DMParticle {
private:
    static DMParticleXPseudoScalar * theInstance;
    DMParticleXPseudoScalar();
    ~DMParticleXPseudoScalar();
public:
    static DMParticleXPseudoScalar * Definition(G4double MassIn, G4double epsilIn);
    static DMParticleXPseudoScalar * Definition();
};
