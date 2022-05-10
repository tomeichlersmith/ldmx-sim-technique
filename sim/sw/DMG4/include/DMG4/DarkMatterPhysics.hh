#pragma once

#include <G4VPhysicsConstructor.hh>

class DarkMatter;

class DarkMatterPhysics : public G4VPhysicsConstructor {
  public:

    DarkMatterPhysics();
    ~DarkMatterPhysics();
    bool DarkMatterPhysicsConfigure();

    //A.C. I introduced this method to allow to pass at run-time the A' mass (in GeV), the mChi/mA mass ratio, and the bias
    //The units for Amass must be GeV!!!
    DarkMatterPhysics(double Amass,double ratio,double alphaD,double Bias);
    bool DarkMatterPhysicsConfigureWithPars(double Amass,double ratio,double alphaD,double Bias);


    // Should call initial constructor of particle singletons
    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    DarkMatter* GetDarkMatterPointer() {return myDarkMatter;}
    G4double GetBiasSigmaFactor() {return BiasSigmaFactor;}
  private:
    DarkMatterPhysics(const DarkMatterPhysics &) = delete;
    DarkMatterPhysics & operator=(const DarkMatterPhysics &) = delete;
  private:
    DarkMatter* myDarkMatter;
    G4double BiasSigmaFactor;
};
