#include "DarkMatterPhysics.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkZ.hh"
#include "ALP.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkAxials.hh"

#include "G4SystemOfUnits.hh"


// BiasSigmaFactor Invisible mode Vector EThresh=35
// 900.  9.e12
// 16.7  8.e8
//  5.   1.75e8
//  2.   3.5e7
// 0.5   1.2e7
// 0.1   3.6e6
// 0.01  1.55e6
// 0.002 1.38e6

// BiasSigmaFactor Invisible mode Scalar EThresh=35
// 16.7  2.3e9

// BiasSigmaFactor Visible mode Vector EThresh=18
// 16.7  3.4e8

bool DarkMatterPhysics::DarkMatterPhysicsConfigure() 
{
  G4double BiasSigmaFactor0 = 9.e8;

  G4double EThresh = 35.; // for sensitivity calculations invisible mode
  //G4double EThresh = 18.; // for sensitivity calculations visible mode
  //G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions

  double DMMass=0.0167;
  myDarkMatter = new DarkPhotons(DMMass, EThresh); // Initialize by default for Pb with eps=0.0001
  //myDarkMatter = new DarkScalars(DMMass, EThresh); // Initialize by default for Pb with eps=0.0001
  //myDarkMatter = new DarkPseudoScalars(DMMass, EThresh); // Initialize by default for Pb with eps=0.0001
  //myDarkMatter = new DarkAxials(DMMass, EThresh); // Initialize by default for Pb with eps=0.0001

  // Below is the code for the visible mode

  //G4double Epsilon = 0.0009; // benchmark value 0.000316

  //myDarkMatter = new DarkPhotons(DMMass, EThresh, 1., 184., 74., 19.25, Epsilon, 2); // Initialize for W
  //myDarkMatter = new DarkScalars(DMMass, EThresh, 1., 184., 74., 19.25, Epsilon, 2); // Initialize for W
  //myDarkMatter = new DarkPhotons(DMMass, EThresh, 1., 207., 82., 11.35, Epsilon, 2); // Initialize for Pb
  //myDarkMatter = new DarkPseudoScalars(DMMass, EThresh, 1., 207., 82., 11.35, Epsilon, 2); // Initialize for Pb
  //myDarkMatter = new ALP(DMMass, EThresh, 1., 207., 82., 11.35, Epsilon, 2); // Initialize for Pb
  //myDarkMatter = new DarkZ(DMMass, EThresh, 1., 207., 82., 11.35, Epsilon, 1); // Initialize for Pb

  // End of code for the visible mode

  BiasSigmaFactor = BiasSigmaFactor0 * 0.0001 * 0.0001 / (myDarkMatter->Getepsil()*myDarkMatter->Getepsil());

  if(!myDarkMatter) return false;
  return true;
}
