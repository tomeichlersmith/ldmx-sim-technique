#include "globals.hh"

#include "G4ios.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkScalars.hh"
#include "ALP.hh"
#include "DarkZ.hh"

#include "Randomize.hh"


int main() {

  G4double MA = 0.017;
  G4double SigmaNorm = 1.;

  //G4double EThresh = 2.*MA; // for full shape
  //G4double EThresh = 35.; // for sensitivity calculations
  G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions

  DarkMatter* myDarkMatter = new DarkPhotons(MA, EThresh, SigmaNorm); // Initialize by default for Pb with eps=0.0001
  myDarkMatter->PrepareTable();

  double ekin = 100.;
  double angles[2];

  G4cout << G4endl;
  G4cout << "Test of the DarkMatter package: DM emission simulation, energy = " << ekin << " GeV, mass = " << MA << " GeV" << G4endl;
  G4cout << G4endl;

  int NTry=100;
  int ITry;
  for(int i=0; i<NTry; i++) {

    ITry = myDarkMatter->Emission(ekin, 11.35, 1.);
    //double XAcc = myDarkMatter->SimulateEmission(ekin, angles);            // used in invisible mode for electrons
    //double XAcc = myDarkMatter->SimulateEmissionWithAngle(ekin, angles);   // one-step sampling, for electrons, not used by default
    double XAcc = myDarkMatter->SimulateEmissionWithAngle2(ekin, angles);    // two-step sampling, used by default for electrons if decays are enabled
    //double XAcc = myDarkMatter->SimulateEmissionByMuon(ekin, angles);      // one-step sampling, used by default for muons
    //double XAcc = myDarkMatter->SimulateEmissionByMuon2(ekin, angles);     // two-step sampling, for muons, problematic because of dsdx

    if(XAcc > 0.0000001) {
      G4cout << "Emission simulated, X = " << XAcc << " Theta = " << angles[0] << G4endl;
    }

  }
  (void)ITry; // to avoid warning

  G4cout << G4endl;
  G4cout << "Cross section in pb for eps=0.0001 cs = " << myDarkMatter->GetAccumulatedProbability()/((double)NTry) << G4endl;

  return 0;
}
