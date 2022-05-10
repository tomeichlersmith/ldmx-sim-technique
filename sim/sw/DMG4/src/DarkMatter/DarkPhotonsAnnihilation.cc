/*
 * DarkPhotonsAnnihilation.cc
 *
 *  Created on: Oct 6, 2020
 *      Author: celentan
 *  Fixed: Nov 2, 2020
 */

#include "DarkMatter.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "Utils.hh"

#include "G4Electron.hh" // to get CLHEP constants

#include <iostream>
#include <cmath>


DarkPhotonsAnnihilation::DarkPhotonsAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn,
                                                 double DensityIn, double epsilIn, int IDecayIn, double rIn, double alphaDIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn),
  r(rIn), alphaD(alphaDIn)
{
  DMType = 1; //A.C.
  ParentPDGID = -11;
  DaughterPDGID = 11;
  mChi = MA * r;
  std::cout << "Initialized DarkPhotonsAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
  std::cout << std::endl;
}


DarkPhotonsAnnihilation::~DarkPhotonsAnnihilation()
{;}


//Input: E0, positron energy in GeV
//output: total annihilation cross-section in pbarn.
//Since the framework assumes this method is returning the total cross section per nucleous, for the moment I scale this by Z.
double DarkPhotonsAnnihilation::TotalCrossSectionCalc(double E0)
{
  double ss = 2. * Mel * E0;
  if (sqrt(ss) < 2.*mChi) return 0.;   // A.C. e+e- -> A' -> chi chi can happen also for an A' and chi with large mass,
                                       // i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed
  double qq = sqrt(ss) / 2. * sqrt(1 - 4 * mChi * mChi / (ss));
  double gg = this->Width();

  double sigma = 4 * M_PI * alphaEW * epsil * epsil * alphaD;
  sigma = sigma * qq / sqrt(ss);
  sigma = sigma / ((ss - MA * MA) * (ss - MA * MA) + MA * MA * gg * gg);

  sigma = sigma * (ss - 4. / 3. * qq * qq); // A.C. this is for final state fermions (default)
  //sigma = sigma * (2./3. * qq * qq);      // A.C. this is for final state scalars

  sigma = sigma * GeVtoPb;

  //A.C. correct here for atomic effects
  sigma = sigma * ZNucl;
  return sigma;
}


double DarkPhotonsAnnihilation::GetSigmaTot(double E0) {
  return TotalCrossSectionCalc(E0);
}


bool DarkPhotonsAnnihilation::EmissionAllowed(double E0, double DensityMat) // Different kinematic limit here
{
  if (sqrt(2.*Mel*E0) < 2.*mChi) return false;
  if(E0 < EThresh) return false;
  if(NEmissions) return false; // For G4 DM classes
  if(fabs(DensityMat - Density) > 0.1) return false;
  return true;
}


double DarkPhotonsAnnihilation::CrossSectionDSDX(double XEv, double E0) {
  if (XEv > 0.9999) return 1.;
  return 0.;
}


double DarkPhotonsAnnihilation::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) {
  if (XEv > 0.9999) return 1.;
  return 0.;
}


double DarkPhotonsAnnihilation::Width() {
  double ret;
  ret = MA * epsil * epsil * alphaEW * 1. / 3;
  if (MA/2. > mChi) {
    ret += MA * alphaD * (1 + 2 * mChi * mChi / (MA * MA)) * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 3.;
  }
  return ret;
}


void DarkPhotonsAnnihilation::SetMA(double MAIn) {
  std::cout << "DarkPhotonsAnnihilation::SetMA was called with MAIn = " << MAIn << std::endl;
  MA = MAIn;
  mChi = MA * r;
}
