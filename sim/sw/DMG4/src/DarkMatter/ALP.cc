// This is a class to simulate conversion of gammas to ALP in matter gammaN -> AN
// Description to follow
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "ALP.hh"
#include "Utils.hh"

#include "G4Electron.hh" // to get CLHEP constants

#include <iostream>


ALP::ALP(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 21;
  ParentPDGID = 22;
  DaughterPDGID = 22;
  std::cout << "Initialized ALP (gamma conversion to ALP) for material density = " << DensityIn << std::endl;
  std::cout << std::endl;
}


ALP::~ALP()
{;}


double ALP::TotalCrossSectionCalc(double E0)
{
  if(E0 < 2.*MA) return 0.; // TODO: the exact threshold is not so big (note that it is present also in the DarkMatter methods)

  double tmin = MA*MA*MA*MA/(4.*E0*E0);
  //double tmax = MA*MA;

  double GaggBench=epsilBench; // Dimensional coupling of ALP in 1/GeV
  double atomFFcoeff=111.0*pow(ZNucl,-1.0/3.0)/Mel;
  double tAtom=1.0/(atomFFcoeff*atomFFcoeff); // atomic coefficient of form-factor
  double tNucl=0.164*pow(ANucl,-2.0/3.0); // nuclear form-factor coefficient in GeV**2
  double LogFactor=log((tNucl+tmin)/(tAtom+tmin))-2.0; // see Note_ALP.pdf and check it
  double sigmaALPtotal=1.0/8.0*GaggBench*GaggBench*alphaEW*ZNucl*ZNucl*LogFactor*GeVtoPb;

  //G4cout << "Total CS calc, E, M, cs = " << E0 << " " << MA << " " << sigmaALPtotal << G4endl;

  return sigmaALPtotal; 
}


double ALP::GetSigmaTot(double E0)
{
  return TotalCrossSectionCalc(E0);
}


double ALP::CrossSectionDSDX(double XEv, double E0)
{
  if(XEv > 0.999) return 1.;
  return 0.;
}


double ALP::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(XEv > 0.999) return 1.;
  return 0.;
}


double ALP::Width()
{
  return 1./(64.*3.1415926)*MA*MA*MA*epsil*epsil;
}
