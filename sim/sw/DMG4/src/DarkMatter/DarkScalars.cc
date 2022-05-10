// This is a class to simulate dark scalars production by electrons in matter eN -> eNA
// Description is in the base class DarkMatter code
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "DarkScalars.hh"
#include "Utils.hh"

#include "G4Electron.hh" // to get CLHEP constants

#include <iostream>

#include "KFactorsScalars.code"


#define  nMALowM 18 // number of MA grid divisions

double TotCSScalarParticle(double MAtest) // CS in GeV^-2 for masses below 1 MeV for epsilon=1
{
  // These are total cross sections of scalar DM production in Brem. processes calculated at ETL.
  // The lower X limit of integration is 0.01. It must be the same in the differential cross sections table, then the correct
  // cutoff will be made in sampling.

  double  MMAA[nMALowM] = {0.000001, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.0001, 0.00015, 0.0002,
                           0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009}; // mass of A' in GeV
  double TotCSList[nMALowM] = {387019.,278547.,219644.,185952.,163112.,146194.,132967.,122242.,99180.3,76215.6,62157.5,45520.8, 35870.9,29518.7,25002.2,21618.9,18987.2,16881.1};
                              
  return parinv(MAtest, MMAA, TotCSList, nMALowM); // This is to be converted to pb and multiplied by eps^2
}



#include "DiffCS2DInterpScalars.code"


DarkScalars::DarkScalars(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 2;
  ParentPDGID = 11;
  DaughterPDGID = 11;
  std::cout << "Initialized DarkScalars off electrons and positrons for material density = " << DensityIn << std::endl;
  std::cout << std::endl;
}


DarkScalars::~DarkScalars()
{;}



double DarkScalars::TotalCrossSectionCalc(double E0)
{
  //double ThetaMaxA;
  //double ThetaMaxEl;
  double sigmaTot;

  if(E0 < 2.*MA) return 0.;

  if(MA > 0.001) { // analytical calculation above 1 MeV

    double tmin = MA*MA*MA*MA/(4.*E0*E0);
    double tmax = MA*MA;

    double aa = 111.*pow(ZNucl,-1./3.)/Mel;
    double d = 0.164*pow(ANucl,-2./3.);
    double ta = pow(1./aa,2.);
    double td = d;
    double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax)
                          - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));

    double beta = sqrt(1. - MA*MA/(E0*E0));
    double cutoff1 = Mel/MA;
    double cutoff2 = MA/E0;
    double cutoff = cutoff2;
    if(cutoff1 > cutoff2) cutoff = cutoff1;

    sigmaTot= GeVtoPb*(2./3.)*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*fluxAnalytical*beta*log(1./(cutoff*cutoff))/(MA*MA);
    if(sigmaTot < 0.) sigmaTot=0.;

    double KFactor = KfactorScalarsApproximate(MA, E0);

    double result = sigmaTot / KFactor; // This K-factor decreases the cross section for MA > ~2 MeV;

    std::cout << "Total CS calc, E = " << E0 << "  M = " << MA << "  KFactor = " << KFactor << " CS = " << result << std::endl;

    return result;
    
  } else {// below MA = 0.001 only ETL tabulated cross sections

    double XMin = 0.01; // to be taken from the table
    if(MA/E0 > XMin) XMin = MA/E0;
    if(XMin > 1.) return 0.;
    int NSteps=10000;
    double StepSize = (1.-XMin)/((double)NSteps);
    double TotCS=0.;
    double TotCSCut=0.;
    for(int ix=0; ix<NSteps; ix++) {
      double xi = XMin + ((double)ix + 0.5)*StepSize;
      if(xi > 1. || xi < XMin) continue;
      TotCS += CrossSectionDSDX(xi, E0);
      if(E0*xi > EThresh) TotCSCut += CrossSectionDSDX(xi, E0);
    }
    double result = 0.;
    if(TotCS > 0.) {
      result = (TotCSCut/TotCS) * TotCSScalarParticle(MA)*(ZNucl*ZNucl/(82.*82.))*GeVtoPb*epsilBench*epsilBench; // ETL calculations are made for Pb
                                                                                                                 // The dependency Z^2 is approximate!    
    }
    std::cout << "Total CS calc for masses < 1 MeV, E = " << E0 << "  M = " << MA << " Cut reduction factor = " << TotCSCut/TotCS 
              << " CS = " << result << std::endl;
    return result;
  }
}


double DarkScalars::GetSigmaTot(double E0)
{
  if(MA > 0.001) {
    return GetSigmaTot0(E0);
  } else {
    return TotCSScalarParticle(MA);
  }
}


double DarkScalars::CrossSectionDSDX(double XEv, double E0)
{
  if(MA > 0.001) {
    double momentumOfDP=sqrt(XEv*XEv*E0*E0-MA*MA);
    double umaxtilde = -MA*MA*(1.0-XEv)/XEv - Mel*Mel*XEv;
    double Numerator = Mel*Mel*(2. -XEv)*(2. -XEv) - 2.*umaxtilde*XEv;
    double Denominator = 3.*umaxtilde*umaxtilde;
    double sigma = momentumOfDP*Numerator/Denominator;
    return sigma;
  } else {
    return DsDxBilinearInterpScalars(MA, XEv);
  }
}


double DarkScalars::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(MA > 0.001) {
    double Uxtheta = 2.*E0*E0*UThetaEv*XEv + MA*MA*(1. - XEv)/XEv + Mel*Mel*XEv;
    double AA = XEv*XEv / (2*Uxtheta*Uxtheta);
    double BB = (1. - XEv)*(1. - XEv)*(MA*MA - 4.0*Mel*Mel)/(Uxtheta*Uxtheta*Uxtheta*Uxtheta);
    double CC = MA*MA - Uxtheta*XEv/(1. - XEv) + Mel*Mel*XEv*XEv/(1. - XEv);
    double sigma = sqrt(XEv*XEv - MA*MA/(E0*E0)) * (AA + BB*CC);
    return sigma;
  } else {
    return DsDxBilinearInterpScalars(MA, XEv);
  }
}


double DarkScalars::Width()
{
  return 1./2.*alphaEW*MA*epsil*epsil*sqrt(1.-4.*Mel*Mel/(MA*MA))*(1.-4.*Mel*Mel/(MA*MA));
}
