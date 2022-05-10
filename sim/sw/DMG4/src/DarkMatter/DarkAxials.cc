// This is a class to simulate Axials production by electrons in matter eN -> eNA
// Description is in the base class DarkMatter code
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "DarkAxials.hh"
#include "Utils.hh"

#include "G4Electron.hh" // to get CLHEP constants

#include <iostream>

#include "KFactorsAxials.code"


#define  nMALowM 18 // number of MA grid divisions

double TotCSAxialParticle(double MAtest) // CS in GeV^-2 for masses below 1 MeV for epsilon=1
{
  // These are total cross sections of axial-vector DM production in Brem. processes calculated at ETL.
  // The lower X limit of integration is 0.01. It must be the same in the differential cross sections table, then the correct
  // cutoff will be made in sampling.

  double  MMAA[nMALowM] = {0.000001, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.0001, 0.00015, 0.0002,
                           0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009}; // mass of A' in GeV
  double TotCSList[nMALowM] = {39009285714.,390765781.,98066229.,43807779.,24790270.,15972196.,11172093.,8270899.,4153578.,1928757.,1134564.,550112.,335527.,231262.,171841.,134298.,108807.,90558.2};
                              
  return parinv(MAtest, MMAA, TotCSList, nMALowM); // This is to be converted to pb and multiplied by eps^2
}



#include "DiffCS2DInterpAxials.code"


DarkAxials::DarkAxials(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 3;
  ParentPDGID = 11;
  DaughterPDGID = 11;
  std::cout << "Initialized DarkAxials off electrons and positrons for material density = " << DensityIn << std::endl;
  std::cout << std::endl;
}


DarkAxials::~DarkAxials()
{;}


double DarkAxials::TotalCrossSectionCalc(double E0)
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

    sigmaTot= GeVtoPb*(4./3.)*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*fluxAnalytical*beta*log(1./(cutoff*cutoff))/(MA*MA);
    if(sigmaTot < 0.) sigmaTot=0.;

    double KFactor = KfactorAxialsApproximate(MA, E0);

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
      result = (TotCSCut/TotCS) * TotCSAxialParticle(MA)*(ZNucl*ZNucl/(82.*82.))*GeVtoPb*epsilBench*epsilBench; // ETL calculations are made for Pb
                                                                                                                 // The dependency Z^2 is approximate!    
    }
    std::cout << "Total CS calc for masses < 1 MeV, E = " << E0 << "  M = " << MA << " Cut reduction factor = " << TotCSCut/TotCS 
              << " CS = " << result << std::endl;
    return result;

  }
}


double DarkAxials::GetSigmaTot(double E0)
{
  if(MA > 0.001) {
    return GetSigmaTot0(E0);
  } else {
    return TotCSAxialParticle(MA);  }
}


double DarkAxials::CrossSectionDSDX(double XEv, double E0)
{
  if(MA > 0.001) {
    double momentumOfDP=sqrt(XEv*XEv*E0*E0-MA*MA);
    double umaxtilde = -MA*MA*(1.0-XEv)/XEv - Mel*Mel*XEv;
    double numAA = Mel*Mel*XEv*(2.-XEv)*(2.-XEv) - 2*(3-3*XEv+XEv*XEv)*umaxtilde;
    double denAA = 3*XEv*umaxtilde*umaxtilde;
    double numBB = 2*Mel*Mel*(1.-XEv);
    double denBB = umaxtilde*(umaxtilde + Mel*Mel*XEv);
    double sigma = 2*momentumOfDP*(numAA/denAA + numBB/denBB);
    return sigma;
  } else {
    return DsDxBilinearInterpAxials(MA, XEv);  }
}


double DarkAxials::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(MA > 0.001) {
    double Uxtheta = 2.*E0*E0*UThetaEv*XEv + MA*MA*(1.0-XEv)/XEv + Mel*Mel*XEv;
    double AA = Mel*Mel*XEv*XEv/(MA*MA*(1-XEv)) + (1. - XEv + XEv*XEv/2.) / (Uxtheta*Uxtheta);
    double BB = (1. - XEv)*(1. - XEv)*(MA*MA - 4.0*Mel*Mel)/(Uxtheta*Uxtheta*Uxtheta*Uxtheta);
    double CC = MA*MA - Uxtheta*XEv/(1. - XEv) + Mel*Mel*XEv*XEv/(1.0-XEv);
    double sigma = sqrt(XEv*XEv-MA*MA/(E0*E0)) * (AA + BB*CC);
    return sigma;
  } else {
    return DsDxBilinearInterpAxials(MA, XEv);
  }
}


double DarkAxials::Width()
{
  return 1./3.*alphaEW*MA*epsil*epsil*sqrt(1.-4.*Mel*Mel/(MA*MA))*(1.-4.*Mel*Mel/(MA*MA));
}
