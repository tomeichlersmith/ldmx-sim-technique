// This is a class to simulate dark Z' boson production by muons in matter muN -> muNZ'
// Description is in the base class DarkMatter code
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "DarkVector.hh"
#include "Utils.hh"

#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

#include <iostream>
#include "G4ios.hh"


// Auxiliary structures and functions:


// Additional structure holding E0 value to be passed into GSL callbacks.
struct BoundParms {
    DarkVector * this_;
    double E0;
};


// A callback wrapping function for DarkVector::CrossSectionDSDXDTheta()
static double _DarkVectorDsDxDTheta(double x[], size_t dim, void * parms_) {
    BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
    // Forward invocation to target method
    return parms->this_->CrossSectionDSDXDTheta( x[0], x[1], parms->E0 );
}


// Class methods:  ---------------


DarkVector::DarkVector(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 11;
  ParentPDGID = 11;
  DaughterPDGID = 0;
  IApprox = 2; // Approximation: 1 - IWW; 2 - WW
  EMinMessenger = 0.001; // Energy cutoff [GeV], default 0.001
  ThetaMax = 0.1; // Max. angle of Messenger
  std::cout << "Initialized Dark Vector for material density = " << DensityIn << std::endl;
  if(IApprox == 1) std::cout << "Using IWW approximation" << std::endl;
  if(IApprox == 2) std::cout << "Using WW approximation" << std::endl;
  std::cout << "Energy cutoff = " << EMinMessenger << " GeV" << std::endl;
  std::cout << std::endl;
}


DarkVector::~DarkVector()
{;}


double DarkVector::TotalCrossSectionCalc(double E0)
{
  if(IApprox == 2) return TotalCrossSectionCalc_WW2(E0);
  std::cout << "DarkVector: wrong value of IApprox, exiting" << std::endl;
  exit(1);
}


// Below is the 2 - dimensional integration of WW ds/dxdTheta
//
double DarkVector::TotalCrossSectionCalc_WW2(double E0)
{
  if(E0 < 2.*MA) return 0.;

  double Xmin1 = MA/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  if(Xmax1 < Xmin1) return 0.;

  double PrefactorEpsilonAlphaEWE0 = 2.0*epsilBench*epsilBench*alphaEW*alphaEW*alphaEW*E0*E0;

  double xl[2] = { Xmin1, 0.};
  double xu[2] = { Xmax1, ThetaMax};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G;
  BoundParms parms = {this, E0};
  G.f = _DarkVectorDsDxDTheta;
  G.dim = 2;
  G.params = &parms;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  double res, err, sigmaTot;

  // monte_miser: adaptive MC integration

  size_t calls = 5000000;
  gsl_monte_miser_state* stat = gsl_monte_miser_alloc(2);
  gsl_monte_miser_integrate(&G, xl, xu, 2, calls, r, stat, &res, &err);
  gsl_monte_miser_free(stat);
  sigmaTot = GeVtoPb*res*PrefactorEpsilonAlphaEWE0;

  // Other vegas integration methods:
  // gsl_monte_plain : plain MC
  // gsl_monte_vegas

#if 0
  gsl_monte_vegas_state* stat = gsl_monte_vegas_alloc(2);
  gsl_monte_vegas_integrate(&G, xl, xu, 2, 10000, r, s, &res, &err);
  do {
    gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s, &res, &err);
  } while (fabs (gsl_monte_vegas_chisq(stat) - 1.0) > 0.5);
  double sigmaTot = GeVtoPb*res*PrefactorEpsilonAlphaEWE0;
  gsl_monte_vegas_free(stat);
#endif

  gsl_rng_free (r);
  G4cout << "Total CS calc, E = " << E0 << "  M = " << MA << " CS = " << sigmaTot << G4endl;
  return sigmaTot;
}


double DarkVector::GetSigmaTot(double E0)
{
  return GetSigmaTot0(E0);
}


double DarkVector::CrossSectionDSDX(double XEv, double E0)
{
  return CrossSectionDSDX_IWW(XEv, E0);
}


// See ArXiv:1705.01633, second line of Eq.(30)
//
double DarkVector::CrossSectionDSDX_IWW(double XEv, double E0)
{
  if(XEv*E0 <= MA) return 0.;
  double momentumOfDP = sqrt(XEv*XEv*E0*E0-MA*MA);
  double thetamax2 = ThetaMax*ThetaMax;
  double umintilde = -XEv*E0*E0*thetamax2 - MA*MA*(1.0-XEv)/XEv - MParent*MParent*XEv;
  double umaxtilde = -MA*MA*(1.0-XEv)/XEv - MParent*MParent*XEv;
  double NumeratorMax = MParent*MParent*XEv*(-2. + 2.*XEv + XEv*XEv) - 2.*umaxtilde*(3. - 3.*XEv + XEv*XEv);
  double DenominatorMax = 3.*XEv*umaxtilde*umaxtilde;
  double NumeratorMin = MParent*MParent*XEv*(-2. + 2.*XEv + XEv*XEv) - 2.*umintilde*(3. - 3.*XEv + XEv*XEv);
  double DenominatorMin = 3.*XEv*umintilde*umintilde;
  double sigma = momentumOfDP*(NumeratorMax/DenominatorMax - NumeratorMin/DenominatorMin); // with value at thetamax
  //double sigma = momentumOfDP*(NumeratorMax/DenominatorMax); // neglecting value at thetamax
  return sigma;
}


// Below is the IWW formula for the double differential cross-section
// from 1705.01633, see e.g. corresponding second line of Eq.(25)
//
double DarkVector::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(XEv*E0 <= MA) return 0.;
  double Uxtheta = E0*E0*2.*UThetaEv*XEv + MA*MA*(1.0-XEv)/XEv + MParent*MParent*XEv;
  double AA = (1. - XEv + XEv*XEv/2.) / (Uxtheta*Uxtheta);
  double BB = (1. - XEv)*(1. - XEv)*(MA*MA + 2.0*MParent*MParent)/(Uxtheta*Uxtheta*Uxtheta*Uxtheta);
  double CC = MA*MA - Uxtheta*XEv/(1. - XEv);
  double sigma = sqrt(XEv*XEv - MA*MA/(E0*E0)) * (AA + BB*CC);
  return sigma;
}


// WW cross section
//
double DarkVector::CrossSectionDSDXDTheta(double XEv, double ThetaEv, double E0)
{
  if(E0*XEv < EMinMessenger) return 0.;
  double x2 = XEv*XEv;
  double theta2 = ThetaEv*ThetaEv;
  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double MA2 = MA*MA;
  double Mmu2 = MParent*MParent;
  double E02= E0*E0;
  double utilde = -XEv*E02*theta2 - MA2*(1.0-XEv)/XEv - Mmu2*XEv;
  double utilde2=utilde*utilde;
  double ta = 1.0/(aa*aa);
  double td = d;
  double tmax = MA2 + Mmu2;
  double tmin= utilde2/(4.0*E02*(1.0-XEv)*(1.0-XEv));
  // I've calculated ChiWWAnalytical by using mathematica's "Integrate[...]" function
  // and converted the resulted expression to C-like form
  double ChiWWAnalytical = -ZNucl*ZNucl*((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax)
                           - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td)));
  //double ChiWWAnalytical= ZZZ*ZZZ*(log((td+tmin)/(tmin+ta))-2.0); // <--- this approximate expression is valid for MA < 500 MeV and E0=150 GeV
  double Factor1= 2.0*(2.0-2.0*XEv+x2)/(1.0-XEv);
  double Factor2= 4.0*(MA2 + 2.0*Mmu2)/utilde2;
  double Factor3= utilde*XEv + MA2*(1.0-XEv) + Mmu2*x2;
  double AmplZpr2WWVEGAS = Factor1 + Factor2*Factor3;
  // one should multiply the prefactor written below by factor
  // 2.0*epsilon^2*alphaEW^3*E0^2 to get diff_CS_WW_Z' in GeV^(-2)
  // (see e.g. calling VEGAS MC  function)
  double PrefactorWithoutE0EpsilonAlphaEW=sqrt(x2-MA2/E02)*(1.0-XEv)/utilde2;
  double DsDxDthetaWithoutE0EpsilonAlphaEW=sin(ThetaEv)*PrefactorWithoutE0EpsilonAlphaEW*AmplZpr2WWVEGAS*ChiWWAnalytical;

  double ResTemporary;
  if (DsDxDthetaWithoutE0EpsilonAlphaEW < 0.0 ) {
    ResTemporary = 0.0;
  } else {
    ResTemporary = DsDxDthetaWithoutE0EpsilonAlphaEW;
  }
  return ResTemporary;
}


double DarkVector::Width()
{
  const double muMass = G4MuonMinus::MuonMinusDefinition()->GetPDGMass()/GeV;
  const double tauMass = G4TauMinus::TauMinusDefinition()->GetPDGMass()/GeV;
  const double massRatio2 = muMass*muMass/(MA*MA);
  double width          = 0.;
  double nuWidth        = 0.;
  double muWidth        = 0.;
  if (MA < 2.*tauMass) {
    nuWidth = epsil*epsil*alphaEW*(1./3.)*MA;
    if (MA > 2.*muMass) {
      double factor = (1.+2.*massRatio2)*sqrt(1.-4.*massRatio2);
      muWidth = nuWidth*factor;
    }
    width = nuWidth+muWidth; // in GeV
  } else {
    G4cout << "DarkZ: width for the mass above 2 Mtau is not implemented, exiting" << G4endl;
    exit(1);
  }
  return width;
}
