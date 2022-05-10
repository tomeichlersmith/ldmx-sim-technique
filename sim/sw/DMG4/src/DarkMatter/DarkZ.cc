// This is a class to simulate dark Z' boson production by muons in matter muN -> muNZ'
// Description is in the base class DarkMatter code
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "DarkZ.hh"
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
#include <gsl/gsl_sf_dilog.h>

#include <iostream>
#include "G4ios.hh"


// Auxiliary structures and functions:


// Additional structure holding E0 value to be passed into GSL callbacks.
struct BoundParms {
    DarkZ * this_;
    double E0;
};


// A callback wrapping function for DarkZ::CrossSectionDSDX_IWW()
static double _DarkZDsDxMuon(double x1, void * parms_) {
    //BoundParms * parms = (BoundParms*) parms_;  // or, equivalently, in C++ style
    BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
    // Forward invocation to target method
    return parms->this_->CrossSectionDSDX_IWW( x1, parms->E0 );
}


// A callback wrapping function for DarkZ::CrossSectionDSDX_WW()
static double _DarkZDsDxMuon_WW(double x1, void * parms_) {
    //BoundParms * parms = (BoundParms*) parms_;  // or, equivalently, in C++ style
    BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
    // Forward invocation to target method
    return parms->this_->CrossSectionDSDX_WW( x1, parms->E0 );
}


// A callback wrapping function for DarkZ::CrossSectionDSDXDTheta()
static double _DarkZDsDxDThetaMuon(double x[], size_t dim, void * parms_) {
    BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
    // Forward invocation to target method
    return parms->this_->CrossSectionDSDXDTheta( x[0], x[1], parms->E0 );
}

// A callback wrapping function for DarkZ::CrossSectionDSDXDpsi()
static double _DarkZDsDxDPsiMuon(double x[], size_t dim, void * parms_) {
    BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
    // Forward invocation to target method
    return parms->this_->CrossSectionDSDXDPSI_WW( x[0], x[1], parms->E0 );
}


// Class methods:  ---------------


DarkZ::DarkZ(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 11;
  ParentPDGID = 13;
  DaughterPDGID = 0;

  IApprox =        2;   // Approximation: 1 - IWW; 2 - WW (default is 2)
  IMethodTotalCS = 1;   // Method for total CS: 1 - ds/dxdTheta; 2 - ds/dxdPsi; 3 - ds/dx (default is 1)
  tMax =       10000.;  // tmax initial; Value 10000. means that tmax = E0*E0 will be taken
  ThetaMax =     0.3;   // Max. angle of Z
  PsiMax =       1.0;   // Max. angle of recoil muon

  std::cout << "Initialized Dark Z boson for material density = " << DensityIn << std::endl;
  if(IApprox == 1) std::cout << "Using IWW approximation" << std::endl;
  if(IApprox == 2) std::cout << "Using WW approximation" << std::endl; 
  if(IMethodTotalCS == 1) std::cout << "ds/dxdTheta is used for total CS" << std::endl;
  if(IMethodTotalCS == 2) std::cout << "ds/dxdPsi is used for total CS" << std::endl;
  if(IMethodTotalCS == 3) std::cout << "ds/dx is used for total CS" << std::endl;
  std::cout << "Energy cutoff = " << EThresh << " GeV" << std::endl;
  std::cout << std::endl;
}


DarkZ::~DarkZ()
{;}


double DarkZ::TotalCrossSectionCalc(double E0)
{
  if(IApprox == 1) return TotalCrossSectionCalc_IWW(E0); 
  if(IApprox == 2) {
    if(IMethodTotalCS == 1) return TotalCrossSectionCalc_WW2(E0); // Integral of ds/dxdTheta
    if(IMethodTotalCS == 2) return TotalCrossSectionCalc_WW3(E0); // Integral of ds/dxdPsi
    if(IMethodTotalCS == 3) return TotalCrossSectionCalc_WW(E0);  // Integral of ds/dx
    std::cout << "DarkZ: wrong value of IMethodTotalCS, exiting" << std::endl;
    exit(1);
  } else {
    std::cout << "DarkZ: wrong value of IApprox, exiting" << std::endl;
  }
  exit(1);
}


// We integrate below the differential cross-section over x and t (see e.g. second line of Eq.(30) in 1705.01633,
// where x=E_Z'/E_0, t is a transfer momentum squared
//
double DarkZ::TotalCrossSectionCalc_IWW(double E0)
{
  if(E0 < 2.*MA) return 0.;

  gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000);
  double result1, error1;
  double tmin = MA*MA*MA*MA/(4.*E0*E0);
  //double tmax = MA*MA+Mmu*Mmu;
  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;
  double Xmin1=MA/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;

  gsl_function F1;
  BoundParms parms = { this, E0 };
  F1.function = _DarkZDsDxMuon;
  F1.params = &parms;

  gsl_integration_qags (&F1, Xmin1, Xmax1, 0, 1e-7, 1000, w1, &result1, &error1);

  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double ta = pow(1./aa,2.);
  double td = d;
  double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax)
                          - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));
  double IntDsDx = result1;
  gsl_integration_workspace_free (w1);

  double PrefactorMuonZTotCS= 2.0*epsilBench*epsilBench*alphaEW*alphaEW*alphaEW/E0;

  double sigmaTot= GeVtoPb*PrefactorMuonZTotCS*fluxAnalytical*IntDsDx;
  G4cout << "Total CS calc, E = " << E0 << "  M = " << MA << " CS = " << sigmaTot << G4endl;
  return sigmaTot;
}


// Below is the integration of WW ds/dx
//
double DarkZ::TotalCrossSectionCalc_WW(double E0)
{
  double sigmaTot;

  //if(E0 < 2.*MA) return 0.;
  //if(E0 < MA*99.) return 0.; // Some approximations probably don't work for smaller energy

  //to switch off default error handler, store old error handler in old_handler:
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

  gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000);
  double result1, error1;
  double Xmin1 = MA/E0;
  if(EThresh/E0 > Xmin1) Xmin1 = EThresh/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  if(Xmax1 < Xmin1) return 0.;

  gsl_function F1;
  BoundParms parms = {this, E0};
  F1.function = _DarkZDsDxMuon_WW;
  F1.params = &parms;

  //gsl_integration_qags (&F1, Xmin1, Xmax1, 0, 1e-7, 1000, w1, &result1, &error1);
  double relerr=1.0e-7;   //initial error tolerance (relative error)
  int status=1;
  while(status) {
    status=gsl_integration_qags (&F1, Xmin1, Xmax1, 0, relerr, 1000, w1, &result1, &error1);
    relerr *= 1.2;
    if(status) G4cout << "Increased tolerance=" << relerr << G4endl;
  }
  //if integration routine returns error code, integration is repeated
  //using increased error tolerance, message is printed out
  gsl_set_error_handler(old_handler); //reset error handler (might be unneccessary.)

  double IntDsDx = result1;
  gsl_integration_workspace_free (w1);

  //double PrefactorMuonZTotCS = 2.0*epsilBench*epsilBench*alphaEW*alphaEW*alphaEW/E0; // old prefactor
  double PrefactorMuonZTotCS = epsilBench*epsilBench*alphaEW*alphaEW*alphaEW*ZNucl*ZNucl;
  sigmaTot= GeVtoPb*PrefactorMuonZTotCS*IntDsDx;

  G4cout << "Total CS calc, E = " << E0 << "  M = " << MA << " CS = " << sigmaTot << G4endl;

  return sigmaTot;
}


// Below is the 2 - dimensional integration of WW ds/dxdTheta
//
double DarkZ::TotalCrossSectionCalc_WW2(double E0)
{
  if(E0 < 2.*MA) return 0.;

  double Xmin1 = MA/E0;
  if(EThresh/E0 > Xmin1) Xmin1 = EThresh/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  if(Xmax1 < Xmin1) return 0.;

  double PrefactorEpsilonAlphaEWE0 = 2.0*epsilBench*epsilBench*alphaEW*alphaEW*alphaEW*E0*E0;

  double xl[2] = { Xmin1, 0.};
  double xu[2] = { Xmax1, ThetaMax};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G;
  BoundParms parms = {this, E0};
  G.f = _DarkZDsDxDThetaMuon;
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

// Below is the 2 - dimensional integration of WW ds/dxdPsi
//
double DarkZ::TotalCrossSectionCalc_WW3(double E0)
{
  if(E0 < 2.*MA) return 0.;

  double Xmin1 = MA/E0;
  if(EThresh/E0 > Xmin1) Xmin1 = EThresh/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  if(Xmax1 < Xmin1) return 0.;

  double PrefactorEpsilonAlphaEWE0 = 8.0*epsilBench*epsilBench*alphaEW*alphaEW*alphaEW*E0*E0*ZNucl*ZNucl;

  double xl[2] = { Xmin1, 0.};
  double xu[2] = { Xmax1, PsiMax};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G;
  BoundParms parms = {this, E0};
  G.f = _DarkZDsDxDPsiMuon;
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


double DarkZ::GetSigmaTot(double E0)
{
  return GetSigmaTot0(E0);
}


double DarkZ::CrossSectionDSDX(double XEv, double E0)
{
  if(IApprox == 1) return CrossSectionDSDX_IWW(XEv, E0);
  if(IApprox == 2) return CrossSectionDSDX_WW(XEv, E0);
  std::cout << "DarkZ: wrong value of IApprox, exiting" << std::endl;
  exit(1);
}


// See ArXiv:1705.01633, second line of Eq.(30)
//
double DarkZ::CrossSectionDSDX_IWW(double XEv, double E0)
{
  if(XEv*E0 <= MA) return 0.;
  double momentumOfDP = sqrt(XEv*XEv*E0*E0-MA*MA);
  double thetamax2 = ThetaMax*ThetaMax;
  double umintilde = -XEv*E0*E0*thetamax2 - MA*MA*(1.0-XEv)/XEv - Mmu*Mmu*XEv;
  double umaxtilde = -MA*MA*(1.0-XEv)/XEv - Mmu*Mmu*XEv;
  double NumeratorMax = Mmu*Mmu*XEv*(-2. + 2.*XEv + XEv*XEv) - 2.*umaxtilde*(3. - 3.*XEv + XEv*XEv);
  double DenominatorMax = 3.*XEv*umaxtilde*umaxtilde;
  double NumeratorMin = Mmu*Mmu*XEv*(-2. + 2.*XEv + XEv*XEv) - 2.*umintilde*(3. - 3.*XEv + XEv*XEv);
  double DenominatorMin = 3.*XEv*umintilde*umintilde;
  double sigma = momentumOfDP*(NumeratorMax/DenominatorMax - NumeratorMin/DenominatorMin); // with value at thetamax
  //double sigma = momentumOfDP*(NumeratorMax/DenominatorMax); // neglecting value at thetamax
  return sigma;
}


#if 0
// Below is cross section obtained by integration of ds/dxdpsi with additional simplification of flux.
// for references see e.g. https://www.overleaf.com/8686251169nmvbnqfszfxs.
// It is unstable, works only for very small M/E
//
double DarkZ::CrossSectionDSDX_WW(double XEv, double E0)
{
  double TempThetaMax= 0.07; //is the typical maximum angle from NA64mu design (We should disciss that value)
  double TempThetaMax2=TempThetaMax*TempThetaMax;
  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double MA2 = MA*MA;
  double Mmu2 = Mmu*Mmu;
  double E02 = E0*E0;
  double umax = -MA2*(1.0-XEv)/XEv - Mmu2*XEv;
  double x2 = XEv*XEv;
  double umax2 = umax*umax;
  double umax3= umax2*umax;
  double umin = -XEv*E02*TempThetaMax2 + umax;
  double umin2 = umin*umin;
  double umin3 = umin2*umin;
  double betad2 = 4.0*E02*(1.0-XEv)*(1.0-XEv)*d/7.38905609893; // here e^2 = (2.71828182846)^2 = 7.38905609893
  double betaa2 = 4.0*E02*(1.0-XEv)*(1.0-XEv)/(aa*aa);
  double I2max = -1.0/umax*log(betad2/(betaa2 + umax2)) - 2.0/sqrt(betaa2)*atan(umax/sqrt(betaa2));
  double I2min = -1.0/umin*log(betad2/(betaa2 + umin2)) - 2.0/sqrt(betaa2)*atan(umin/sqrt(betaa2));
  double I2 = I2max - I2min;
  double I3max = -1.0/(2.0*umax2)*log(betad2/(betaa2 + umax2))-1.0/(2.0*betaa2)*log(umax2/(betaa2 + umax2));
  double I3min = -1.0/(2.0*umin2)*log(betad2/(betaa2 + umin2))-1.0/(2.0*betaa2)*log(umin2/(betaa2 + umin2));
  double I3 = I3max - I3min;
  double I4max = -1.0/(3.0*umax3)*log(betad2/(betaa2 + umax2)) + 2.0/(3.0*umax*betaa2) + 2.0/(3.0*betaa2*sqrt(betaa2))*atan(umax/sqrt(betaa2));
  double I4min = -1.0/(3.0*umin3)*log(betad2/(betaa2 + umin2)) + 2.0/(3.0*umin*betaa2) + 2.0/(3.0*betaa2*sqrt(betaa2))*atan(umin/sqrt(betaa2));
  double I4 = I4max - I4min;
  double AtWWAnalyticalIntegral = 2.0*((2.0 - 2.0*XEv + x2)/(1.0 - XEv))*I2 + 4.0*(MA2 + 2.0*Mmu2)*(I3*XEv + I4*(MA2*(1.0 - XEv) + Mmu2*x2));
  double PrefactorAnalyticalIntegral = sqrt(x2 - MA2/E02)*(1.0 - XEv)/XEv;

  //For absolute value of the diff. cs one should multiply dsdxWW by FactorForTrueCS written below
  //double epsilon=0.001;
  //double alphaQED=1.0/137.0;
  //double FactorForTrueCS = epsilon*epsilon*alphaQED*alphaQED*alphaQED*ZZ*ZZ;

  double dsdxWW;
  if (AtWWAnalyticalIntegral < 0.0 ) {
    dsdxWW =0.0;
  } else {
    dsdxWW = AtWWAnalyticalIntegral*PrefactorAnalyticalIntegral;
  }
  return dsdxWW;
}
#endif


// Below is cross section obtained by integration of ds/dxdTheta
//
double DarkZ::CrossSectionDSDX_WW(double XEv, double E0)
{
  // return ds/sx = 0 if energy Z' less mass rest
  if( XEv*E0 <= MA ){ return 0.0; }
  
  // Get sq for varible and constant
  double XEv2 = XEv*XEv, Mmu2 = Mmu*Mmu, MA2 = MA*MA, E02 = E0*E0;

  // Limits on u(x)
  double uMax   = - MA2 * (1.0 - XEv) / XEv - Mmu2*XEv
       , uMin   = -  XEv*E02*ThetaMax*ThetaMax -  MA2*(1.0 - XEv)/XEv - Mmu2*XEv
       , uMax2  = uMax*uMax, uMin2 = uMin*uMin; 
  
  // Properties nucleus
  double aa     = 111.*pow(ZNucl,-1./3.)/Mel;
  double d      = 0.164*pow(ANucl,-2./3.); 
  // Transmitted momentum, GeV^2
  double t_scr  = pow(1./aa,2.) // nuclear shielding
       , t_size = d;            // nuclear size

  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;

  // Conversion coefficient from tmin to u varible
  double gZ   = 1.0 / ( 2.0*E0*(1.0 - XEv) ), gZ2 = gZ*gZ;
  // Coefficients in amplitude
  double JZ   = 2.0 * ( (2.0 - 2.0*XEv + XEv2) / (1.0 - XEv) )
       , K    = 4.0 * (MA2 + 2.0*Mmu2) * XEv
       , LZ   = 4.0 * (MA2 + 2.0*Mmu2) * (MA2*(1.0 - XEv) + Mmu2*XEv2 );
  
  // Coefficients in photon flux
  double 
  genCoef   = t_size*t_size / ( pow( t_scr - t_size, 3.0 ) ), 
        D   =  genCoef * ( (t_scr - t_size) * t_scr / (tmax + t_scr) 
                         + (t_scr - t_size) * t_size / (tmax + t_size) 
                         - 2.0 * ( t_scr - t_size ) 
                         + ( t_scr + t_size ) 
                           * std::log( (tmax + t_size)/(tmax + t_scr) )  
                         ), 
        F   = genCoef *gZ2 * (
                               (t_scr - t_size) / (tmax + t_scr) 
                             + (t_scr - t_size) / (tmax + t_size) 
                             + 2.0 * std::log( (tmax + t_size)/(tmax + t_scr) )
                             ), 
        H   = - genCoef*(t_size + t_scr), 
        I   = - 2.0*genCoef*gZ2;
  
  // Helper functions
  double lnuMax       = std::log( (gZ2*uMax2 + t_scr) / (gZ2*uMax2 + t_size) )
       , lnuMin       = std::log( (gZ2*uMin2 + t_scr) / (gZ2*uMin2 + t_size) )
       , aTantSizeMax = std::atan( gZ*uMax / sqrt(t_size) )
       , aTantSizeMin = std::atan( gZ*uMin / sqrt(t_size) )
       , aTantScrMax  = std::atan( gZ*uMax / sqrt(t_scr) )
       , aTantScrMin  = std::atan( gZ*uMin / sqrt(t_scr) );
  
  // Result integration on u
  double
  Ing1 =   JZ*F * (uMax - uMin) + K*F  * ( log(uMax / uMin) )
         - ( JZ*D + LZ*F ) / uMax 
         - K*D / ( 2.0*uMax2 ) 
         - LZ*D / ( 3.0*pow(uMax, 3.0) )
         + ( JZ*D + LZ*F ) / uMin 
         + K*D / ( 2.0*uMin2 ) 
         + LZ*D / ( 3.0*pow(uMin, 3.0) ),

  Ing2 = LZ * H
           * ( lnuMax / ( 3.0*pow( uMax, 3.0 ) )
             - 2.0*gZ2 / ( 3.0*t_size*uMax )
             + 2.0*gZ2 / ( 3.0*t_scr*uMax )
             -  ( 2.0*pow( gZ, 3.0 ) / 3.0 )
                * ( pow(t_size, -3.0/2.0)*aTantSizeMax 
                  - pow(t_scr, -3.0/2.0)*aTantScrMax 
                  )
             - lnuMin / ( 3.0*pow( uMin, 3.0 ) )
             + 2.0*gZ2 / ( 3.0*t_size*uMin )
             - 2.0*gZ2 / ( 3.0*t_scr*uMin )
             +  ( 2.0*pow( gZ, 3.0 ) / 3.0 )
                * ( pow(t_size, -3.0/2.0)*aTantSizeMin
                  - pow(t_scr, -3.0/2.0)*aTantScrMin 
                  )
             ),
  
  Ing3 = (K * H / 2.0)
           * ( lnuMax / (uMax2)
             + gZ2
               * ( std::log( (uMax2) / (gZ2*uMax2 + t_size) ) / t_size
                 - std::log( (uMax2) / (gZ2*uMax2 + t_scr) ) / t_scr
                 )
             - lnuMin / (uMin2)
             - gZ2
               * ( std::log( (uMin2) / (gZ2*uMin2 + t_size) ) / t_size
                 - std::log( (uMin2) / (gZ2*uMin2 + t_scr) ) / t_scr
                 )
             ),
  
  Ing4 = ( JZ*H + LZ*I )
         * ( lnuMax / uMax
           + 2.0* gZ * ( std::pow(t_size, -1.0/2.0) * aTantSizeMax
                      - std::pow(t_scr, -1.0/2.0) * aTantScrMax )
           - lnuMin / uMin
           - 2.0* gZ * ( std::pow(t_size, -1.0/2.0) * aTantSizeMin
                      - std::pow(t_scr, -1.0/2.0) * aTantScrMin )
           ),
  
  Ing5 = K * I
           * ( std::log( std::abs(uMax) ) * std::log( t_size / t_scr )
             + (1.0/2.0) * ( gsl_sf_dilog(- uMax2*gZ2/t_scr )
                           - gsl_sf_dilog(- uMax2*gZ2/t_size ) )
             - std::log( std::abs(uMin) ) * std::log( t_size / t_scr )
             - (1.0/2.0) * ( gsl_sf_dilog(- uMin2*gZ2/t_scr )
                           - gsl_sf_dilog(- uMin2*gZ2/t_size ) )
             ),
  
  Ing6 = JZ * I
           * ( - ( lnuMax * uMax )
             + (2/gZ)
               * ( std::pow(t_size, 1./2) * aTantSizeMax
                 - std::pow(t_scr, 1./2) * aTantScrMax )
             + ( lnuMin * uMin )
             - (2/gZ)
               * ( std::pow(t_size, 1./2) * aTantSizeMin
                 - std::pow(t_scr, 1./2) * aTantScrMin )
             );

  // Sum integrals
  double sumIng = Ing1 + Ing2 + Ing3 + Ing4 + Ing5 + Ing6;
  // Prefactor analytical integral
  double coef = sqrt( XEv*XEv - (MA2) / (E02) )*( 1.0 - XEv ) / ( XEv );
  return sumIng*coef;
}


// Below is the IWW formula for the double differential cross-section
// from 1705.01633, see e.g. corresponding second line of Eq.(25)
//
double DarkZ::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(XEv*E0 <= MA) return 0.;
  double Uxtheta = E0*E0*2.*UThetaEv*XEv + MA*MA*(1.0-XEv)/XEv + Mmu*Mmu*XEv;
  double AA = (1. - XEv + XEv*XEv/2.) / (Uxtheta*Uxtheta);
  double BB = (1. - XEv)*(1. - XEv)*(MA*MA + 2.0*Mmu*Mmu)/(Uxtheta*Uxtheta*Uxtheta*Uxtheta);
  double CC = MA*MA - Uxtheta*XEv/(1. - XEv);
  double sigma = sqrt(XEv*XEv - MA*MA/(E0*E0)) * (AA + BB*CC);
  return sigma;
}


double DarkZ::CrossSectionDSDXDPSI(double XEv, double auxpsi, double E0)
{
  if(IApprox == 1) return CrossSectionDSDXDPSI_IWW(XEv, auxpsi, E0);
  if(IApprox == 2) return CrossSectionDSDXDPSI_WW(XEv, auxpsi, E0);
  std::cout << "DarkZ: wrong value of IApprox, exiting" << std::endl;
  exit(1);
}


double DarkZ::CrossSectionDSDXDPSI_IWW(double XEv, double auxpsi, double E0)
{
  // In IWW approach we suppose that flux factor \chi doesn't depend on y and psi
  // in that method y=Emu'/Emu is muon energy fraction and auxpsi=psi^2/2 is an auxiliar variable
  // to sample muon deflection angle psi at first step, the accepted angle will be then psiAcc = sqrt(2.0*auxpsi)
  //
  double y = 1. - XEv;
  double t2 = -(y*E0*E0*2.*auxpsi + Mmu*Mmu*(1.-y)/y+Mmu*Mmu*y) + Mmu*Mmu;
  double t = MA*MA - t2;
  // ds/dpprime
  double fac1 = (1.-y)/(t*t);
  double fac2 = 1./(2.*y) + y/2.;
  double fac3 = (MA*MA + 2.0*Mmu*Mmu)*(1.-y)*(1.-y)/(t*t*y);
  double fac4 = Mmu*Mmu*(1.-y)*(1.-y)/y + MA*MA - t;
  double part1 = fac1*(fac2+fac3*fac4);
  //beta factor
  double beta = sqrt(y*y - Mmu*Mmu/(E0*E0));
  return part1*beta;
}


double DarkZ::CrossSectionDSDXDPSI_WW(double XEv, double auxpsi, double E0)
{
  if(E0*XEv < EThresh) return 0.;
  double Xmin = MA/E0;
  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  if(XEv < Xmin || XEv > Xmax) return 0.;
  double y = 1. - XEv;
  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double ta = pow(1./aa,2.);
  double td = d;
  double t2 = -(y*E0*E0*2.*auxpsi + Mmu*Mmu*(1.-y)/y + Mmu*Mmu*y) + Mmu*Mmu;
  double t  = MA*MA - t2;
  double q = t/(2.*E0*(1.0-y));
  double tmin = q*q;
  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;
  if(tmax < tmin) return 0.;
  //double flux = log(td/(tmin + ta)) - 2.0;
  double flux = -((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax)
                  - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td)));
  // MODIFIED
  //flux = -ZNucl*ZNucl*((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax)
    //              - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td)));

  //
  if(flux < 0.) return 0.;
  // ds/dpprime
  double fac1 = (1.-y)/(t*t);
  double fac2 = 1./(2.*y)+y/2.;
  double fac3 = (MA*MA + 2.0*Mmu*Mmu)*(1.-y)*(1.-y)/(t*t*y);
  double fac4 = Mmu*Mmu*(1.-y)*(1.-y)/y + MA*MA - t;
  double part1 = fac1*(fac2 + fac3*fac4);

  //beta factor
  double beta = sqrt(y*y - Mmu*Mmu/(E0*E0));

  double rescs = flux*part1*beta;
  if(std::isnan(rescs) || rescs < 0.) {
    std::cout << "DSDXDPSI: bad result = " << rescs << " XEv = " << XEv << " E0 = " << E0 << std::endl;
    rescs = 0.;
  }
  return rescs;
}


// WW cross section for the total cross section
//
double DarkZ::CrossSectionDSDXDTheta(double XEv, double ThetaEv, double E0)
{
  if(E0*XEv < EThresh) return 0.;
  double x2=XEv*XEv;
  double theta2=ThetaEv*ThetaEv;
  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double MA2= MA*MA;
  double Mmu2= Mmu*Mmu;
  double E02= E0*E0;
  double utilde = -XEv*E02*theta2-MA2*(1.0-XEv)/XEv-Mmu2*XEv;
  double utilde2=utilde*utilde;
  double ta = 1.0/(aa*aa);
  double td = d;
  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;
  double tmin= utilde2/(4.0*E02*(1.0-XEv)*(1.0-XEv));
  // I've calculated ChiWWAnalytical by using mathematica's "Integrate[...]" function
  // and converted the resulted expression to C-like form
  double ChiWWAnalytical = -ZNucl*ZNucl*((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax)
                           - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td)));
  //double ChiWWAnalytical= ZZZ*ZZZ*(log((td+tmin)/(tmin+ta))-2.0); // <--- this approximate expression is valid for MA < 500 MeV and E0=150 GeV
  double Factor1= 2.0*(2.0-2.0*XEv+x2)/(1.0-XEv);
  double Factor2= 4.0*(MA2+2.0*Mmu2)/utilde2;
  double Factor3= utilde*XEv+MA2*(1.0-XEv)+Mmu2*x2;
  double AmplZpr2WWVEGAS = Factor1+Factor2*Factor3;
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


double DarkZ::Width()
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
