// This is a base class to simulate dark photons and dark scalars A production by electrons in matter eN -> eNA
// To be used in a Geant4 application.
//
// Reference: arXiV:1604.08432 [hep-ph] Phys. Rev. D 94, 095025 (2016);  arXiV:1712.05706 [hep-ph] Phys.Lett. B782 (2018) 406-411
//
// Concrete implementations: DarkPhotons, DarkScalars, ...
//
#include "DarkMatter.hh"
#include "Utils.hh"

#include "Randomize.hh"

#include "DarkMatterSampler.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

#include <iostream>

#define NPTAB 16


DarkMatter::DarkMatter(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                       double epsilIn, int IDecayIn)
:MA(MAIn), EThresh(EThreshIn), SigmaNorm(SigmaNormIn),
ANucl(ANuclIn), ZNucl(ZNuclIn), Density(DensityIn), epsilBench(0.0001), epsil(epsilIn), IDecay(IDecayIn),
AccumulatedProbability(0.), NEmissions(0)
{
  if(MA > 3.) {std::cout << "Maximal allowed mass is 3 GeV, exiting" << std::endl; exit(1);}
  nptable = NPTAB;
  double epi[NPTAB]={0.008, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10., 15., 25., 50., 80., 150., 200.};
  for(int ip=0; ip < nptable; ip++) {ep[ip] = epi[ip];}
}


DarkMatter::~DarkMatter()
{;}


void DarkMatter::PrepareTable()
{
  ISampler = 0;
  if(DMType == 1) ISampler = 1; // available only for vector; not yet fully tested
  //ISampler = 0; // Force work without sampler 
  MParent = 0.;
  if(fabs(ParentPDGID) == 11) MParent = Mel;
  if(fabs(ParentPDGID) == 13) MParent = Mmu;
  //if(fabs(ParentPDGID) == 11 && MA < 0.001) return;
  for(int ip=0; ip < nptable; ip++) {
    sigmap[ip] = TotalCrossSectionCalc(ep[ip]);
    sigmax[ip] = MaxCrossSectionCalc(ep[ip]);
    if(MA >= 0.001) sigmaxa[ip] = MaxCrossSectionAngleCalc(ep[ip]);
    if(fabs(ParentPDGID) == 13) sigmaxpsi[ip] = MaxCrossSectionPsiCalc(ep[ip]);
    if(fabs(ParentPDGID) == 13) sigmaxtheta[ip] = MaxCrossSectionThetaCalc(ep[ip]);
  }
}


double DarkMatter::GetSigmaTot0(double E0)
{
  return parinv(E0, ep, sigmap, nptable);
}


double DarkMatter::GetSigmaMax(double E0)
{
  return parinv(E0, ep, sigmax, nptable);
}


double DarkMatter::GetSigmaAngleMax(double E0)
{
  return parinv(E0, ep, sigmaxa, nptable);
}


double DarkMatter::GetSigmaPsiMax(double E0)
{
  return parinv(E0, ep, sigmaxpsi, nptable);
}


double DarkMatter::GetSigmaThetaMax(double E0)
{
  return parinv(E0, ep, sigmaxtheta, nptable);
}


bool DarkMatter::EmissionAllowed(double E0, double DensityMat)
{
  if(E0 < 1.001*MA) return false;
  if(E0 < EThresh) return false;
  if(NEmissions) return false; // For G4 DM classes
  if(fabs(DensityMat - Density) > 0.1) return false;
  return true;
}


bool DarkMatter::Emission(double E0, double DensityMat, double StepLength)
{
  if(E0 < 1.001*MA) return false;
  if(E0 < EThresh) return false;
  if(fabs(DensityMat - Density) > 0.1) return false;
  double prob = SigmaNorm*GetSigmaTot(E0)*StepLength;
  AccumulatedProbability += prob;
  double tmprandom = G4UniformRand();
  if(tmprandom < prob) return true;
  return false;
}


double DarkMatter::MaxCrossSectionCalc(double E0)
{
  if(E0 < 2.*MA) return 0.;
  if(ParentPDGID == 13 && E0 < EThresh) return 0.;

  double Xmin = MA/E0;
  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;

  double csmax = 0.;

  if(ParentPDGID == 22 || ParentPDGID == -11) Xmax = 0.99999;

  if(MA >= 0.001) csmax = CrossSectionDSDX(Xmax, E0); // preliminary
  for(int i=0; i<10000; i++) {
    double xi = 0.00005 + 0.0001*((double)i);
    if(xi >= Xmin && xi <= Xmax) {
      double csi = CrossSectionDSDX(xi, E0);
      if(MA < 0.001 && xi < EThresh/E0) csi = 0.;        // we cut DM at EThresh for these masses
      if(ParentPDGID == 13 && xi < EThresh/E0) csi = 0.; // we cut DM at EThresh for the muon beam
      if(csi > csmax) csmax = csi;
    }
  }
  std::cout << " E0 = " << E0 << "  Max cross section = " << csmax << std::endl;
  return 1.1*csmax;
}


double DarkMatter::MaxCrossSectionAngleCalc(double E0)
{
  double Xmin = MA/E0;
  double Xmax;

  double csmax = 0.;

  if(E0 < 2.*MA) return 0.;

  Xmax = 1.0-Xmin; // Incorrect limit, but for MA > 100 MeV works only like this without sampler
  if(MA <= 0.02 || ISampler) Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  if(ParentPDGID == 22) Xmax = 0.99999;

  csmax = CrossSectionDSDXDU(Xmax, 0., E0);
  double csmax1 = CrossSectionDSDXDU(Xmax, 0.000005 , E0);
  if(csmax1 > csmax) csmax = csmax1;

  for(int i=0; i<1000; i++) {
    double xi = 0.0005 + 0.001*((double)i);
    if(xi >= Xmin && xi <= Xmax) {
      double csi = CrossSectionDSDXDU(xi, 0., E0);
      if(csi > csmax) csmax = csi;
      csi = CrossSectionDSDXDU(xi, 0.000005, E0);
      if(csi > csmax) csmax = csi;
    }
  }
  return 1.1*csmax;
}


double DarkMatter::MaxCrossSectionPsiCalc(double E0)
{
  double Xmin = MA/E0;
  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  double psimax = 160.*MA/E0;
  if(psimax > 1.) psimax = 1.;

  double csmax = 0.;

  if(E0 < 2.*MA) return 0.;

  double xi, psii, auxpsi, csi;
  for(int i=0; i < 5000000; i++) {
    xi = G4UniformRand() * (Xmax-Xmin) + Xmin;
    psii = G4UniformRand() * psimax;
    auxpsi = 0.5 * psii*psii;
    csi = CrossSectionDSDXDPSI(xi, auxpsi, E0);
    if(csi > csmax) csmax = csi;
  }
  return 1.1*csmax;
}

double DarkMatter::MaxCrossSectionThetaCalc(double E0)
{
  double Xmin = MA/E0;
  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  double thetamax = 0.1;

  double csmax = 0.;

  if(E0 < 2.*MA) return 0.;

  double xi, thetai, csi;
  for(int i=0; i < 5000000; i++) {
    xi = G4UniformRand() * (Xmax-Xmin) + Xmin;
    thetai = G4UniformRand() * thetamax;
    csi = CrossSectionDSDXDTheta(xi, thetai, E0);
    if(csi > csmax) csmax = csi;
  }
  return 1.1*csmax;
}


double DarkMatter::SimulateEmission(double E0, double* angles)
{
  double Xmin = MA/E0;
  if(MA < 0.001 && EThresh/E0 > Xmin) Xmin = EThresh/E0;
  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;

  if(ParentPDGID == 22 || ParentPDGID == -11) {
    Xmin = 0.999;
    Xmax = 0.99999;
  }

  double sigmaMax = GetSigmaMax(E0);
  if(MA < 0.001 && sigmaMax > 1.2) { // diff. cross section normalized to 1.
    std::cout << "Strange too big sigma max for the mass below 0.001, exiting" << std::endl; exit(1);
  }

  int maxiter = 1000000;

  double XAcc, ThetaAcc, PhiAcc;

  for( int iii = 1; iii < maxiter; iii++) {

    double XEv  =  G4UniformRand() * (Xmax-Xmin) + Xmin;
    double UThetaEv = 0.; // we set angles to zero, generate only x

    //Now we sample only Diff. c.s. for X below 

    if(XEv*E0 < MA) return 0.;

    double sigma = CrossSectionDSDX(XEv, E0);

    double UU = G4UniformRand() * sigmaMax;

    if(sigma > sigmaMax) printf ("Maximum violated: ratio = % .18f\n", sigma/sigmaMax);

    if(sigma >= UU) {
      XAcc = XEv;
      ThetaAcc = sqrt(2.0*UThetaEv);
      //PhiAcc = G4UniformRand() * 2. * 3.1415926;
      PhiAcc = 0.;

      /*
      printf ("Accepted at iteration %d\n", iii);
      printf( "EParent = %e XAcc = %e ThetaAcc = %e\n ", E0, XAcc, ThetaAcc);
      */

      angles[0] = ThetaAcc;
      angles[1] = PhiAcc;
      return XAcc;
    }
  }
  printf ("Simulation of emission failed after N iterations = %d\n", maxiter);
  return 0.;
}


double DarkMatter::SimulateEmissionWithAngle(double E0, double* angles)
{
  double Xmin = MA/E0;
  if(MA < 0.001 && EThresh/E0 > Xmin) Xmin = EThresh/E0;

  if(ParentPDGID == 22) {
    std::cout << "ALP: Error: double differential cross section DSDXDU is not implemented, exiting" << std::endl;
    exit(1);
  }
  if(ParentPDGID == 13) {
    std::cout << "DarkZ: Error: procedure of simulation Z with angle not implemented, exiting" << std::endl;
    exit(1);
  }
  if(ParentPDGID == -11) {
    std::cout << "DarkPhotonsAnnihilation: Error: annihilation simulation with angle not implemented, exiting" << std::endl;
    exit(1);
  }

  if(MA < 0.001) {
    std::cout << "Error: mass < 0.001, don't use SimulateEmissionWithAngle" << std::endl;
    exit(1);
  }

  if(!ISampler) { // Don't use external sampler DarkMatterSampler

    double Xmax = 1. - Xmin; // Incorrect limit, but for MA > 100 MeV works only like this without sampler
    if(MA <= 0.02) Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;

    //double ThetaMaxA = 0.0001*sqrt((MA/E0)/(0.001/100.));
    double ThetaMaxA = 0.0001*sqrt(MA/0.001)*(100./E0);
    double UThetaMaxA = 0.5*ThetaMaxA*ThetaMaxA; // Nota Bene !!! this is maximum of u= 0.5*theta^2 variable!!
    if(MA < 0.001) UThetaMaxA = 0.; // Angle is simulated only for MA > 0.001 GeV
    double sigmaMax = GetSigmaAngleMax(E0);
    int maxiter = 3000000;

    double XAcc, ThetaAcc, PhiAcc;

    for( int iii = 1; iii < maxiter; iii++) {

      double XEv, FactorSigma=1.;
      if(MA >= 0.001) {
        double XFactor = 1.5*sqrt(0.001/MA);
        double AlphaX = exp(-(1. - Xmax)/XFactor);
        double BetaX = exp(-(1. - Xmin)/XFactor);
        double DeltaX = - XFactor * log(BetaX+G4UniformRand()*(AlphaX-BetaX));
        XEv = 1. - DeltaX;
        FactorSigma = exp(DeltaX/XFactor);
      } else {
        XEv  =  G4UniformRand() * (Xmax-Xmin) + Xmin;
      }

      double UThetaEv, FactorSigmaU=1.;
      if(MA >= 0.001) {
        double UFactor = 0.3*UThetaMaxA;
        double BetaU = exp(-UThetaMaxA/UFactor);
        UThetaEv = - UFactor * log(BetaU+G4UniformRand()*(1.-BetaU));
        FactorSigmaU = exp(UThetaEv/UFactor);
      } else {
        UThetaEv = G4UniformRand() * UThetaMaxA; // this is a u = 0.5*theta^2 variable!!!
      }

      if(XEv*E0 < MA) return 0.;

      double sigma = FactorSigma * FactorSigmaU * CrossSectionDSDXDU(XEv, UThetaEv, E0);

      double UU = G4UniformRand() * sigmaMax;

      if(sigma > sigmaMax) {
        printf ("Maximum violated: ratio = % .18f\n", sigma/sigmaMax);
        sigmaMax = 1.05*sigma;
      }

      if(sigma >= UU) {
        XAcc = XEv;
        ThetaAcc =sqrt(2.0*UThetaEv); // this is just a theta accepted!!!
        PhiAcc = G4UniformRand() * 2. * 3.1415926;

        /*
        printf ("Accepted at iteration %d\n", iii);
        printf( "EParent = %e XAcc = %e ThetaAcc = %e\n ", E0, XAcc, ThetaAcc);
        */

        angles[0] = ThetaAcc;
        angles[1] = PhiAcc;
        return XAcc;
      }
    }
    printf ("Simulation of emission failed after N iterations = %d\n", maxiter);
    return 0.;

  } else { // Use external sampler

    // NOTE: the engine later must be supplied by the G4-physics instance
    dphmc_URandomState state = {CLHEP::HepRandom::getTheEngine()};

    double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;

    Sampler s(E0, MA);
    double x, theta;
    double accProb = s.sample_x_theta(&state, x, theta, Xmin, Xmax);
    (void)accProb; // to avoid warning
    if(x < Xmin || x > Xmax) {
      std::cout << "SimulateEmissionWithAngle: error, X from sampler beyond limits, exiting" << std::endl;
      exit(1);
    } else {
      angles[0] = theta;
      angles[1] = G4UniformRand() * 2. * 3.1415926;
      return x;
    }
    std::cout << "SimulateEmissionWithAngle: simulation of emission with a sampler failed" << std::endl;
    return 0.;
  }
}


double DarkMatter::SimulateEmissionWithAngle2(double E0, double* angles)
{
  double Xmin = MA/E0;
  if(MA < 0.001 && EThresh/E0 > Xmin) Xmin = EThresh/E0;

  if(ParentPDGID == 22) {
    std::cout << "ALP: Error: double differential cross section DSDXDU is not implemented, exiting" << std::endl;
    exit(1);
  }
  if(ParentPDGID == 13) {
    std::cout << "DarkZ: Error: procedure of simulation Z with angle not implemented, exiting" << std::endl;
    exit(1);
  }
  if(ParentPDGID == -11) {
    std::cout << "DarkPhotonsAnnihilation: Error: annihilation simulation with angle not implemented, exiting" << std::endl;
    exit(1);
  }

  angles[0] = 0.;
  angles[1] = 0.;

  if(!ISampler || MA < 0.001) { // Don't use external sampler DarkMatterSampler

    double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
    if(Xmin > Xmax) return 0.;

    if(ParentPDGID == 22 || ParentPDGID == -11) {
      Xmin = 0.999;
      Xmax = 0.99999;
    }
    double sigmaMax = GetSigmaMax(E0);
    if(MA < 0.001 && sigmaMax > 1.2) { // diff. cross section normalized to 1.
      std::cout << "Strange too big sigma max for the mass below 0.001, exiting" << std::endl; exit(1);
    }
    int maxiterX = 1000000;

    double XAcc = 0., ThetaAcc = 0., PhiAcc = 0.;
    int NIterX = 0;

    for(int iii = 1; iii < maxiterX; iii++) { // X simulation loop

      double XEv  =  G4UniformRand() * (Xmax-Xmin) + Xmin;
      if(XEv*E0 < MA) return 0.;

      double sigma = CrossSectionDSDX(XEv, E0);
      if(sigma > sigmaMax) printf ("Maximum of single diff. CS violated: ratio = % .18f\n", sigma/sigmaMax);

      double UU = G4UniformRand() * sigmaMax;

      if(sigma >= UU) {
        XAcc = XEv;
        NIterX = iii;
        break;
      }
    }
    if(XAcc < 0.5*Xmin) {
      printf ("Simulation of X failed after N iterations = %d\n", maxiterX);
      return 0.;
    }
    if(MA < 0.001) { // No angle sampling for these masses
      std::cout << "Accepted after " << NIterX << " iterations for X " << std::endl;
      return XAcc;
    }

    double ThetaMaxA = 0.0002*pow((MA/0.001), 0.7)*(100./E0);
    if(XAcc > 0.999) ThetaMaxA *= 0.5;
    if(XAcc > 0.9999) ThetaMaxA *= 0.5;
    if(ThetaMaxA > 1.) ThetaMaxA = 1.;
    double UThetaMaxA = 0.5*ThetaMaxA*ThetaMaxA; // Nota Bene !!! this is maximum of u= 0.5*theta^2 variable!!
    double UThetaEv, sigma;

    int NIterMax = 100000.;
    sigmaMax = 0.;
    for(int iii = 0; iii < NIterMax; iii++) {
      UThetaEv = ((double)iii) * (0.1 * UThetaMaxA / (double)NIterMax);
      sigma = CrossSectionDSDXDU(XAcc, UThetaEv, E0);
      if(sigma > sigmaMax) sigmaMax = sigma;
    }
    sigmaMax *= 1.5;

    int maxiterA = 2000000;
    for(int iii = 1; iii < maxiterA; iii++) { // Angle simulation loop

      UThetaEv = UThetaMaxA * G4UniformRand();
      sigma = CrossSectionDSDXDU(XAcc, UThetaEv, E0);

      if(sigma > sigmaMax) {
        printf ("Maximum violated: ratio = % .18f\n", sigma/sigmaMax);
        sigmaMax = 1.05*sigma;
      }

      double UU = G4UniformRand() * sigmaMax;

      if(sigma >= UU) {
        ThetaAcc =sqrt(2.0*UThetaEv); // this is just a theta accepted!!!
        PhiAcc = G4UniformRand() * 2. * 3.1415926;

        /*
        std::cout << "Accepted after " << NIterX << " iterations for X and " << iii << " iterations for Angle" << std::endl;
        printf( "EParent = %e XAcc = %e ThetaAcc = %e\n ", E0, XAcc, ThetaAcc);
        */

        angles[0] = ThetaAcc;
        angles[1] = PhiAcc;
        return XAcc;
      }
    }
    std::cout << "Simulation of Angle failed after N iterations = " << maxiterA << " ,X = " << XAcc << std::endl;
    return 0.;

  } else { // Use external sampler

    // NOTE: the engine later must be supplied by the G4-physics instance
    dphmc_URandomState state = {CLHEP::HepRandom::getTheEngine()};

    double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;

    Sampler s(E0, MA);
    double x, theta;
    double accProb = s.sample_x_theta(&state, x, theta, Xmin, Xmax);
    (void)accProb; // to avoid warning
    if(x < Xmin || x > Xmax) {
      std::cout << "SimulateEmissionWithAngle: error, X from sampler beyond limits, exiting" << std::endl;
      exit(1);
    } else {
      angles[0] = theta;
      angles[1] = G4UniformRand() * 2. * 3.1415926;
      return x;
    }
    std::cout << "SimulateEmissionWithAngle: simulation of emission with a sampler failed" << std::endl;
    return 0.;
  }
}


// Z' sampling in 2 steps using DSDX and DSDXDPSI
double DarkMatter::SimulateEmissionByMuon2(double E0, double* angles)
{
  double Xmin = MA/E0;
  if(EThresh/E0 > Xmin) Xmin = EThresh/E0;

  if(abs(ParentPDGID) != 13) {
    std::cout << "Error: SimulateEmissionByMuon2: this is only for muons, exiting" << std::endl;
    exit(1);
  }

  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;

  double sigmaMax = GetSigmaMax(E0);
  int maxiterX = 1000000;

  double XAcc = 0., PsiAcc, PhiAcc;
  int NIterX = 0;

  for(int iii = 1; iii < maxiterX; iii++) { // X simulation loop

    double XEv  =  G4UniformRand() * (Xmax-Xmin) + Xmin;
    if(XEv*E0 < MA) return 0.;

    double sigma = CrossSectionDSDX(XEv, E0);
    if(sigma < -0.5) {std::cout << "Error: CrossSectionDSDX is not properly defined, exiting" << std::endl; exit(1);}
    if(sigma > sigmaMax) printf ("Maximum of single diff. CS violated: ratio = % .18f\n", sigma/sigmaMax);

    double UU = G4UniformRand() * sigmaMax;

    if(sigma >= UU) {
      XAcc = XEv;
      NIterX = iii;
      break;
    }
  }
  if(XAcc < 0.5*Xmin) {
    printf ("Simulation of X failed after N iterations = %d\n", maxiterX);
    return 0.;
  }

  //double PsiMaxA = 0.001*(MA/0.001)*(100./E0);
  double PsiMaxA = 700.*MA/E0;
  if(PsiMaxA > 1.) PsiMaxA = 1.;
  if(PsiMaxA < 0.007*100./E0) PsiMaxA = 0.007*100./E0;
  double UPsiMaxA = 0.5*PsiMaxA*PsiMaxA; // Nota Bene !!! this is maximum of u= 0.5*theta^2 variable!!
  double UPsiEv, sigma;

  int NIterMax = 20000.;
  sigmaMax = 0.;
  for(int iii = 0; iii < NIterMax; iii++) {
    UPsiEv = ((double)iii) * (0.1 * UPsiMaxA / (double)NIterMax);
    sigma = CrossSectionDSDXDPSI(XAcc, UPsiEv, E0);
    if(sigma > sigmaMax) sigmaMax = sigma;
  }
  sigmaMax *= 1.5;

  int maxiterA = 20000000;
  for(int iii = 1; iii < maxiterA; iii++) { // Angle simulation loop

    UPsiEv = UPsiMaxA * G4UniformRand();
    sigma = CrossSectionDSDXDPSI(XAcc, UPsiEv, E0);

    if(sigma > sigmaMax) {
      printf ("Maximum violated: ratio = % .18f\n", sigma/sigmaMax);
      sigmaMax = 1.05*sigma;
    }

    double UU = G4UniformRand() * sigmaMax;

    if(sigma >= UU) {
      PsiAcc =sqrt(2.0*UPsiEv);
      PhiAcc = G4UniformRand() * 2. * 3.1415926;

      /*
      std::cout << "Accepted after " << NIterX << " iterations for X and " << iii << " iterations for Angle" << std::endl;
      printf( "EParent = %e XAcc = %e PsiAcc = %e\n ", E0, XAcc, PsiAcc);
      */

      angles[0] = PsiAcc;
      angles[1] = PhiAcc;
      return XAcc;

    }
  }
  std::cout << "Simulation of Angle failed after N iterations = " << maxiterA << " ,X = " << XAcc << std::endl;
  return 0.;
}


// Z' sampling using DSDXDPSI
double DarkMatter::SimulateEmissionByMuon(double E0, double* angles)
{
  if(ParentPDGID != 13) {
    std::cout << "Error: SimulateEmissionByMuon: this is only for muons, exiting" << std::endl;
    exit(1);
  }

  double Xmin = MA/E0;
  if(EThresh/E0 > Xmin) Xmin = EThresh/E0;
  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  double PsiMax = 700.*MA/E0;
  if(PsiMax > 1.) PsiMax = 1.;
  if(PsiMax < 0.007*100./E0) PsiMax = 0.007*100./E0;
  double UPsiMax = 0.5*PsiMax*PsiMax;

  // Parameters of adaptive sampling
  // X
  double XFactor = 1.;
  if(MA < 0.001) XFactor = 0.3;
  if(MA >= 0.001 && MA < 0.004) XFactor = 0.6;
  double AlphaX = exp(-Xmin/XFactor);
  double BetaX = exp(-Xmax/XFactor);
  double FactorSigmaX;
  // Angle
  double ARatio = 0.2; // becomes straight if = 1.
  double Ratio = 20.;
  if(MA <= 0.001) Ratio = 30.;
  double Frac = (ARatio*Ratio)/(ARatio*Ratio + (1.-ARatio));
  double FactorSigmaA;

  double sigmaMax = GetSigmaPsiMax(E0);

  int maxiter = 25000000;
  double XEv, UPsiEv, sigma, XAcc, PsiAcc, PhiAcc;

  for( int iii = 1; iii < maxiter; iii++) {

    FactorSigmaX = 1.;
    if(MA <= 0.01) { // adaptive sampling of X
      XEv = - XFactor * log(BetaX+G4UniformRand()*(AlphaX-BetaX));
      FactorSigmaX = exp(XEv/XFactor);
    } else { // straight sampling of X
      XEv  =  G4UniformRand() * (Xmax-Xmin) + Xmin;
    }
//    std::cout << "X, FactorSigmaX = " << XEv << " " << FactorSigmaX << std::endl;
    FactorSigmaA = 1.;
    if(MA <= 0.008) { // adaptive sampling of angle
      if(G4UniformRand() < Frac) {
        UPsiEv = UPsiMax * ARatio * G4UniformRand();
      } else {
        UPsiEv = UPsiMax * (ARatio + (1.-ARatio)*G4UniformRand());
        FactorSigmaA = Ratio;
      }
    } else { // straight sampling of angle
      UPsiEv = UPsiMax * G4UniformRand();
    }

    sigma = FactorSigmaX * FactorSigmaA * CrossSectionDSDXDPSI(XEv, UPsiEv, E0);

    if(sigma > sigmaMax) {
      printf ("Maximum violated: ratio = % .18f\n", sigma/sigmaMax);
      sigmaMax = 1.05*sigma; // actually will not not work, the values in the table to be scaled 
    }

    double UU = G4UniformRand() * sigmaMax;

    if(sigma >= UU) {
      XAcc = XEv;
      PsiAcc = sqrt(2.0*UPsiEv);
      PhiAcc = G4UniformRand() * 2. * 3.1415926;

      /*
      std::cout << "Accepted after " << iii << " iterations" << std::endl;
      printf( "EParent = %e XAcc = %e PsiAcc = %e\n ", E0, XAcc, PsiAcc);
      */

      angles[0] = PsiAcc;
      angles[1] = PhiAcc;
      return XAcc;
    }
  }
  printf ("Simulation of emission failed after N iterations = %d\n", maxiter);
  return 0.;
}


// Vector sampling using DSDXDTheta
//
double DarkMatter::SimulateEmissionVector(double E0, double* angles)
{
  double Xmin = MA/E0;
  double Xmax = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MParent/E0;
  double ThetaMax = 0.1;

  double sigmaMax = GetSigmaThetaMax(E0);

  int maxiter = 20000000;
  double XEv, ThetaEv, sigma, XAcc, ThetaAcc, PhiAcc;

  for( int iii = 1; iii < maxiter; iii++) {

    XEv  =  G4UniformRand() * (Xmax-Xmin) + Xmin;
    ThetaEv = ThetaMax * G4UniformRand();

    sigma = CrossSectionDSDXDTheta(XEv, ThetaEv, E0);

    if(sigma > sigmaMax) {
      printf ("Maximum violated: ratio = % .18f\n", sigma/sigmaMax);
      sigmaMax = 1.05*sigma;
    }

    double UU = G4UniformRand() * sigmaMax;

    if(sigma >= UU) {
      XAcc = XEv;
      ThetaAcc = ThetaEv;
      PhiAcc = G4UniformRand() * 2. * 3.1415926;

      /*
      std::cout << "Accepted after " << iii << " iterations" << std::endl;
      printf( "EParent = %e XAcc = %e ThetaAcc = %e\n ", E0, XAcc, ThetaAcc);
      */

      angles[0] = ThetaAcc;
      angles[1] = PhiAcc;
      return XAcc;
    }
  }
  printf ("Simulation of emission failed after N iterations = %d\n", maxiter);
  return 0.;
}
