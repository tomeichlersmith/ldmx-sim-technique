#ifndef DarkPhotons_h
#define DarkPhotons_h 1

#define Mel 5.11E-4 // electron mass (GeV)
#define Mmu 0.1056 // muon mass (GeV)
#define alphaEW 1.0/137.0
#define MUp 2.79 // protonMu
#define Mpr 0.938 // proton mass
#define max_uint 4294967296.0l
#define GeVtoPb 3.894E+08

#define NPTAB 15
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <bits/stdc++.h>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

struct ParamsForChi {double AA; double ZZ; double MMA; double EE0;};
struct momentum {double E0; double Theta; double Phi;};
struct frame {TLorentzVector* fEl; TLorentzVector* cm; double E;};

class DarkPhotons
{

  public:

    DarkPhotons(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                double epsilIn, std::string fname);

    ~DarkPhotons();

    void SetSigmaNorm(double SigmaNormIn);
    void ParseLHE(std::string fname);
    void ParseROOT(std::string fname);
    void MakePlaceholders();
    double DMG4CrossSectionCalc(double E0, int IApprox);
    double TotalCrossSectionCalc_IWW(double E0);
    double TotalCrossSectionCalc_WW2(double E0);
    double TotalCrossSectionCalc_WW3(double E0);
    double TotalCrossSectionCalc_WW(double E0);
    double CrossSectionDSDXDPSI_WW(double XEv, double auxpsi, double E0);
    double CrossSectionDSDXDTheta(double XEv, double ThetaEv, double E0);
    double CrossSectionDSDX_IWW(double XEv, double E0);
    double CrossSectionDSDX_WW(double XEv, double E0);
    double TotalCrossSectionCalc(double E0, bool KFactor=false);
    double TotalMuCrossSectionCalc(double E0);
    double MaxCrossSectionCalc(double E0);
    void PrepareTable();
    double GetMA() {return MA;}
    double GetEThresh() {return EThresh;}
    double GetSigmaNorm() {return SigmaNorm;}
    double Getepsil() {return epsil;}
     // usage of normalization below:   Nsign = (Naccepted/Nsimulated)*Normalization*EOT
    double GetNormalization() {return 3.0e-15 * (Density/11.35) * (207./ANucl) *
                                      epsil * epsil / (SigmaNorm * epsilBench * epsilBench);}
    double GetsigmaTot(double E0);
    double GetsigmaMax(double E0);
    bool Emission(double E0, double DensityMat, double StepLength); // E0 in GeV, density in g/cm3, StepLength in mm
    TLorentzVector* SimulateEmission(double E0, std::string type);
    TLorentzVector* MuSimulateEmission(double E0, std::string type);
    double GetAccumulatedProbability() {return AccumulatedProbability;}
    double MParent;
    frame GetMadgraphData(double E0);

  private:
    std::map< double , std::vector < frame > > mgdata;
    std::vector < std::pair < double, int > > energies;
    double MA;
    double EThresh;
    double SigmaNorm;
    double ANucl;
    double ZNucl;
    double Density;
    double epsilBench;
    double epsil;
    int nptable;
    double ep[15];
    double sigmap[15];
    double sigmax[15];

    double AccumulatedProbability;

};

#endif
