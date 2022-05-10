#define Mel 5.109989461E-04 // electron mass in GeV
#define Mmu 0.1056583745 // muon mass in GeV
//#define alphaEW 1./137.
#define alphaEW CLHEP::fine_structure_const
#define MUp 2.79 // protonMu
#define Mpr 0.938 // proton mass
#define max_uint 4294967296.0l
#define GeVtoPb 3.894E+08

#include <stdlib.h>

struct ParamsForChi {double AA; double ZZ; double MMA; double EE0;};

struct ParamsForMuonTotCS {double AA; double ZZ; double MMA; double EE0;};


class DarkMatter
{
  friend class DarkPhotons;
  friend class DarkScalars;
  friend class DarkPseudoScalars;
  friend class DarkAxials;
  friend class DarkZ;
  friend class DarkVector;
  friend class ALP;
  friend class DarkPhotonsAnnihilation;
  friend class DarkScalarsAnnihilation;
  friend class DarkPseudoScalarsAnnihilation;
  friend class DarkAxialsAnnihilation;
  public:

    DarkMatter(double MAIn, double EThreshIn, double SigmaNormIn=1., double ANuclIn=207., double ZNuclIn=82., double DensityIn=11.35,
               double epsilIn=0.0001, int IDecayIn=0);

    virtual ~DarkMatter();

    void SetSigmaNorm(double SigmaNormIn);
    void ResetNEmissions() {NEmissions = 0;} // For G4 DM classes
    void EmissionSimulated() {NEmissions++;} // For G4 DM classes; in future do it automatically in SimulateEmission
    virtual double TotalCrossSectionCalc(double E0) = 0;
    void PrepareTable();
    double GetMA() {return MA;}
    virtual void SetMA(double MAIn) {MA = MAIn;}
    double GetEThresh() {return EThresh;}
    double GetSigmaNorm() {return SigmaNorm;}
    double Getepsil() {return epsil;}
    double GetDensity() {return Density;}
    int GetDMType() {return DMType;}
    int GetParentPDGID() {return ParentPDGID;}
    int GetDaughterPDGID() {return DaughterPDGID;}
    int Decay() {return IDecay;}
     // usage of normalization below:   Nsign = (Naccepted/Nsimulated)*Normalization*EOT
    double GetNormalization() {return 3.0e-15 * (Density/11.35) * (207./ANucl) *
                                      epsil * epsil / (SigmaNorm * epsilBench * epsilBench);}
    double GetMeanFreePathFactor() {return 1./(GetNormalization());}
    double GetSigmaTot0(double E0);
    virtual double GetSigmaTot(double E0) = 0;
    double GetSigmaMax(double E0);
    double GetSigmaAngleMax(double E0);
    double GetSigmaPsiMax(double E0);
    double GetSigmaThetaMax(double E0);
    virtual bool EmissionAllowed(double E0, double DensityMat); // E0 in GeV, density in g/cm3
    bool Emission(double E0, double DensityMat, double StepLength); // E0 in GeV, density in g/cm3, StepLength in mm
    virtual double CrossSectionDSDX(double XEv, double E0) = 0;
    virtual double CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) = 0;
    virtual double CrossSectionDSDXDPSI(double XEv, double auxpsi, double E0) {(void)XEv; (void)auxpsi; (void)E0; abort();}
    virtual double CrossSectionDSDXDPSI_IWW(double XEv, double auxpsi, double E0) {(void)XEv; (void)auxpsi; (void)E0; abort();}
    virtual double CrossSectionDSDXDPSI_WW(double XEv, double auxpsi, double E0) {(void)XEv; (void)auxpsi; (void)E0; abort();}
    virtual double CrossSectionDSDXDTheta(double XEv, double auxpsi, double E0) {(void)XEv; (void)auxpsi; (void)E0; abort();}
    virtual double CrossSectionDSDX_WW(double XEv, double E0) {(void)XEv; (void)E0; abort();}
    virtual double CrossSectionDSDX_IWW(double XEv, double E0) {(void)XEv; (void)E0; abort();}
    virtual double Width() = 0;
    double MaxCrossSectionCalc(double E0);
    double MaxCrossSectionAngleCalc(double E0);
    double MaxCrossSectionPsiCalc(double E0);
    double MaxCrossSectionThetaCalc(double E0);
    double SimulateEmission(double E0, double* angles);
    double SimulateEmissionWithAngle(double E0, double* angles);
    double SimulateEmissionWithAngle2(double E0, double* angles);
    double SimulateEmissionByMuon2(double E0, double* angles);
    double SimulateEmissionByMuon(double E0, double* angles);
    double SimulateEmissionVector(double E0, double* angles);
    double GetAccumulatedProbability() {return AccumulatedProbability;}

  private:

    double MA;
    double EThresh;
    double SigmaNorm;
    double ANucl;
    double ZNucl;
    double Density;
    double epsilBench;
    double epsil;
    int DMType; // 1 - Dark Photon; 2 - Dark Scalar; 3 - Dark Axials; 4 - Dark Pseudoscalars; 11 - Z'; 21 - ALP
    int ParentPDGID;
    double MParent;
    int DaughterPDGID;
    int IDecay; // 0 - DM particle does not decay; 1 - DM particle decays; 2 - Force DM particle decay
    int ISampler; // if use a special sampler to simulate X, angle

    int nptable;
    double ep[16];
    double sigmap[16];
    double sigmax[16];
    double sigmaxa[16];
    double sigmaxpsi[16];
    double sigmaxtheta[16];

    double AccumulatedProbability;

    int NEmissions;
};
