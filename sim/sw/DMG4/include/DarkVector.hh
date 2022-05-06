class DarkVector : public DarkMatter
{

  public:

    DarkVector(double MAIn, double EThreshIn, double SigmaNormIn=1., double ANuclIn=207., double ZNuclIn=82., double DensityIn=11.35,
                double epsilIn=0.0001, int IDecayIn=0);

    virtual ~DarkVector();

    virtual double TotalCrossSectionCalc(double E0);
    double TotalCrossSectionCalc_WW2(double E0);
    virtual double GetSigmaTot(double E0);
    virtual double CrossSectionDSDX(double XEv, double E0);
    virtual double CrossSectionDSDX_IWW(double XEv, double E0);
    virtual double CrossSectionDSDXDU(double XEv, double UThetaEv, double E0);
    virtual double CrossSectionDSDXDTheta(double XEv, double Theta, double E0);
    virtual double Width();

  private:

    int IApprox;
    double EMinMessenger;
    double ThetaMax;
};
