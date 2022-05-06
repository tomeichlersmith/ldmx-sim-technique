class DarkPhotons : public DarkMatter
{

  public:

    DarkPhotons(double MAIn, double EThreshIn, double SigmaNormIn=1., double ANuclIn=207., double ZNuclIn=82., double DensityIn=11.35,
                double epsilIn=0.0001, int IDecayIn=0);

    virtual ~DarkPhotons();

    virtual double TotalCrossSectionCalc(double E0);
    virtual double GetSigmaTot(double E0);
    virtual double CrossSectionDSDX(double Xev, double E0);
    virtual double CrossSectionDSDXDU(double Xev, double UThetaEv, double E0);
    virtual double Width();

  private:

};
