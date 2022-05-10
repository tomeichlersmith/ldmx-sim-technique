/*
 * DarkAxialsAnnihilation.hh
 *
 *  Created on: Oct 6, 2020
 *      Author: celentan
 */

#ifndef INCLUDE_DARKAXIALSANNIHILATION_HH_
#define INCLUDE_DARKAXIALSANNIHILATION_HH_

//class DarkMatter;


class DarkAxialsAnnihilation : public DarkMatter
{

  public:

    DarkAxialsAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn=1., double ANuclIn=207., double ZNuclIn=82.,
                            double DensityIn=11.35, double epsilIn=0.0001, int IDecayIn=0,double rIn=1./3,double alphaDIn=0.5);

    virtual ~DarkAxialsAnnihilation();

    virtual double TotalCrossSectionCalc(double E0);
    virtual double GetSigmaTot(double E0);
    virtual bool EmissionAllowed(double E0, double DensityMat); // E0 in GeV, density in g/cm3
    virtual double CrossSectionDSDX(double Xev, double E0);
    virtual double CrossSectionDSDXDU(double Xev, double UThetaEv, double E0);
    virtual double Width();
    virtual void SetMA(double MAIn);

  private:

    double r;
    double alphaD;
    double mChi;
};

#endif /* INCLUDE_DARKAXIALSANNIHILATION_HH_ */
