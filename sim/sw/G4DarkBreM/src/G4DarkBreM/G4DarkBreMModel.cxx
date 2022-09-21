
#include "G4DarkBreM/G4DarkBreMModel.h"

#include "Framework/Exception/Exception.h"
#include "Framework/Logger.h"
#include "G4DarkBreM/G4APrime.h"

// Geant4
#include "G4Electron.hh"
#include "G4MuonMinus.hh"
#include "G4EventManager.hh"  //for EventID number
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"  //for VerboseLevel
#include "G4SystemOfUnits.hh"

// Boost
#include <boost/math/quadrature/gauss_kronrod.hpp>

// STL
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace simcore {
namespace darkbrem {

/**
 * The integration method we will be using for our numerical integrals
 *
 * The Gauss-Kronrod method was chosen due to its ability to limit the
 * number of calls to the function representing the integrand which
 * should help improve performance for us due to the complexity of our
 * integrand. The order of the GK method was chosen after some 
 * experimentation, starting at a high value (61) and then lowering
 * it to achieve better performance while checking the accuracy of
 * the results.
 *
 * As explained in the [Boost GK Docs](https://www.boost.org/doc/libs/master/libs/math/doc/html/math_toolkit/gauss_kronrod.html),
 * generally the error estimation technique for this method is
 * overly pessimistic, so we can confidently set the maximum
 * depth low and the desired relative error high compared
 * to other methods. We have followed the examples in the docs
 * where we use max_depth to 5 and relative error to 1e-9.
 */
using int_method = boost::math::quadrature::gauss_kronrod<double, 61>;

/**
 * numerically integrate the value of the flux factory chi
 *
 * The integration of the form factor into the flux factor can
 * be done analytically with a tool like mathematica, but when
 * including the inelastic term, it produces such a complicated 
 * result that the numerical integration is actually *faster*
 * than the analytical one.
 *
 * The form factors are copied from Appendix A (Eq A18 and A19) of
 * https://journals.aps.org/prd/pdf/10.1103/PhysRevD.80.075018
 */
static double flux_factor_chi_numerical(G4double A, G4double Z, double tmin, double tmax) {
  /**
   * bin = (mu_p^2 - 1)/(4 m_pr^2)
   * mel = mass of electron in GeV
   */
  static const double bin = (2.79*2.79 - 1)/(4*0.938*0.938),
                      mel = 0.000511;
  const double ael = 111.0*pow(Z,-1./3.)/mel,
               del = 0.164*pow(A,-2./3.),
               ain = 773.0*pow(Z,-2./3.)/mel,
               din = 0.71;

  auto integrand = [&](double t) {
    double ael_factor = (ael*ael*t)/(1 + ael*ael*t),
           del_factor = 1./(1+t/del),
           ain_factor = (ain*ain*t)/(1 + ain*ain*t),
           din_factor = 1./(1+t/din),
           nucl = (1 + t*bin);
    
    return (pow(ael_factor*del_factor*Z, 2)
            +
            Z*pow(ain_factor*nucl*din*din*din*din, 2)
           )*(t-tmin)/pow(t,2.);
  };

  return int_method::integrate(integrand,tmin,tmax,5,1e-9);
}

/**
 * analytic flux factor chi integrated and simplified by DMG4 authors
 *
 * This only includes the elastic form factor term
 */
static double flux_factor_chi_analytic(G4double A, G4double Z, double tmin, double tmax) {
  static const double mel = 0.000511;
  const double a_el = 111.*pow(Z,-1./3)/mel,
               d_el = 0.164*pow(A,-2./3);
  double ta = 1.0/(a_el*a_el);
  double td = d_el;
  return -Z*Z*((td*td*(
              ((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) 
              + (ta + td + 2.0*tmin)*(log(ta + tmax) - log(td + tmax) - log(ta + tmin) + log(td + tmin))
             ))/((ta-td)*(ta-td)*(ta-td)));
}

G4DarkBreMModel::G4DarkBreMModel(framework::config::Parameters &params, bool muons)
    : G4DarkBremsstrahlungModel(params, muons), method_(DarkBremMethod::Undefined) {
  method_name_ = params.getParameter<std::string>("method");
  if (method_name_ == "forward_only") {
    method_ = DarkBremMethod::ForwardOnly;
  } else if (method_name_ == "cm_scaling") {
    method_ = DarkBremMethod::CMScaling;
  } else if (method_name_ == "undefined") {
    method_ = DarkBremMethod::Undefined;
  } else {
    EXCEPTION_RAISE("InvalidMethod", "Invalid dark brem simulation method '" +
                                         method_name_ + "'.");
  }

  threshold_ = std::max(
      params.getParameter<double>("threshold"),
      2. * G4APrime::APrime()->GetPDGMass() / CLHEP::GeV  // mass A' in GeV
  );

  epsilon_ = params.getParameter<double>("epsilon");

  library_path_ = params.getParameter<std::string>("library_path");
  if (params.getParameter("load_library",true))
    SetMadGraphDataLibrary(library_path_);
}

void G4DarkBreMModel::PrintInfo() const {
  G4cout << " Dark Brem Vertex Library Model" << G4endl;
  G4cout << "   Threshold [GeV]: " << threshold_ << G4endl;
  G4cout << "   Epsilon:         " << epsilon_ << G4endl;
  G4cout << "   Scaling Method:  " << method_name_ << G4endl;
  G4cout << "   Vertex Library:  " << library_path_ << G4endl;
}

void G4DarkBreMModel::RecordConfig(ldmx::RunHeader &h) const {
  h.setFloatParameter("Minimum Threshold to DB [GeV]", threshold_);
  h.setFloatParameter("DB Xsec Epsilon", epsilon_);
  h.setStringParameter("Vertex Scaling Method", method_name_);
  h.setStringParameter("Vertex Library", library_path_);
}

G4double G4DarkBreMModel::ComputeCrossSectionPerAtom(
    G4double lepton_ke, G4double A, G4double Z) {
  static const double MA = G4APrime::APrime()->GetPDGMass() / GeV;
  static const double MA2 = MA*MA;

  const double lepton_mass{
    (muons_ ? G4MuonMinus::MuonMinus()->GetPDGMass() : G4Electron::Electron()->GetPDGMass()) / GeV};
  const double lepton_mass_sq{lepton_mass*lepton_mass};

  // the cross section is zero if the lepton does not have enough
  // energy to create an A'
  // the threshold_ can also be set by the user to a higher value
  // to prevent dark-brem within inaccessible regions of phase
  // space
  if (lepton_ke < keV or lepton_ke < threshold_*GeV) return 0.;

  // Change energy to GeV.
  double lepton_e = lepton_ke/GeV + lepton_mass;
  double lepton_e_sq = lepton_e*lepton_e;

  /**
   * IWW
   *
   * assume theta = 0, and x = 1 for form factor integration
   * i.e. now chi is a constant pulled out of the integration
   */
  double chi = flux_factor_chi_numerical(A,Z,MA2*MA2/(4*lepton_e_sq),MA2+lepton_mass_sq);

  /**
   * Differential cross section with respect to x and theta
   *
   * Equation (16) from Appendix A of https://arxiv.org/pdf/2101.12192.pdf
   */
  auto diff_cross = [&](double x, double theta) {
    static const G4double alphaEW = 1.0 / 137.0;
    if (x*lepton_e < threshold_) return 0.;

    double theta_sq = theta*theta;
    double x_sq = x*x;

    double utilde = -x*lepton_e_sq*theta_sq - MA2*(1.-x)/x - lepton_mass_sq*x;
    double utilde_sq = utilde*utilde;

    if (muons_) {
      /**
       * WW
       *
       * Since muons are so much more massive than electrons, the
       * keep form factor integration limits dependent on x and theta
       */
      // non-zero theta and non-zero m_l
      double tmin = utilde_sq/(4.0*lepton_e_sq*(1.0-x)*(1.0-x));
      // zero theta
      //double tmin = pow((MA2*(1-x)/x-lepton_mass_sq*x),2)/(4*lepton_e_sq*(1-x)*(1-x));
      // zero theta and zero lepton mass
      //double tmin = MA2*MA2/(4*lepton_e_sq*x_sq);
  
      // maximum t is the incident energy ?
      double tmax = lepton_e_sq;
  
      // require 0 < tmin < tmax to procede
      if (tmin < 0) return 0.;
      if (tmax < tmin) return 0.;
  
      /**
       * numerically integrate to calculate chi ourselves
       * this _has not_ been well defined probably due to the extreme values of t
       * that must be handled
      double chi = flux_factor_chi_numerical(A,Z, tmin, tmax);
       */
  
      /**
       * use analytic elastic-only chi derived for DMG4
       */
      chi = flux_factor_chi_analytic(A,Z,tmin,tmax);
    }

    /**
     * Amplitude squared is taken from 
     * Equation (17) from Appendix A of https://arxiv.org/pdf/2101.12192.pdf
     * with X = V
     */
    double factor1 = 2.0*(2.0 - 2.*x + x_sq)/(1. - x);
    double factor2 = 4.0*(MA2 + 2.0*lepton_mass_sq)/utilde_sq;
    double factor3 = utilde*x + MA2*(1. - x) + lepton_mass_sq*x_sq;
    double amplitude_sq = factor1 + factor2*factor3;

    return 2.*pow(epsilon_,2.)*pow(alphaEW,3.)
             *sqrt(x_sq*lepton_e_sq - MA2)*lepton_e*(1.-x)
             *(chi/utilde_sq)*amplitude_sq*sin(theta);
  };

  // deduce integral bounds
  double xmin = 0;
  double xmax = 1;
  if ((lepton_mass / lepton_e) > (MA / lepton_e))
    xmax = 1 - lepton_mass / lepton_e;
  else
    xmax = 1 - MA / lepton_e;

  /**
   * max recoil angle of A'
   *
   * The wide angle A' are produced at a negligible rate
   * so we enforce a hard-coded cut-off to stay within
   * the small-angle regime.
   *
   * We choose the same cutoff as DMG4.
   */
  double theta_max{0.3};

  /**
   * Numerical 2D integration
   *
   * Its pretty simple, we cast the 2D integration down to 2 1D integrations
   * by defining the theta_integral function to numerically integrate
   * the integrand at a fixed x and then we put that through the same integration
   * procedure going through the possible x.
   */
  auto theta_integral = [&](double x) {
    auto theta_integrand = [&](double theta) {
      return diff_cross(x, theta);
    };
    // integrand, min, max, max_depth, tolerance, error, pL1
    return int_method::integrate(theta_integrand, 0., theta_max, 5, 1e-9);
  };

  double error;
  double integrated_xsec = int_method::integrate(theta_integral, xmin, xmax, 5, 1e-9, &error);

  G4double GeVtoPb = 3.894E08;

  /**
   * The integrated_xsec should be the correct value, we are just
   * converting it to Geant4's pb units here
   */
  G4double cross = integrated_xsec * GeVtoPb * picobarn;

  if (cross < 0.) return 0.;  // safety check all the math

  return cross;
}

void G4DarkBreMModel::GenerateChange(
    G4ParticleChange &particleChange, const G4Track &track,
    const G4Step &step) {
  // mass A' in GeV
  static const double MA = G4APrime::APrime()->GetPDGMass() / CLHEP::GeV;

  // mass of incident lepton (usually electrons, muons experimental)
  double Mel = track.GetDefinition()->GetPDGMass() / CLHEP::GeV;

  // convert to energy units in LHE files [GeV]
  G4double incidentEnergy = step.GetPostStepPoint()->GetTotalEnergy()/CLHEP::GeV;

  OutgoingKinematics data = GetMadgraphData(incidentEnergy);
  double EAcc = (data.electron.E() - Mel) *
                    ((incidentEnergy - Mel - MA) / (data.E - Mel - MA)) +
                Mel;
  double Pt = data.electron.Pt();
  double P = sqrt(EAcc * EAcc - Mel * Mel);
  double PhiAcc = data.electron.Phi();
  if (method_ == DarkBremMethod::ForwardOnly) {
    unsigned int i = 0;
    while (Pt * Pt + Mel * Mel > EAcc * EAcc) {
      // Skip events until the transverse energy is less than the total energy.
      i++;
      data = GetMadgraphData(incidentEnergy);
      EAcc = (data.electron.E() - Mel) *
                 ((incidentEnergy - Mel - MA) / (data.E - Mel - MA)) +
             Mel;
      Pt = data.electron.Pt();
      P = sqrt(EAcc * EAcc - Mel * Mel);
      PhiAcc = data.electron.Phi();

      if (i > maxIterations_) {
        ldmx_log(warn)
            << "Could not produce a realistic vertex with library energy "
            << data.electron.E() << " GeV.\n"
            << "Consider expanding your libary of A' vertices to include a "
               "beam energy closer to "
            << incidentEnergy << " GeV.";
        break;
      }
    }
  } else if (method_ == DarkBremMethod::CMScaling) {
    TLorentzVector el(data.electron.X(), data.electron.Y(), data.electron.Z(),
                      data.electron.E());
    double ediff = data.E - incidentEnergy;
    TLorentzVector newcm(data.centerMomentum.X(), data.centerMomentum.Y(),
                         data.centerMomentum.Z() - ediff,
                         data.centerMomentum.E() - ediff);
    el.Boost(-1. * data.centerMomentum.BoostVector());
    el.Boost(newcm.BoostVector());
    double newE = (data.electron.E() - Mel) *
                      ((incidentEnergy - Mel - MA) / (data.E - Mel - MA)) +
                  Mel;
    el.SetE(newE);
    EAcc = el.E();
    Pt = el.Pt();
    P = el.P();
  } else if (method_ == DarkBremMethod::Undefined) {
    EAcc = data.electron.E();
    P = sqrt(EAcc * EAcc - Mel * Mel);
    Pt = data.electron.Pt();
  }

  // What we need:
  //  - EAcc
  //  - P and Pt for ThetaAcc
  //  - PhiAcc
  // Basically we need the 3-momentum of the recoil electron
  
  // Change the energy back to MeV, the internal GEANT unit.
  EAcc = EAcc * CLHEP::GeV;  
  Mel  = Mel  * CLHEP::GeV;

  // outgoing lepton momentum
  G4double recoilElectronMomentumMag = sqrt(EAcc * EAcc - Mel*Mel);
  G4ThreeVector recoilElectronMomentum;
  double ThetaAcc = std::asin(Pt / P);
  recoilElectronMomentum.set(std::sin(ThetaAcc) * std::cos(PhiAcc),
                             std::sin(ThetaAcc) * std::sin(PhiAcc),
                             std::cos(ThetaAcc));
  recoilElectronMomentum.rotateUz(track.GetMomentumDirection());
  recoilElectronMomentum.setMag(recoilElectronMomentumMag);

  // create g4dynamicparticle object for the dark photon.
  // define its 3-momentum so we conserve 3-momentum with primary and recoil
  // electron NOTE: does _not_ take nucleus recoil into account
  G4ThreeVector darkPhotonMomentum =
      track.GetMomentum() - recoilElectronMomentum;
  G4DynamicParticle *dphoton =
      new G4DynamicParticle(G4APrime::APrime(), darkPhotonMomentum);
  // energy of primary
  G4double finalKE = EAcc - Mel;

  // stop tracking and create new secondary instead of primary
  if (alwaysCreateNewElectron_) {
    // TODO copy over all other particle information from track I am killing
    G4DynamicParticle *el = new G4DynamicParticle(
        track.GetDefinition(), recoilElectronMomentum);
    particleChange.SetNumberOfSecondaries(2);
    particleChange.AddSecondary(dphoton);
    particleChange.AddSecondary(el);
    particleChange.ProposeTrackStatus(fStopAndKill);
    // continue tracking
  } else {
    // just have primary lose energy (don't rename to different track ID)
    // TODO untested this branch, not sure if it works as expected
    particleChange.SetNumberOfSecondaries(1);
    particleChange.AddSecondary(dphoton);
    particleChange.ProposeMomentumDirection(recoilElectronMomentum.unit());
    particleChange.ProposeEnergy(finalKE);
  }
}

void G4DarkBreMModel::SetMadGraphDataLibrary(std::string path) {
  // Assumptions:
  //  - Directory passed is a flat directory (no sub directories) containing LHE
  //  files
  //  - LHE files are events generated with the correct mass point
  // TODO automatically select LHE files of the correct mass point?

  bool foundOneFile = false;
  DIR *dir;            // handle to opened directory
  struct dirent *ent;  // handle to entry inside directory
  if ((dir = opendir(path.c_str())) != NULL) {
    // directory can be opened
    while ((ent = readdir(dir)) != NULL) {
      std::string fp = path + '/' + std::string(ent->d_name);
      if (fp.substr(fp.find_last_of('.') + 1) == "lhe") {
        // file ends in '.lhe'
        ParseLHE(fp);
        foundOneFile = true;
      }
    }
    closedir(dir);
  }

  if (not foundOneFile) {
    EXCEPTION_RAISE("DirDNE", "Directory '" + path +
                                  "' was unable to be opened or no '.lhe' "
                                  "files were found inside of it.");
  }

  MakePlaceholders();  // Setup the placeholder offsets for getting data.

  ldmx_log(info) << "MadGraph Library of Dark Brem Vertices:\n";
  for (const auto &kV : madGraphData_) {
    ldmx_log(info) << "\t" << std::setw(8) << kV.first << " GeV Beam -> "
                   << std::setw(6) << kV.second.size() << " Events";
  }

  return;
}

void G4DarkBreMModel::ParseLHE(std::string fname) {
  static const double MA =
      G4APrime::APrime()->GetPDGMass() / CLHEP::GeV;  // mass A' in GeV

  // TODO: use already written LHE parser?
  ldmx_log(info) << "Parsing LHE file '" << fname << "'... ";

  std::ifstream ifile;
  ifile.open(fname.c_str());
  if (!ifile) {
    EXCEPTION_RAISE("LHEFile", "Unable to open LHE file '" + fname + "'.");
  }

  std::string line;
  while (std::getline(ifile, line)) {
    std::istringstream iss(line);
    int ptype, state;
    double skip, px, py, pz, E, M;
    if (iss >> ptype >> state >> skip >> skip >> skip >> skip >> px >> py >>
        pz >> E >> M) {
      if ((ptype == 11 or ptype == 13) && (state == -1)) {
        double ebeam = E;
        double e_px, e_py, e_pz, a_px, a_py, a_pz, e_E, a_E, e_M, a_M;
        for (int i = 0; i < 2; i++) {
          std::getline(ifile, line);
        }
        std::istringstream jss(line);
        jss >> ptype >> state >> skip >> skip >> skip >> skip >> e_px >> e_py >>
            e_pz >> e_E >> e_M;
        if ((ptype == 11 or ptype == 13) && (state == 1)) {  // Find a final state electron.
          for (int i = 0; i < 2; i++) {
            std::getline(ifile, line);
          }
          std::istringstream kss(line);
          kss >> ptype >> state >> skip >> skip >> skip >> skip >> a_px >>
              a_py >> a_pz >> a_E >> a_M;
          if (ptype == 622 and state == 1) {
            if (abs(1. - a_M / MA) > 1e-3) {
              EXCEPTION_RAISE("BadMGEvnt",
                              "A MadGraph imported event has a different "
                              "APrime mass than the model has (MadGraph = " +
                                  std::to_string(a_M) + "GeV; Model = " +
                                  std::to_string(MA) + "GeV).");
            }
            OutgoingKinematics evnt;
            double cmpx = a_px + e_px;
            double cmpy = a_py + e_py;
            double cmpz = a_pz + e_pz;
            double cmE = a_E + e_E;
            evnt.electron = TLorentzVector(e_px, e_py, e_pz, e_E);
            evnt.centerMomentum = TLorentzVector(cmpx, cmpy, cmpz, cmE);
            evnt.E = ebeam;
            madGraphData_[ebeam].push_back(evnt);
          }  // get a prime kinematics
        }    // check for final state
      }      // check for particle type and state
    }        // able to get momentum/energy numbers
  }          // while getting lines
  // Add the energy to the list, with a random offset between 0 and the total
  // number of entries.
  ifile.close();
  ldmx_log(info) << "done parsing.";
}

void G4DarkBreMModel::MakePlaceholders() {
  currentDataPoints_.clear();
  maxIterations_ = 10000;
  for (const auto &iter : madGraphData_) {
    currentDataPoints_[iter.first] = int(G4UniformRand() * iter.second.size());
    if (iter.second.size() < maxIterations_)
      maxIterations_ = iter.second.size();
  }
}

G4DarkBreMModel::OutgoingKinematics
G4DarkBreMModel::GetMadgraphData(double E0) {
  OutgoingKinematics cmdata;  // data frame to return

  // Cycle through imported beam energies until the closest one above is found,
  // or the max is reached.
  double samplingE = 0.;
  for (const auto &keyVal : currentDataPoints_) {
    samplingE = keyVal.first;  // move samplingE up
    // check if went under the sampling energy
    //  the map is sorted by key, so we can be done right after E0 goes under
    //  samplingE
    if (E0 < samplingE) break;
  }
  // now samplingE is the closest energy above E0 or the maximum energy imported
  // from mad graph

  // Need to loop around if we hit the end, when the size of
  // madGraphData_[samplingE] is smaller than
  //  the number of events we want
  if (currentDataPoints_.at(samplingE) >= madGraphData_.at(samplingE).size()) {
    currentDataPoints_[samplingE] = 0;
  }

  // Get the lorentz vectors from the index given by the placeholder.
  cmdata = madGraphData_.at(samplingE).at(currentDataPoints_.at(samplingE));

  // Increment the current index
  currentDataPoints_[samplingE]++;

  return cmdata;
}

}  // namespace darkbrem
}  // namespace simcore
