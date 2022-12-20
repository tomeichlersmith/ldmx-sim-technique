/**
 * g++ -o xsec-dbg xsec_calc.cxx
 */

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Gauss-Kronrod
 *  - has been working for a range of incident lepton energies
 */
#include <boost/math/quadrature/gauss_kronrod.hpp>
template<typename Integrand>
double integrate(Integrand f, double low, double high) {
  return boost::math::quadrature::gauss_kronrod<double, 61>::integrate(
      f, low, high, 5, 1e-9);
}

/*
 * tanh-sinh
 *  - promoted as a method to handle pesky integrals with divergences at
 *    or near limits
 *  - produces same results as GK for higher energies,
 *    takes a /really long/ time to converge for lower energies
 *    (probably due to failure to match criteria of integration method)
#include <boost/math/quadrature/tanh_sinh.hpp>
template<typename Integrand>
double integrate(Integrand f, double low, double high) {
  static boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, low, high);
}
*/

/**
 * print out how to use g4db-xsec-calc
 */
void usage() {
  std::cout <<
    "USAGE:\n"
    "  xsec-dbg [options]\n"
    "\n"
    "Calculate muon dark brem cross sections and write them out to a CSV table\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help    : produce this help and exit\n"
    "  -t,--total   : output file to write total xsec to\n"
    "  --2d         : output file to write ds/dxdtheta to\n"
    "                 WARNING: this file gets very big if doing multiple energy samples\n"
    "  --1d         : output file to write ds/dx to\n"
    "                 WARNING: this file gets very big if doing multiple energy samples\n"
    "  -M,--ap-mass : mass of dark photon in GeV\n"
    "                 default is 1.0GeV\n"
    "  --energy     : single energy or python-like arange for input energies in GeV\n"
    "                 (single_energy OR start stop OR start stop step)\n"
    "                 default start is 2.0 and default step is 1.0 GeV\n"
    "  --target     : define target material with three parameters Z, A\n"
    "                 default is copper 29.0 63.55\n"
    << std::flush;
}


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
static double flux_factor_chi_numerical(double A, double Z, double tmin, double tmax) {
  /*
   * bin = (mu_p^2 - 1)/(4 m_pr^2)
   * mel = mass of electron in GeV
   */
  static const double bin = (2.79*2.79 - 1)/(4*0.938*0.938),
                      mel = 0.000511;
  const double ael = 111.0*pow(Z,-1./3.)/mel,
               del = 0.164*pow(A,-2./3.),
               ain = 773.0*pow(Z,-2./3.)/mel,
               din = 0.71,
               ael_inv2 = pow(ael, -2),
               ain_inv2 = pow(ain, -2);

  /**
   * We've manually expanded the integrand to cancel out the 1/t^2 factor
   * from the differential, this helps the numerical integration converge
   * because we aren't teetering on the edge of division by zero
   *
   * The `auto` used in the integrand definition represents a _function_ 
   * whose return value is a `double` and which has a single input `t`. 
   * This lambda expression saves us the time of having to re-calculate 
   * the form factor constants that do not depend on `t` because it 
   * can inherit their values from the environment. 
   * The return value is a double since it is calculated
   * by simple arithmetic operations on doubles.
   */
  auto integrand = [&](double t) {
    double ael_factor = 1./(ael_inv2 + t),
           del_factor = 1./(1+t/del),
           ain_factor = 1./(ain_inv2 + t),
           din_factor = 1./(1+t/din),
           nucl = (1 + t*bin);
    
    return (pow(ael_factor*del_factor*Z, 2)
            +
            Z*pow(ain_factor*nucl*din_factor*din_factor*din_factor*din_factor, 2)
           )*(t-tmin);
  };

  return integrate(integrand,tmin,tmax);
}

/**
 * analytic flux factor chi integrated and simplified by DMG4 authors
 *
 * This only includes the elastic form factor term
 */
static double flux_factor_chi_analytic(double A, double Z, double tmin, double tmax) {
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

/**
 * definition of g4db-xsec-calc
 */
int main(int argc, char* argv[]) try {
  double ap_mass{1.0};
  double min_energy{2.};
  double max_energy{1000.};
  double energy_step{1.0};
  double target_Z{29.};
  double target_A{63.55};
  std::string total_xsec{"total_xsec.csv"};
  std::string dsdxdtheta_fn{"/dev/null"};
  std::string dsdx_fn{"/dev/null"};
  for (int i_arg{1}; i_arg < argc; ++i_arg) {
    std::string arg{argv[i_arg]};
    if (arg == "-h" or arg == "--help") {
      usage();
      return 0;
    } else if (arg == "-t" or arg == "--total") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      total_xsec = argv[++i_arg];
    } else if (arg == "--1d") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      dsdx_fn = argv[++i_arg];
    } else if (arg == "--2d") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      dsdxdtheta_fn = argv[++i_arg];
    } else if (arg == "-M" or arg == "--ap-mass") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      ap_mass = std::stod(argv[++i_arg]);
    } else if (arg == "--energy") {
      std::vector<std::string> args;
      while (i_arg+1 < argc and argv[i_arg+1][0] != '-') {
        args.push_back(argv[++i_arg]);
      }
      if (args.size() == 0) {
        std::cerr << arg << " requires arguments after it" << std::endl;
        return 1;
      } else if (args.size() == 1) {
        min_energy = std::stod(args[0]);
        max_energy = min_energy;
      } else if (args.size() == 2) {
        min_energy = std::stod(args[0]);
        max_energy = std::stod(args[1]);
      } else if (args.size() == 3) {
        min_energy = std::stod(args[0]);
        max_energy = std::stod(args[1]);
        energy_step = std::stod(args[2]);
      }
    } else if (arg == "--target") {
      std::vector<std::string> args;
      while (i_arg+1 < argc and argv[i_arg+1][0] != '-') {
        args.push_back(argv[++i_arg]);
      }
      if (args.size() != 2) {
        std::cerr << arg << " requires three arguments: DENSITY Z A" << std::endl;
        return 1;
      }
      target_Z       = std::stod(args[0]);
      target_A       = std::stod(args[1]);
    } else {
      std::cout << arg << " is an unrecognized option" << std::endl;
      return 1;
    }
  }

  std::ofstream total_file(total_xsec);
  if (!total_file.is_open()) {
    std::cerr << "File '" << total_xsec << "' was not able to be opened." << std::endl;
    return 2;
  }

  std::ofstream dsdxdtheta_f(dsdxdtheta_fn);
  if (!dsdxdtheta_f.is_open()) {
    std::cerr << "File '" << dsdxdtheta_fn << "' was not able to be opened." << std::endl;
    return 2;
  }

  std::ofstream dsdx_f(dsdx_fn);
  if (!dsdx_f.is_open()) {
    std::cerr << "File '" << dsdx_fn << "' was not able to be opened." << std::endl;
    return 2;
  }

  total_file << "energy,xsec\n";
  dsdxdtheta_f << "x,theta,dsdxdtheta\n";
  dsdx_f << "x,dsdx\n";

  const double epsilon_ = 1.0;
  const double MA = ap_mass;
  const double MA2 = MA*MA;
  const double alphaEW = 1.0 / 137.0;
  const double lepton_mass{0.10566};
  const double lepton_mass_sq = lepton_mass*lepton_mass;
  const double threshold_ = 2*MA;

  double current_energy = min_energy;
  int bar_width = 80;
  int pos = 0;
  while (current_energy < max_energy + energy_step) {
    double lepton_ke = current_energy;
    // the cross section is zero if the lepton does not have enough
    // energy to create an A'
    // the threshold_ can also be set by the user to a higher value
    // to prevent dark-brem within inaccessible regions of phase
    // space
    if (lepton_ke < threshold_) {
      total_file << current_energy << "," << 0 << "\n";
      continue;
    }

    double lepton_e = lepton_ke + lepton_mass;
    double lepton_e_sq = lepton_e*lepton_e;

    /*
     * "Hyper-Improved" WW
     *
     * assume theta = 0, and x = 1 for form factor integration
     * i.e. now chi is a constant pulled out of the integration
     */
    double chi_hiww = flux_factor_chi_numerical(target_A,target_Z,
        MA2*MA2/(4*lepton_e_sq),MA2+lepton_mass_sq);
  
    /*
     * Differential cross section with respect to x and theta
     *
     * Equation (16) from Appendix A of https://arxiv.org/pdf/2101.12192.pdf
     *
     * This `auto` represents a lambda-expression function, inheriting many
     * pre-calculated constants (like lepton_e and chi) while also calculating
     * the variables dependent on the integration variables. The return value
     * of this function is a double since it is calculated by arithmetic
     * operations on doubles.
     */
    auto diff_cross = [&](double x, double theta) {
      if (x*lepton_e < threshold_) return 0.;
  
      double theta_sq = theta*theta;
      double x_sq = x*x;
  
      double utilde = -x*lepton_e_sq*theta_sq - MA2*(1.-x)/x - lepton_mass_sq*x;
      double utilde_sq = utilde*utilde;
  
      /*
       * WW
       *
       * Since muons are so much more massive than electrons, we keep 
       * the form factor integration limits dependent on x and theta
       */
  
      // non-zero theta and non-zero m_l
      //double tmin = utilde_sq/(4.0*lepton_e_sq*(1.0-x)*(1.0-x));
      // maximum t kinematically limited to the incident lepton energy
      //double tmax = lepton_e_sq;
  
      /*
       * The chi integrand limits given by
       *
       * Eqs (3.20) and (A6) of
       * https://journals.aps.org/prd/pdf/10.1103/PhysRevD.8.3109
       * OR
       * Eqs (3.2) and (3.6) of 
       * https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.46.815
       *
       * to be
       *
       * tmax = m^2(1+l)^2
       * tmin = tmax / (2*E*x*(1-x))^2
       *
       * where
       *
       *  l = E^2x^2theta^2/m^2
       *  m is mass of dark photon
       *  E is the incident lepton energy
       * 
       * were investigated in an attempt to control the numerical integration
       * of chi in the hopes that cutting the integral away from odd places
       * would be able to avoid the funky business. This was not successful,
       * but we are leaving them here in case a typo is found in the future
       * or the search is chosen to resume.
       */
      double el = lepton_e_sq*x_sq*theta_sq/MA2;
      double tmax = MA2*pow(1 + el,2);
      double tmin = MA2*tmax / pow(2*lepton_e*x*(1-x),2);
    
      // require 0 < tmin < tmax to procede
      if (tmin < 0) return 0.;
      if (tmax < tmin) return 0.;
    
      /*
       * numerically integrate to calculate chi ourselves
       * this _has not_ been well behaved due to the extreme values
       * of t that must be handled
      double chi = flux_factor_chi_numerical(A,Z, tmin, tmax);
       */
    
      /*
       * use analytic elastic-only chi derived for DMG4
       * and double-checked with Mathematica
       *
       * The inelastic integral contains some 4000 terms
       * according to Mathematica so it is expensive to
       * compute and only an O(few) percent change.
       */
      double chi_analytic_elastic_only = flux_factor_chi_analytic(target_A,target_Z,tmin,tmax);
      
      /*
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
               *(chi_analytic_elastic_only/utilde_sq)*amplitude_sq*sin(theta);
    };

    // deduce integral bounds
    double xmin = 0;
    double xmax = 1;
    if ((lepton_mass / lepton_e) > (MA / lepton_e))
      xmax = 1 - lepton_mass / lepton_e;
    else
      xmax = 1 - MA / lepton_e;
  
    /*
     * max recoil angle of A'
     *
     * The wide angle A' are produced at a negligible rate
     * so we enforce a hard-coded cut-off to stay within
     * the small-angle regime.
     *
     * We choose the same cutoff as DMG4.
     */
    double theta_max{0.3};
  
    /*
     * Integrand for integral over x
     *
     * For muons, we want to include the variation over theta from the chi
     * integral, so we calculate the x-integrand by numerically integrating
     * over theta in the differential cross section defined above.
     *
     * For electrons, we are using the Improved WW method where the theta
     * integral has already been done analytically and we can use the
     * numerical Chi (including both inelastic and elastic form factors)
     * calculated above.  
     *
     * This is the final lambda expression used here. Its one argument is a double
     * and it returns a double.
     */
    auto theta_integral = [&](double x) {
      auto theta_integrand = [&](double theta) {
        double dsdxdtheta = diff_cross(x, theta);
        dsdxdtheta_f << x << "," << theta << "," << dsdxdtheta << "\n";
        return dsdxdtheta;
      };
      //theta_max = 3*(1-x)/lepton_e;
      // integrand, min, max, max_depth, tolerance, error, pL1
      double dsdx = integrate(theta_integrand, 0., theta_max);
      dsdx_f << x << "," << dsdx << "\n";
      return dsdx;
      /*
       * electron fork which we are ignoring while debugging muon
      if (muons) {
      } else {
        if (x*lepton_e < threshold_) {
          diff_file << x << "," << 0. << "\n";
          return 0.;
        }
        double beta = sqrt(1 - MA2/lepton_e_sq),
               nume = 1. - x + x*x/3.,
               deno = MA2*(1-x)/x + lepton_mass_sq;
        double dsdx = 4*pow(epsilon_,2)*pow(alphaEW,3)*chi_hiww*beta*nume/deno;
        return dsdx;
      }
       */
    };
  
    double error;
    double integrated_xsec = integrate(theta_integral, xmin, xmax);
  
    double GeVtoPb = 3.894E08;
    /*
     * The integrated_xsec should be the correct value, we are just
     * converting it to Geant4's pb units here
     */
    double cross = integrated_xsec * GeVtoPb;
    total_file << current_energy << "," << cross << "\n";

    current_energy += energy_step;
    int old_pos{pos};
    pos = bar_width * current_energy / max_energy;
    if (pos != old_pos) {
      std::cout << "[";
      for (int i{0}; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "] " << int(current_energy / max_energy * 100.0) << " %\r" << std::flush;
    }
  }
  std::cout << std::endl;

  total_file.close();
  dsdxdtheta_f.close();
  dsdx_f.close();

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
