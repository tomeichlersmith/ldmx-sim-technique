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
 *  - larger depth of 15 gave same results as depth 5
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
    "  --chi        : output file to write chi integrand values to\n"
    "                 WARNING: this file gets very big if doing multiple energy samples\n"
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
 * numerically integrate the value of the flux factor chi
 *
 * The integration of the form factor into the flux factor can
 * be done analytically with a tool like mathematica, but when
 * including the inelastic term, it produces such a complicated 
 * result that the numerical integration is actually *faster*
 * than the analytical one.
 *
 * The form factors are copied from Appendix A (Eq A18 and A19) of
 * https://journals.aps.org/prd/pdf/10.1103/PhysRevD.80.075018
 *
 * Here, the equations are listed for reference.
 * \f{equation}{
 * \chi(x,\theta) = \int^{t_{max}}_{t_{min}} dt \left( \frac{Z^2a^4t^2}{(1+a^2t)^2(1+t/d)^2}+\frac{Za_p^4t^2}{(1+a_p^2t)^2(1+t/0.71)^8}\left(1+\frac{t(\mu_p^2-1)}{4m_p^2}\right)^2\right)\frac{t-t_{min}}{t^2}
 * \f}
 * where
 * \f{equation}{
 * a = \frac{111.0}{m_e Z^{1/3}}
 * \quad
 * a_p = \frac{773.0}{m_e Z^{2/3}}
 * \quad
 * d = \frac{0.164}{A^{2/3}}
 * \f}
 * - \f$m_e\f$ is the mass of the electron in GeV
 * - \f$m_p = 0.938\f$ is the mass of the proton in GeV
 * - \f$\mu_p = 2.79\f$ is the proton \f$\mu\f$
 * - \f$A\f$ is the atomic mass of the target nucleus in amu
 * - \f$Z\f$ is the atomic number of the target nucleus
 *
 * @param[in] A atomic mass of the target nucleus in amu
 * @param[in] Z atomic number of target nucleus
 * @param[in] tmin lower limit of integration over t
 * @param[in] tmax upper limit of integration over t
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
   * whose return value is a `double` and which has a single input `lnt`. 
   * This lambda expression saves us the time of having to re-calculate 
   * the form factor constants that do not depend on `t` because it 
   * can inherit their values from the environment. 
   * The return value is a double since it is calculated
   * by simple arithmetic operations on doubles.
   *
   * The integrand is so sharply peaked at t close to tmin,
   * it is very helpful to do the integration in the variable
   * u = ln(t) rather than t itself.
   */
  auto integrand = [&](double lnt) {
    double t = exp(lnt);
    double ael_factor = 1./(ael_inv2 + t),
           del_factor = 1./(1+t/del),
           ain_factor = 1./(ain_inv2 + t),
           din_factor = 1./(1+t/din),
           nucl = (1 + t*bin);
    
    return (pow(ael_factor*del_factor*Z, 2)
            + Z*pow(ain_factor*nucl*din_factor*din_factor*din_factor*din_factor, 2)
           )*(t-tmin)*t;
  };

  return integrate(integrand,log(tmin),log(tmax));
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
  std::string chi_fn{"/dev/null"};
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
    } else if (arg == "--chi") {
      if (i_arg+1 >= argc) {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      chi_fn = argv[++i_arg];
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

  std::ofstream chi_f(chi_fn);
  if (!chi_f.is_open()) {
    std::cerr << "File '" << chi_fn << "' was not able to be opened." << std::endl;
    return 2;
  }

  total_file << "energy,xsec\n";
  dsdxdtheta_f << "x,theta,dsdxdtheta\n";
  dsdx_f << "x,dsdx\n";
  chi_f << "x,theta,tmin,tmax,chi\n";

  const double epsilon_ = 1.0;
  const double MA = ap_mass;
  const double MA2 = MA*MA;
  const double alphaEW = 1.0 / 137.0;
  const bool muons = false; //true;
  const double lepton_mass{0.000511}; //0.10566};
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

    double xmax = 1 - std::min(lepton_mass,MA) / lepton_e;

    double integrated_xsec{-1};
    if (true) { // Full WW
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
      //double theta_max{2*3.14159};

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
    
        double x_sq = x*x;
    
        double utilde = -x*lepton_e_sq*theta*theta - MA2*(1.-x)/x - lepton_mass_sq*x;
        //double utilde = -x*lepton_e_sq*pow(sin(theta),2) - MA2*(1.-x)/x - lepton_mass_sq*x;
        //double utilde = -x*lepton_e_sq*2*(1 - cos(theta)) - MA2*(1.-x)/x - lepton_mass_sq*x;
        double utilde_sq = utilde*utilde;
    
        // non-zero theta and non-zero m_l
        double tmin = utilde_sq/(4.0*lepton_e_sq*(1.0-x)*(1.0-x));
        // maximum t kinematically limited to the incident lepton energy
        double tmax = lepton_e_sq;

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
        double el = lepton_e_sq*x_sq*theta_sq/MA2;
        double tmax = MA2*pow(1 + el,2);
        double tmin = MA2*tmax / pow(2*lepton_e*x*(1-x),2);
        tmax = lepton_e_sq;
         */
  
        // require 0 < tmin < tmax to procede
        if (tmin < 0 or tmax < tmin) {
          return 0.;
        }
      
        double chi = flux_factor_chi_numerical(target_A,target_Z,tmin, tmax);
        chi_f << x << "," << theta << "," << tmin << "," << tmax << "," << chi << "\n";
      
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
                 *(chi/utilde_sq)*amplitude_sq*sin(theta);
      };

      integrated_xsec = integrate(
          [&](double x) {
            auto theta_integrand = [&](double theta) {
              double dsdxdtheta = diff_cross(x, theta);
              dsdxdtheta_f << x << "," << theta << "," << dsdxdtheta << "\n";
              return dsdxdtheta;
            };
            double dsdx = integrate(theta_integrand, 0., theta_max);
            dsdx_f << x << "," << dsdx << "\n";
            return dsdx;
          }, 0, xmax);
    } else if (false) { // Improved WW
      /*
       * do the theta integral analytically by neglecting
       * all theta terms in the integrand.
       */
      integrated_xsec = integrate(
          [&](double x) {
            if (x*lepton_e < threshold_) return 0.;
            double utilde = -MA2*(1.-x)/x -lepton_mass_sq*x;
            double utilde_sq = utilde*utilde;
            // non-zero theta and non-zero m_l
            double tmin = utilde_sq/(4.0*lepton_e_sq*(1.0-x)*(1.0-x));
            // maximum t kinematically limited to the incident lepton energy
            double tmax = lepton_e_sq;
            // require 0 < tmin < tmax to procede
            if (tmin < 0) return 0.;
            if (tmax < tmin) return 0.;
            double chi = flux_factor_chi_numerical(target_A,target_Z,tmin,tmax);
            chi_f << x << "," << 0 << "," << tmin << "," << tmax << "," << chi << "\n";
            double beta = sqrt(1 - MA2/lepton_e_sq),
                   nume = 1. - x + x*x/3.,
                   deno = MA2*(1-x)/x + lepton_mass_sq;
            double dsdx = 4*pow(epsilon_,2)*pow(alphaEW,3)*chi*beta*nume/deno;
            dsdx_f << x << "," << dsdx << "\n";
            return dsdx;
          }, 0, xmax);
    } else if (false) { // HyperImproved WW
      /*
       * calculate chi once at x=1, theta=0 and then use it
       * everywhere in the integration over dsigma/dx
       *
       * cut off the integration earlier than the lepton energy squared
       * so that this overestimate isn't too much of an overestimate.
       */
      double tmin = MA2*MA2/(4*lepton_e_sq);
      double tmax = MA2+lepton_mass_sq;
      double chi_hiww = flux_factor_chi_numerical(target_A,target_Z,tmin,tmax);
      chi_f << 1 << "," << 0 << "," << tmin << "," << tmax << "," << chi_hiww << "\n";
  
      integrated_xsec = integrate(
          [&](double x) {
            if (x*lepton_e < threshold_) return 0.;
            double beta = sqrt(1 - MA2/lepton_e_sq),
                   nume = 1. - x + x*x/3.,
                   deno = MA2*(1-x)/x + lepton_mass_sq;
            double dsdx = 4*pow(epsilon_,2)*pow(alphaEW,3)*chi_hiww*beta*nume/deno;
            dsdx_f << x << "," << dsdx << "\n";
            return dsdx;
          }, 0, xmax);
    } else {
      std::cerr << "Did not choose a cross section calculation method.";
      return -127;
    }

    static const double GeVtoPb = 3.894E08;

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
  chi_f.close();

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}
