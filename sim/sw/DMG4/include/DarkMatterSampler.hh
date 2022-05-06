#ifndef H_MJ1_SAMPLER_H
#define H_MJ1_SAMPLER_H

#undef NDEBUG  // XXX
#include <cassert>

// NOTE: this file will be apparently moved to the library further

// BGN U-random stub __________
#ifndef URANDOM

// Stores random generator state
struct dphmc_URandomState {
    //void * uRandState;
    CLHEP::HepRandomEngine * hepREng;
};
// Generates random number in [0:1]
//#define URANDOM(u) u = (G4UniformRand());
#define URANDOM(u) u = uRandom->hepREng->flat()

const double me = 5.10998950e-4;  // XXX electron mass, GeV

#endif
// END U-random stub ^^^^^^^^^^

class Sampler {
public:
    const double E0, ma;
//protected:
    double E02  ///< initiating particle energy squared, \f$E_0^2\f$
         , ma2  ///< squared A' mass, \f$m_{A'}^2\f$
         , E02ma2  ///< squared ratio of \f$E_0\f$ and \f$m_{A'}\f$
         , ma4  ///< \f$\m_{A'}^4\F$
         , mema2  ///< \f$\frac{me^2}{ma^2}\f$
         , pi2  ///< \f$\pi^2\f$
         , mj2Norm  ///< norm of \f$M_{2,x}(x)\f$
         ;
public:
    Sampler( double E0_, double ma_ ) : E0(E0_), ma(ma_) {
        E02 = E0*E0;
        ma2 = ma*ma;
        E02ma2 = (E0/ma)*(E0/ma);
        ma4 = ma2*ma2;
        mema2 = (me/ma)*(me/ma);
        pi2 = M_PI*M_PI;
        mj2Norm = log(1/mema2)/(2*E02*(1-mema2)*ma2);
    }

    /** This optimistic estimation for a random \f$x\f$ is derived by integrating
     * the majorant function:
     *
     * \f[
     * M_{x,2} (x) = (m_e^2 x + m_{A'}^2 (1-x) )^(-1),
     * \f]
     *
     * which yields the following reverse function for x:
     *
     * \f[
     * x = \frac{1}{2 E_0} \frac{1-(m_e/m_{A'})^{2 u}}{1-(m_e/m_{A'})^{2}}
     * \f]
     *
     * Used to derive X in two-step modified von Neumann method on A' with WW.
     * */
    double mj2_rev_x( double u ) const {
        return ( (1. - pow(mema2, u))
               / (1. - mema2)
               );
    }

    /** The \f$M_{2,x}(x) > \f$\int\limits_{0}{\pi} M_{1,x} d \theta\f$ for
     * $x \in [0:1]$ being thus a "majorant for integrated majorant". */
    double mj2( double x ) const {
        return (1/(2 * E02 * ma2 * x * x))
              / ( (1-x)/(x*x) + mema2)
              ;
    }

    /** Returns value of the "first" x-majorant (2D majorant integrated over theta):
     *
     * \f[
     * M_{x,1} = \frac{\pi^2 }{2 m_{A'}^4 x} /
     *      (( \frac{m_e^2}{m_{A'}}^2 + \frac{1-x}{x^2} )
     *       ( (\frac{m_e^2}{m_{A'}}^2 + \frac{E_0^2}{m_{A'}^2} \pi^2) + \frac{1-x}{x^2} ))
     * \f]
     *
     * */
    double mj1(double x) const {
        double mxx2 = (1-x)/(x*x)
             , factor1 = pi2/(2*ma2*ma2*x)
             , factor2 = (mema2 + mxx2)
             , factor3 = (mema2 + E02ma2 * pi2 + mxx2)
             ;
        return factor1/(factor2*factor3);
    }

    /** Uses \f$M_{x,2}(u)\f$ (defined by `dphmc_aprime_ww_mj2_rev_x()`)
     * to sample \f$\tilde{x}\f$ ("optimistic" \f$x\f$) according to \f$M_{x,1}(u)\f$
     * (defined by `dphmc_aprime_ww_mj1_x()`).
     *
     * Returned \f$\tilde{x}\f$ has to be used further as a parameter to first-order
     * majorant \f$M_{x,theta,1}(x, \theta)\f$ to sample A' kinematics.
     *
     * \warning Since \f$M_{x,1}(u)\f$ and \f$M_{x,2}(u)\f$ are considered to be
     * very close, no limit on iteration is done at this routine. */
    double sample_x_mj1( struct dphmc_URandomState * uRandom
                       , double xMin=0.
                       , double xMax=1.) const {
        assert(xMin < xMax);
        // TODO: code commented out has to be more efficient, but requires
        //       some algebra I'm (RD) is unable to do right away. Need
        //       inverse relation x <- u for M_{x,2} (see imj2(x)).
        // Use reverse transform to derive restrictions on X. This function
        // monotonically increases so min/max mapping has to take place anyway
        //const double rxMin = imj2(xMin)
        //           , rxMax = imj2(xMax)
        //           , rxScale = rxMax - rxMin;
        //assert(fabs(mj2_rev_x(rxMin) - xMin) < 1e-12);
        //assert(fabs(mj2_rev_x(rxMax) - xMax) < 1e-12);
        //assert(rxMin < rxMax);
        // Sample X according to 2nd majorant (M_{2,x}(x))
        double x;
        double m1, m2, probe;
        do {
            // Generate a pair of random numbers
            URANDOM(x);
            URANDOM(probe);
            //x *= rxScale;  // TODO (see above)
            //x += rxMin;
            // Map one uniform random variable from pair to x wrt M_{2,x}(x)
            // using reverse transform method:
            //  x = \frac{1}{2 E_0} \frac{1-(m_e/m_{A'})^{2 u}}{1-(m_e/m_{A'})^{2}}
            x = mj2_rev_x(x);
            if( x < xMin || x > xMax ) continue;  // TODO (see above)
            // Get the desired value M_{1,x}(x) ...
            m1 = mj1(x);
            // ... and current M_{2,x}(x) value.
            m2 = mj2(x);
            if( probe*m2 <= m1 ) {
                return x;
            }
        } while(true);
    }

    /**For given \f$x\f$ returns random \f$theta\f$ distributed wrt
     * \f$M_{1}(x, \theta)\f$:
     * 
     * \f[
     * \theta = \pi \sqrt{\frac{ u \times ((1-x)/x^2 + m_{e}^2/m_{A'}^2) }
     *          { (1-x)/x^2 + m_{e}^2/m_{A'}^2 + E_0^2/m_{A'} \pi^2 (u-1) }}
     * \f]
     * */
    double sample_theta_mj1( double x, struct dphmc_URandomState * uRandom ) const {
        const double mxx2 = (1-x)/(x*x);
        double u;
        URANDOM(u);
        return M_PI*sqrt( u*(mxx2+mema2) / (mxx2 + mema2 + E02ma2*pi2*(1-u)) );
    }

    /** Returns value of \f$M_{1}(x, \theta)\f$ */
    double mj1(double x, double theta) const {
        const double mxx2 = (1-x)/(x*x)
                   , factor1 = ma2*(E02ma2 * theta * theta + mxx2 + mema2)
                   ;
        return theta / (x*factor1*factor1);
    }

    /** X-section function under consideration
     *
     * Note: was modified according to D.Kirpichnikov notes, second time 20.05.2021
     */
    double reference_f( double x, double theta ) const {
        const double xm1 = 1 - x
                   , U = E02*theta*theta*x + ma2*xm1/x + me*me*x
                   , U2 = U*U
                   ;
        return ( (xm1 + x*x/2)/U2
               + (ma2+2*me*me)*xm1*(ma2*xm1 - U*x + me*me*x*x
                                   )/(U2*U2)
               )*sqrt(x*x-ma2/E02)*sin(theta);
    }

    /** Samples point and returns reference value at this point */
    double sample_x_theta( struct dphmc_URandomState * uRandom
                         , double & xRef
                         , double & thetaRef
                         , const double xMin=0.
                         , const double xMax=1.) const {
        // Sample X according to 2nd majorant (M_{2,x}(x))
        double x, theta, probe, reference;
        do {
            // generate a candidate (triplet)
            x = sample_x_mj1(uRandom, xMin, xMax);
            theta = sample_theta_mj1(x, uRandom);
            URANDOM(probe);
            // test triplet
            reference = reference_f(x, theta);
            if( probe*mj1(x, theta) <= reference ) {
                // within a target function -- accept:
                xRef = x;
                thetaRef = theta;
                return reference;
            }
        } while(true);
    }
};

#endif  // H_MJ1_SAMPLER_H
