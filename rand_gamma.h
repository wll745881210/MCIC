#ifndef RAND_GAMMA_H_
#define RAND_GAMMA_H_

#include "rand_base.h"

#include <random>
#include <map>

////////////////////////////////////////////////////////////
// Here the distribution function is integrated to give the
// cumulative distribution, by which a uniform distribution
// is mapped onto desired distribution. The best way of
// doing this is to get the distribution function for
// x := exp( - ( gamma - 1 ) / theta ).
// The normalized PDF for x would then be
// f_x = exp( - 1 / theta ) / K_2( 1 / theta )
//     * ( 1 - theta * ln( x ) )
//     * ( theta**2 ln**2( x ) - 2 * theta * ln( x ) )^(1/2)
// and P = \int_0^x f_x' dx'.
// Then, integration is conducted from ( P = 1, x = 1 ) to
// ( P = 0, x = 0 ).
// BTW, when theta < 1e-2, we shall reduce our PDF to the NR
// Maxwellian case, but let's keep it as it is at this
// moment.

class rand_gamma
{
    ////////// Initializers //////////
private:
    static rand_gamma * singleton;
     rand_gamma(  );
    ~rand_gamma(  );
public:
    static rand_gamma * get_instance(  );
    static void         del_instance(  );
    void set_intg_pts( const int & n_x );
    void set_theta
    ( const double & theta0, const double & theta1,
      const int    & n_theta );

    ////////// Distribution function //////////
private:			// Data
    double theta, ln_theta, norm_theta;
private:			// PDF and conversion
    inline double pdf( const double & x );
    // x := exp( - ( gamma - 1 ) / theta ).

    ////////// Random generator //////////
private:			// Functor
    std::default_random_engine     generator;
    std::uniform_real_distribution<double> t_rand;
    
    ////////// Integration //////////
private:			// Data
    double ln_theta0, ln_theta1, d_ln_theta;
    int  n_theta;
    double d_x;
    int    n_x;
    std::map<double, std::map<double, double> * > theta_map;
private:			// Function
    void rk4( double & p, const double & x );
    void intg_single_theta(  );
public:
    void integrate(  );

    ////////// Interpolation //////////
private:			// Function
    double interp_x
    ( const double & y, std::map<double, double> * p );
    double cdf( const double & y, const double & theta );
public:				// Function
    double get_rand_gamma( const double & theta );
};

#endif
