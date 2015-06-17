#ifndef RAND_GAMMA_H_
#define RAND_GAMMA_H_

#include <random>
#include <map>

////////////////////////////////////////////////////////////
// Random variable generator with one parameter.

class rand_base
{
    ////////// Initializers //////////
private:
    static rand_gamma * singleton;
     rand_gamma(  );
    ~rand_gamma(  );
public:
    static rand_gamma * get_instance(  );
    static void         del_instance(  );

    ////////// Distribution function //////////
private:			// Data
    double theta, ln_theta, norm_theta;
private:			// PDF and conversion
    virtual double pdf( const double & x ) = 0;
    // x := exp( - ( gamma - 1 ) / theta ).

    ////////// Random generator //////////
private:			// Functor
    std::default_random_engine     generator;
    std::uniform_real_distribution<double> t_rand;
    
    ////////// Integration //////////
private:			// Data
    double ln_theta0, ln_theta1, d_ln_theta;
    int  n_theta;
    std::map<double, std::map<double, double> * > theta_map;
private:			// Function
    void rk4( double & p, const double & x,
	      const double & d_x );
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
