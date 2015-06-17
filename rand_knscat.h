#ifndef RAND_KNSCAT_H_
#define RAND_KNSCAT_H_

#include <random>
#include <map>

////////////////////////////////////////////////////////////
// Klein-Nishina scatterings with proper scattering angle
// (and hence energy) distribution.
// Here the "auxiliary variable" x is defined as,
// x := h \nu_0 / ( m_e c^2 ).

class rand_knscat
{
    ////////// Initializers //////////
private:
    static rand_knscat * singleton;
     rand_knscat(  );
    ~rand_knscat(  );
public:
    static rand_knscat * get_instance(  );
    static void          del_instance(  );
    void set_intg_pts( const int & n_x );
    void set_gamma
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
    double get_rand_knscat( const double & theta );
};


#endif
