#ifndef RAND_BASE_H_
#define RAND_BASE_H_

#include <random>
#include <vector>
#include <map>

////////////////////////////////////////////////////////////
// Random variable generator with one parameter.
// Conventions:
// x for the target random variable;
// t for the parameter;
// c for cumulative distribution function (cdf).

class rand_base
{
    ////////// Initializers //////////
public:
     rand_base(  );
    ~rand_base(  );

    ////////// Distribution function //////////
private:			// PDF and conversion
    virtual double pdf( const double & x ) = 0;
    // x := exp( - ( gamma - 1 ) / theta ).

    ////////// Uniform random generator //////////
protected:			// Functor
    std::default_random_engine     generator;
    std::uniform_real_distribution<double> cdf_rand;
    
    ////////// Integration //////////
protected:			// Data
    std::vector<double> x_vec, t_vec;
    std::map<double, std::map<double, double> * > t_map;
    double c_current;
private:			// Function
    void rk4( double & p, const double & x,
	      const double & dx );
protected:
    virtual std::map<double, double> *
    intg_single_t( const double & t );
    virtual void prepare_intg(  );
public:
    void integrate(  );

    ////////// Interpolation //////////
protected:			// Function
    double interp_single_t
    ( const double & c, std::map<double, double> * p );
    double interp_cdf
    ( const double & c, const double & t );
public:				// Function
    double get_rand( const double & t );
};

#endif
