#ifndef ELECTRON_H_
#define ELECTRON_H_

#include <random>
#include <array>

#include "rand_gamma.h"
#include "rand_knscat.h"

class electron
{
    ////////// Initializer //////////
public:
     electron(  );
    ~electron(  );

    ////////// Random related //////////
private:			// Functors
    std::default_random_engine     generator;
    std::uniform_real_distribution<double> uni_rand;
    
    ////////// Scatterings //////////
private:			  // Function
    void lorentz_trans
    ( const double & beta, const double & gamma,
      const double & mu  , std::array<double, 4> & v );
    double get_mu_2d( const std::array< double, 4 > v );
    void rotate_back_mu_2d
    ( std::array< double, 4 > & v, const double & mu );
public:				// Function
    void scatter_ph( std::array< double, 4 > & p_ph,
	             const double & theta );
};
    
#endif
