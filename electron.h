#ifndef ELECTRON_H_
#define ELECTRON_H_

#include <random>
#include <array>

#include "rand_gamma.h"
#include "rand_knscat.h"
#include "seed.h"

class electron
{
    ////////// Initializer //////////
public:
     electron(  );
    ~electron(  );

    ////////// Random related //////////
private:			// Functors
    std::mt19937 generator;
    std::uniform_real_distribution<double> uni_rand;
    
    ////////// Scatterings //////////
public:			  // Function
    void lorentz_trans
    ( const double & beta, const double & gamma,
      const double & mu  , const double & phi,
      std::array<double, 4> & v );
    void rotate_back_ph
    ( std::array< double, 4 > & p_old,
      std::array< double, 4 > & p_new );
public:				// Function
    void scatter_ph( std::array< double, 4 > & p_ph,
	             const double & theta );
};
    
#endif
