#ifndef ELECTRON_H_
#define ELECTRON_H_

#include <random>

#include "photon.h"
#include "my_random.h"

class electron
{
    ////////// Initializer //////////
public:
     electron(  );
    ~electron(  );

    ////////// Random related //////////
private:			// Functors
    static std::default_random_engine     generator;
    std::uniform_real_distribution<double>  mu_rand;
    std::uniform_real_distribution<double> phi_rand;
    
    ////////// ( E, p ) distribution //////////
private:			// Data
    static const double theta_te; // From T_e to Theta
    double px, py, pz;
private:			// Functor
    my_random * 
private:			  // Function
    void get_momentum( const double & t_e );

    ////////// Scatter a photon //////////
public:				// Function
    void scatter_ph( photon & ph );
};
    
#endif
