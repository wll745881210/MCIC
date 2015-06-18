#ifndef PHOTON_H_
#define PHOTON_H_

#include <random>
#include <vector>
#include <array>

#include "rand_planck.h"

class electron;

class photon
{
    ////////// Initializer //////////
public:
     photon(  );
    ~photon(  );
    void reset(  );
    static void set_max_tau( const double   & t_max  );
    static void set_max_itr( const int      & i_max  );
    static void set_n_walk ( const int      & i_walk );

    ////////// Location and momentum //////////
private:			// Data
    static double tau_max;
    double  x,  y,  z;
    std::array< double, 4 > p;	// Momentum
private:			// Function
    void set_init_p(  );	// Initial momentum
    double tau_c(  );		// tau from the center
public:				// Function
    friend class electron;

    ////////// Random related //////////
private:			// Functors
    static std::default_random_engine     generator;
    std::exponential_distribution <double> exp_rand;
    
    ////////// Iteration //////////
private:			// Data
    static int itr_max;
    static int n_walk;
    int n_itr;
private:			// Function
    bool photon_step(  );
    void photon_iterate_internal(  );
public:				// Function
    void walk_photon(  );

    ////////// Statistics //////////
private:			// Data
    std::vector<double> res;
public:				// Data access
    static void write_res( std::vector<photon> & ph_list );
};

#endif
