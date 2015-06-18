#ifndef PHOTON_H_
#define PHOTON_H_

#include <random>
#include <vector>
#include <array>

#include "rand_planck.h"
#include "profile.h"

class electron;

class photon
{
    ////////// Initializer //////////
private:			// Data
    static double theta_bb;
public:
     photon(  );
    ~photon(  );
    void reset(  );
    static void set_theta_bb( const double & t_bb  );
    static void set_max_r   ( const double & s_max );
    static void set_max_scat( const int    & i_sca );
    static void set_n_repeat( const int    & i_rep );

    ////////// Location and momentum //////////
private:			// Data
    static double r_max;
    std::array< double, 3 > x;	// 3-location
    std::array< double, 4 > p;	// 4-momentum
private:			// Function
    void init_loc(  );		// initial position
    void init_mom(  );		// initial momentum
    double radius_c(  );	// reduced r from the center
public:				// Functor
    friend class electron;

    ////////// Random related //////////
private:			// Functors
    static std::default_random_engine     generator;
    std::exponential_distribution <double> exp_rand;
    std::uniform_real_distribution<double> uni_rand;
    
    ////////// Iteration //////////
private:			// Data
    static int scat_max;
    static int n_repeat;
    int n_itr;
    double d_tau_fiducial;
    bool continue_walking;
public:			// Function
    void step_walk( const double & d_tau );
    void iterate  (  );
public:				// Function
    void proceed_photon(  );

    ////////// Utilities //////////
private:
    profile * prof;
    
    ////////// Statistics //////////
private:			// Data
    std::vector<double> res;
public:				// Data access
    
};

#endif
