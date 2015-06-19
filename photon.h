#ifndef PHOTON_H_
#define PHOTON_H_

#include <random>
#include <vector>
#include <array>

#include "rand_planck.h"
#include "electron.h"
#include "profile.h"
#include "input.h"

class photon
{
    ////////// Initializer //////////
private:			// Data
    static double theta_bb;
public:
     photon(  );
    ~photon(  );
    static void init( input & args  );

    ////////// Location and momentum //////////
private:			// Data
    static double r_max;
    std::array< double, 3 > x;	// 3-location
    std::array< double, 4 > p;	// 4-momentum
private:			// Function
    void init_loc(  );		// initial position
    void init_mom(  );		// initial momentum
    void reset(  );
    double radius_c(  );	// reduced r from the center

    ////////// Random related //////////
private:			// Functors
    std::default_random_engine     generator;
    std::exponential_distribution <double> exp_rand;
    std::uniform_real_distribution<double> uni_rand;
    
    ////////// Iteration //////////
private:			// Data
    static int scat_max;
    static int n_repeat;
    int n_itr;
    static double d_tau_fiducial;
    bool continue_walking;
private:			// Functions & functors
    electron elec;
    void step_walk( const double & d_tau );
public:				// Function
    void iterate_photon(  );

    ////////// Utilities //////////
private:
    profile * prof;
    
    ////////// Data dump and access //////////
private:			// Data
    std::vector<double> res;
public:				// Data access
    const std::vector<double> & get_res(  );
};

#endif
