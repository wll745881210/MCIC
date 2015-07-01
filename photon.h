#ifndef PHOTON_H_
#define PHOTON_H_

#include <random>
#include <vector>
#include <array>

#include "rand_planck.h"
#include "electron.h"
#include "profile.h"
#include "input.h"
#include "seed.h"

////////////////////////////////////////////////////////////

class photon
{
    ////////// typedef //////////
public:
    typedef unsigned long long uint;
    
    ////////// Initializer //////////
private:			// Data
    static double theta_bb;
public:
     photon(  );
    ~photon(  );
    static void init( input & args );

    ////////// Location and momentum //////////
private:			// Data
    static double r_max,      r_min;
    static double r_disk_max, r_disk_min;
    std::array< double, 3 > x;	// 3-location
    std::array< double, 4 > p;	// 4-momentum
private:			// Function
    void init_loc(  );		// initial position
    void init_mom(  );		// initial momentum
    void reset(  );
    double radius_c(  );	// reduced r from the center

    ////////// Random related //////////
private:			// Functors
    std::mt19937 generator;
    std::exponential_distribution <double> exp_rand;
    std::uniform_real_distribution<double> uni_rand;
    
    ////////// Iteration //////////
private:			// Data
    static int    scat_max;
    static uint   n_repeat;
    static double d_tau_fiducial;
    bool continue_walking;
private:			// Functions & functors
    electron elec;
    void step_walk( const double & d_tau );
    void locate_bin( const double & eta );
public:				// Function
    void iterate_photon(  );

    ////////// Utilities //////////
private:
    static profile * prof;
    
    ////////// Data dump and access //////////

private:			// Data
    static std::map<double, uint> bin_map;
    static std::vector<double> eta_upper;
    // Upper bound of eta
    std::vector<uint> res;
public:				// Data access
           const std::vector<uint  > & get_res(  );
    static const std::vector<double> & get_eta_upper(  );
};

#endif
