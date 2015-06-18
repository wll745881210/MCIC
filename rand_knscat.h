#ifndef RAND_KNSCAT_H_
#define RAND_KNSCAT_H_

#include "rand_base.h"
#include "input.h"

#include <random>
#include <map>

////////////////////////////////////////////////////////////
// Klein-Nishina scatterings with proper scattering angle
// (and hence energy) distribution.
// Here the parameter t is defined as,
// t := log( eta ); eta := h \nu_0 / ( m_e c^2 ).
// The random variable is simply the scattering cosine mu,
// x := \mu.

class rand_knscat : public rand_base
{
    ////////// Initializers //////////
private:
    static rand_knscat * singleton;
public:
    static rand_knscat * get_instance(  );
    static void          del_instance(  );
    void init( input & args );

    ////////// PDF related //////////
private:			// Data
    double eta, norm;
private:			// Function
    inline double compton_factor
    ( const double & e, const double & mu );    
    double pdf( const double & x );

    ////////// Integration related //////////
private:
    std::map<double, double> *
    intg_single_t( const double & t );
    void prepare_intg(  );
    
public:				// Function
    double get_rand_mu( const double & eta );    
};



#endif
