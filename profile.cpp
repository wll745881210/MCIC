#include "profile.h"

#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
// Static variables

profile * profile::singleton( nullptr );

////////////////////////////////////////////////////////////
// Initializers

profile::profile(  )
{
    return;
}

profile::~profile(  )
{
    return;
}

profile * profile::get_instance(  )
{
    if( singleton == nullptr )
	singleton = new profile;
    return singleton;
}

void profile::del_instance(  )
{
    if( singleton != nullptr )
	delete singleton;
    return;
}

void profile::init( input & args )
{
    // Since n = const in r_core, and due to normalization,
    // tau_core and r_core are identical.
    args.find_key( "tau_core",    r_core,      1e-1 );
    args.find_key( "tau_max",     tau_max,     1.   );
    args.find_key( "theta_r_max", theta_r_max, 0.5  );
    args.find_key( "n_pow",       n_pow,       0.   );
    args.find_key( "theta_pow",   theta_pow,   0.   );
    n_core = 1.;		// Definition...    

    obtain_rmax(  );
    return;
}

////////////////////////////////////////////////////////////
// Get the values.
// If you want to use it for further purposes, modify here.

void profile::obtain_rmax(  )
{
    double tau( 0. );
    std::array<double, 3> x = {};

    double rho_ratio0 = rho_ratio( x );
    static const double d_tau( 1e-2 );
    while( tau < tau_max )
    {
	// Predictor
	const double r_current = x[ 0 ];
	x[ 0 ] += d_tau / rho_ratio0;
	const double rho_ratio1 = rho_ratio( x );
	// Corrector
	rho_ratio0 = 0.5 * ( rho_ratio0 + rho_ratio1 );
	x[ 0 ] = r_current + d_tau / rho_ratio0;
	// Save for the next step
	rho_ratio0 = rho_ratio1;
	tau += d_tau;
    }
    r_max = x[ 0 ];
    return;
}


double profile::radius( const std::array< double, 3 > & x )
{
    double r2( 0. );
    for( int i = 0; i < 3; ++ i )
	r2 += x[ i ] * x[ i ];
    return std::max( sqrt( r2 ), r_core );
}

double profile::rho_ratio
( const std::array< double, 3 > & x )
{
    const double r = radius( x );
    return n_core * pow( r / r_core, n_pow );
}

double profile::theta( const std::array< double, 3 > & x )
{
    const double r = radius( x );
    return theta_r_max * pow( r / r_max, theta_pow );
}

double profile::get_rmax(  )
{
    return r_max;
}
