#include "profile.h"

#include <iostream>
#include <algorithm>
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
    args.find_key( "r_core",     r_core,      1e2  );
    args.find_key( "theta_norm", theta_norm,  0.5  );
    args.find_key( "theta_cap",  theta_cap,   0.5  );
    args.find_key( "n_pow",      n_pow,      -1.25 );
    n_core = 1.;		// Definition...    

    return;
}

////////////////////////////////////////////////////////////
// Get the values.
// If you want to 

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
    const double r   = radius( x );
    const double n_r = n_core * pow( r / r_core, n_pow );
    return std::min( n_core, n_r );
}

double profile::theta( const std::array< double, 3 > & x )
{
    const double r   = radius( x );
    const double t_r = theta_norm / ( r / r_core );
    return std::min( theta_cap, t_r );
}
