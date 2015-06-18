#include "rand_planck.h"

#include <iostream>

////////////////////////////////////////////////////////////
// Static variables

rand_planck * rand_planck::singleton( nullptr );
const double  rand_planck::epsilon  ( 1e-5    );

////////////////////////////////////////////////////////////
// Initializers

rand_planck * rand_planck::get_instance(  )
{
    if( singleton == nullptr )
	singleton = new rand_planck;
    return singleton;
}

void rand_planck::del_instance(  )
{
    if( singleton != nullptr )
	delete singleton;
    return;
}

void rand_planck::init( input & args )
{
    int n_x( 0 );
    args.find_key( "planck_res", n_x, 100 );
    x_vec.clear(  );    
    const double dx = 1. / ( n_x - 1 );
    const int max_x = 20;	// "Magic" -> precision lim
    for( int i = 0; i < max_x * n_x; ++ i )
	x_vec.push_back( max_x - i * dx );
    t_vec.push_back( 0 );	// Just let it go.
    return;
}

////////////////////////////////////////////////////////////
// Override the integration-related things

double rand_planck::pdf( const double & x )
{
    static const double norm = 2.404113806319189;
    // 2 * zeta( 3 )
    return x * x / ( exp( x ) - 1 ) / norm;
}

std::map<double, double> *
rand_planck::intg_single_t( const double & t )
{
    c_current = 1.;
    res = rand_base::intg_single_t( t );
    res->insert( std::make_pair( 0., 0. ) );
    return res;
}

void rand_planck::prepare_intg(  )
{
    return;
}

////////////////////////////////////////////////////////////
// And the final output

double rand_planck::get_rand(  )
{
    const double c = cdf_rand( generator );
    return interp_single_t( c, res );
}
