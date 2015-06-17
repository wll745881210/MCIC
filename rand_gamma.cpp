#include "rand_gamma.h"

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

////////////////////////////////////////////////////////////
// Static variables

rand_gamma * rand_gamma::singleton( nullptr );

////////////////////////////////////////////////////////////
// Initializers

rand_gamma * rand_gamma::get_instance(  )
{
    if( singleton == nullptr )
	singleton = new rand_gamma;
    return singleton;
}

void rand_gamma::del_instance(  )
{
    if( singleton != nullptr )
	delete singleton;
    return;
}

void rand_gamma::set_resolution( const int & n_x )
{
    x_vec.clear(  );
    
    const double dx = 1. / ( n_x - 1 );
    for( int i = 0; i < n_x; ++ i )
	x_vec.push_back( 1. - dx * i );
    return;
}

void rand_gamma::set_theta_range
( const double & theta0, const double & theta1,
  const int    & n_theta )
{
    t_vec.clear(  );
    
    const double t0 = log( theta0 );
    const double t1 = log( theta1 );
    const double dt = ( t1 - t0 ) / ( n_theta - 1 );
    for( int i = 0; i < n_theta; ++ i )
	t_vec.push_back( t0 + dt * i );

    return;
}

////////////////////////////////////////////////////////////
// Override the integration-related things

double rand_gamma::pdf( const double & x )
{
    const double lnx = log( x );
    return norm_theta * ( 1 - theta * lnx ) *
	sqrt( pow( theta * lnx, 2 ) - 2 * theta * lnx );
}

std::map<double, double> *
rand_gamma::intg_single_t( const double & t )
{
    theta      = exp( t );
    ln_theta   = t;
    norm_theta = exp( - 1. / theta )
	/ boost::math::cyl_bessel_k( 2, 1. / theta );

    c_current = 1;
    auto res = rand_base::intg_single_t( t );
    res->insert( std::make_pair( 0, 0 ) );

    return res;
}

void rand_gamma::prepare_intg(  )
{
    return;
}

////////////////////////////////////////////////////////////
// And the final output

double rand_gamma::get_rand_gamma( const double & theta )
{
    const double x = rand_base::get_rand( log( theta ) );
    return 1. - theta * log( x );
}
