#include "photon.h"

#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////
// Static variables
double photon::tau_max( 1. );
int    photon::itr_max( 20 );
int    photon::n_walk ( 20 );
std::default_random_engine photon::generator;

////////////////////////////////////////////////////////////
// Initializer

photon::photon(  ) : exp_rand( 1. ),
    mu_rand( -1., 1. ), phi_rand( 0., 6.2831853 )
{
    this->reset(  );
    return;
}

photon::~photon(  )
{
    return;
}

void photon::reset(  )
{
    x	  = 0.;
    y	  = 0.;
    z	  = 0.;
    n_itr = -1;
    return;
}

void photon::set_max_tau( const double & t_max )
{
    tau_max = t_max;
    return;
}

void photon::set_max_itr( const int & i_max )
{
    itr_max = i_max;
    return;
}

void photon::set_n_walk( const int & i_walk )
{
    n_walk = i_walk;
    return;
}

////////////////////////////////////////////////////////////
// Locations

double photon::tau_c(  )
{
    return sqrt( x * x + y * y + z * z );
}

////////////////////////////////////////////////////////////
// Iteration

bool photon::photon_step(  )
{
    const double d_tau = exp_rand( generator );
    const double mu    =  mu_rand( generator );
    const double phi   = phi_rand( generator );

    const double sin_theta = sqrt( 1. - pow( mu, 2 ) );
    
    x += d_tau * sin_theta * cos( phi );
    y += d_tau * sin_theta * sin( phi );
    z += d_tau * mu;

    ++ n_itr;

    if( this->tau_c(  ) > tau_max )
	return true;
    else
	return false;
}


void photon::photon_iterate_internal(  )
{
    while( n_itr < itr_max )
	if( photon_step(  ) )
	    break;
    return;
}

void photon::walk_photon(  )
{
    res.resize( itr_max, 0. );
    
    for( int i = 0; i < n_walk; ++ i )
    {
	reset(  );
	photon_iterate_internal(  );
	res[ n_itr ] += 1.;
    }
    
    for( unsigned i = 0; i < res.size(  ); ++ i )
	res[ i ] /= n_walk;
    
    return;
}


const std::vector<double> & photon::result(  )
{
    return res;
}
