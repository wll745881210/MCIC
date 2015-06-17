#include "rand_knscat.h"

#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////
// Static variables

rand_knscat * rand_knscat::singleton( nullptr );

////////////////////////////////////////////////////////////
// Initializer

rand_knscat * rand_knscat::get_instance(  )
{
    if( singleton == nullptr )
	singleton = new rand_knscat;
    return singleton;
}

void rand_knscat::del_instance(  )
{
    if( singleton != nullptr )
	delete singleton;
    return;
}

void rand_knscat::set_resolution( const int & n_mu )
{
    x_vec.clear(  );
    
    const double dmu = 2. / ( n_mu - 1 );
    for( int i = 0; i < n_mu; ++ i )
	x_vec.push_back( -1. + dmu * i );
    return;
}

void rand_knscat::set_eta_range
( const double & eta0, const double & eta1,
  const int    & n_eta )
{
    t_vec.clear(  );
    
    const double t0 = log( eta0 );
    const double t1 = log( eta1 );
    const double dt = ( t1 - t0 ) / ( n_eta - 1 );
    for( int i = 0; i < n_eta; ++ i )
	t_vec.push_back( t0 + dt * i );

    return;
}

////////////////////////////////////////////////////////////
// Override the integration-related things

double rand_knscat::compton_factor
( const double & e, const double & mu )
{
    return 1. / ( 1. + e * ( 1. - mu ) );
}

double rand_knscat::pdf( const double & x )
{
    const double cf = compton_factor( eta, x );
    return pow( cf, 2 ) * ( cf + 1. / cf - 1. + x * x )
	/ norm;
}

std::map<double, double> *
rand_knscat::intg_single_t( const double & t )
{
    eta	      = exp( t );
    c_current = 0.;

    // Normalization comes from Mathematica...
    norm = ((2*eta*(2 + eta*(1 + eta)*(8 + eta)))
	    /pow(1 + 2*eta,2) + (-2 + (-2 + eta)*eta)
	*log(1 + 2*eta))/pow(eta,3);
    
    auto res = rand_base::intg_single_t( t );
    res->insert( std::make_pair( 1., 1. ) );
    return res;
}

void rand_knscat::prepare_intg(  )
{
    return;
}

////////////////////////////////////////////////////////////
// Final output

double rand_knscat::get_rand_mu( const double & eta )
{
    return rand_base::get_rand( log( eta ) );
}
