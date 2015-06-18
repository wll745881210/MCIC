#include "photon.h"
#include "rand_planck.h"
#include "profile.h"

#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////
// Static variables
double photon::theta_bb( 2e-6 ); // default: k_B T_bb = 1eV
double photon::r_max   ( 1.   );
int    photon::scat_max ( 20   );
int    photon::n_repeat( 20   );
std::default_random_engine photon::generator;

////////////////////////////////////////////////////////////
// Initializer

photon::photon(  ) : exp_rand( 1. ), uni_rand( 0, 1 )
{
    d_tau_fiducial = 0.01;
    this->reset(  );
    return;
}

photon::~photon(  )
{
    return;
}

void photon::init_loc(  )
{
    x = y = z = 0.;
    return;
}

void photon::init_mom(  )
{
    auto p_planck    = rand_planck::get_instance(  );
    const double e   = theta_bb * p_planck->get_rand(  );

    // Not necessary, but let it be here for future
    // extension onto non-spherical cases.
    const double mu  = -1 + 2. * uni_rand( generator );
    const double cmu = sqrt( 1. - mu * mu );
    const double phi = uni_rand( generator ) * 6.28318531;
    
    p[ 0 ] = e;
    p[ 1 ] = e * cmu * cos( phi );
    p[ 2 ] = e * cmu * sin( phi );
    p[ 3 ] = e * mu;
    
    return;
}

void photon::reset(  )
{
    this->init_loc(  );
    this->init_mom(  );
    return;
}

void photon::set_theta_bb( const double & t_bb )
{
    photon::theta_bb = t_bb;
    return;
}

void photon::set_max_r( const double & s_max )
{
    photon::r_max = s_max;
    return;
}

void photon::set_max_scat( const int & i_sca )
{
    photon::scat_max = i_sca;
    return;
}

void photon::set_n_repeat( const int & i_rep )
{
    photon::n_repeat = i_rep;
    return;
}

////////////////////////////////////////////////////////////
// Locations

double photon::radius_c(  )
{
    double r2( 0. );
    for( int i = 0; i < 3; ++ i )
	r2 += x[ i ] * x[ i ];
    return sqrt( r2 );
}

////////////////////////////////////////////////////////////
// Iteration

void photon::step_walk( const double & tau )
{
    double tau_gone( 0. );
    const double & e = p[ 0 ];

    auto p = profile::get_instance(  ) ;

    double rho_ratio = p->rho_ratio( x );
    double d_tau     = d_tau_fiducial;
    
    while( tau_gone < tau )
    {
	const double d_x = d_tau / rho_ratio;
	for( int i = 0; i < 3; ++ i )
	    x[ i ] += d_x * p[ i + 1 ] / e;
	tau_gone += d_tau;
	
	rho_ratio = p->rho_ratio( x );
	if( tau - tau_gone < d_tau )
	    d_tau = fabs( tau - tau_gone ) * 1.0001;
	// 1.0001: just to make sure that the loop will end
    }

    return;
}

bool photon::iterate(  )
{
    for( int i = 0; i < itr_max; ++ i )
    {
	const double tau = exp_rand( generator );
	step_walk( tau );
	if( this->radius_c(  ) > r_max )
	    break;
    }
    res.push_back(  )
}

void photon::walk_photon(  )
{
    res.resize( itr_max, 0. );
    
    for( int i = 0; i < n_repeat; ++ i )
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
