#include "photon.h"
#include "rand_planck.h"

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
    prof = profile::get_instance(  );
    this->reset(  );
    return;
}

photon::~photon(  )
{
    return;
}

void photon::init_loc(  )
{
    x = { };
    continue_walking = true;
    return;
}

void photon::init_mom(  )
{
    auto p_planck    = rand_planck::get_instance(  );
    const double e   = theta_bb * p_planck->get_rand(  );

    // Not necessary, but let it be here for future
    // extension onto non-spherical cases.
    const double mu  = 2 * uni_rand( generator ) - 1;
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

    double rho_ratio = prof->rho_ratio( x );
    double d_tau     = d_tau_fiducial;

    const double & eta = p[ 0 ]; // photon energy

    // Klein-Nishina factor from Mathematica... messy...
    const double kn_factor
	= ((2*eta*(2 + eta*(1 + eta)*(8 + eta)))
	    /pow(1 + 2*eta,2) + (-2 + (-2 + eta)*eta)
	*log(1 + 2*eta))/pow(eta,3) / ( 8. / 3. );
    
    while( tau_gone < tau )
    {
	// Prediction
        auto x1 = x;
	for( int i = 0; i < 3; ++ i )
	    x1[ i ] += d_tau / rho_ratio / kn_factor
		* p[ i + 1 ] / eta;
	const double rho_ratio1 = prof->rho_ratio( x1 );
	// Correction
	rho_ratio = 0.5 * ( rho_ratio + rho_ratio1 );
	for( int i = 0; i < 3; ++ i )
	    x[ i ] += d_tau / rho_ratio / kn_factor
		* p[ i + 1 ] / eta;
	// Save for the next step
	rho_ratio = rho_ratio1;
	
	tau_gone += d_tau;
	if( tau - tau_gone < d_tau )
	    d_tau = fabs( tau - tau_gone ) * 1.0001;
	// 1.0001: make sure that the loop won't stuck
	else if( this->radius_c(  ) > r_max )
	{
	    continue_walking = false;
	    break;
	}
    }
    return;
}

void photon::iterate(  )
{
    for( int i = 0; i < scat_max; ++ i )
    {
	const double tau = exp_rand( generator );
	step_walk( tau );
	if( ! continue_walking )
	    break;
	const double theta_e = 
    }
    res.push_back( p[ 0 ] );
    return;
}

void photon::proceed_photon(  )
{
    for( int i = 0; i < n_repeat; ++ i )
    {
	reset(  );
	iterate(  );
    }
    
    return;
}

