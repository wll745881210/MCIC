#include "photon.h"
#include "rand_planck.h"

#include <cmath>
#include <iostream>
#include <omp.h>

////////////////////////////////////////////////////////////
// Static variables
double photon::theta_bb( 2e-6 ); // default: k_B T_bb = 1eV
double photon::r_max   ( 1.   );
double photon::d_tau_fiducial( 1e-2 );
int    photon::scat_max ( 20   );

photon::uint photon::n_repeat( 20 );

std::map<double, photon::uint> photon::bin_map;
std::vector<double>        photon::eta_upper;

////////////////////////////////////////////////////////////
// Initializer

photon::photon(  ) : exp_rand( 1. ), uni_rand( 0, 1 )
{
    prof = profile::get_instance(  );
    generator.seed( rand_seed::get_seed( "photon" ) );
    return;
}

photon::~photon(  )
{
    return;
}

void photon::init( input & args )
{
    float n_photon( 0. );
    args.find_key( "theta_bb", theta_bb, 2e-6 );
    args.find_key( "r_max"   , r_max,    1.   );
    args.find_key( "scat_max", scat_max, 20   );
    args.find_key( "n_photon", n_photon, 20   );
    
    int n_thread( 1 );
    args.find_key( "n_thread", n_thread, 1   );
    n_repeat = static_cast<uint>
	( n_photon / abs( n_thread ) + 1 );

    args.find_key( "d_tau", d_tau_fiducial, 1e-2 );

    int n_bin( 0 );
    double eta_min( 0. ), eta_max( 0. );
    args.find_key( "n_photon_bin",   n_bin,   100  );
    args.find_key( "eta_photon_min", eta_min, 1e-7 );
    args.find_key( "eta_photon_max", eta_max, 1e2  );
    const double e0 = log( eta_min );
    const double e1 = log( eta_max );
    const double de = ( e1 - e0 ) / ( n_bin - 1 );
    for( int i = 0; i < n_bin; ++ i )
    {
	const double eta = exp( e0 + de * i );
	eta_upper.push_back( eta );
	bin_map.insert( std::make_pair( eta, i ) );
    }
    return;
}

////////////////////////////////////////////////////////////
// Location and momentum

void photon::init_loc(  )
{
    x = { };
    continue_walking = true;
    return;
}

void photon::init_mom(  )
{
    auto p_planck  = rand_planck::get_instance(  );
    const double e = theta_bb * p_planck->get_rand(  );

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
    static const double tiny( 1e-4 );
    double tau_gone( 0. );

    double rho_ratio = prof->rho_ratio( x );
    double d_tau     = d_tau_fiducial;

    const double & eta = p[ 0 ]; // photon energy

    // Klein-Nishina factor from Mathematica... messy...
    double kn_factor( 0. );
    if( eta < tiny )
	kn_factor = 1.;
    else    
	kn_factor = ((2*eta*(2 + eta*(1 + eta)*(8 + eta)))
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
	    d_tau = fabs( tau - tau_gone ) * ( 1 + tiny );
	// ( 1 + tiny ): the loop won't stuck
	if( this->radius_c(  ) > r_max )
	{
	    continue_walking = false;
	    break;
	}
    }
    return;
}

void photon::locate_bin( const double & eta )
{
    auto p = bin_map.lower_bound( eta );
    if( p == bin_map.end(  ) )
	return;
    ++ res[ p->second ];
    return;
}

void photon::iterate_photon(  )
{
#pragma omp critical
    std::cout << "There are " << n_repeat << " photons "
	      << "being simulated on thread "
	      << omp_get_thread_num(  ) << std::endl;

    res.resize( eta_upper.size(  ), 0 );
    
    for( uint i = 0; i < n_repeat; ++ i )
    {
	reset(  );
	for( int j = 0; j < scat_max; ++ j )
	{
	    const double tau = exp_rand( generator );
	    step_walk( tau );
	    if( ! continue_walking )
		break;
	
	    elec.scatter_ph( this->p, prof->theta( x ) );
	}
	locate_bin( p[ 0 ] );
    }
    
    return;
}

const std::vector<photon::uint> & photon::get_res(  )
{
    return res;
}

const std::vector<double> & photon::get_eta_upper(  )
{
    return eta_upper;
}
