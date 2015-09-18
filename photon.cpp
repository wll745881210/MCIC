#include "photon.h"
#include "rand_planck.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <omp.h>

////////////////////////////////////////////////////////////
// Static variables
double photon::theta_bb      ( 2e-6 ); // 1eV in me c^2
double photon::r_max         ( 1. );
double photon::r_min         ( 1. );
double photon::r_eh          ( 1. );
double photon::r_disk_max    ( 0. );
double photon::r_disk_min    ( 0. );
double photon::d_tau_fiducial( 0. );
int    photon::scat_max      ( 20 );

bool   photon::is_simple_scat_model( false );
photon::uint photon::n_repeat( 20 );

std::map<double, photon::uint> photon::bin_map;
std::vector<double>        photon::eta_upper;

profile * photon::prof;

////////////////////////////////////////////////////////////
// Initializer

photon::photon(  ) : exp_rand( 1. ), uni_rand( 0, 1 )
{
    generator.seed( rand_seed::get_seed( "photon" ) );
    return;
}

photon::~photon(  )
{
    return;
}

void photon::init( input & args )
{
    prof = profile::get_instance(  );

    float n_photon( 0. );
    args.find_key( "theta_bb", theta_bb, 2e-6 );
    args.find_key( "scat_max", scat_max, 20   );
    args.find_key( "n_photon", n_photon, 20   );

    args.find_key( "is_simple_scat_model",
	            is_simple_scat_model, false );
    
    int n_thread( 1 );
    args.find_key( "n_thread", n_thread, 1   );
    n_repeat = static_cast<uint>
	( n_photon / abs( n_thread ) + 1 );

    args.find_key( "d_tau", d_tau_fiducial, 1e-2 );

    // Rebinning parameters
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

    // Disk as radiation source
    args.find_key( "r_disk_min", r_disk_min, 0. );
    args.find_key( "r_disk_max", r_disk_max, 0. );

    args.find_key( "r_max", r_max, 3.1e18       );
    args.find_key( "r_min", r_min, 1e-4 * r_max );
    args.find_key( "r_eh" , r_eh,  1e-4 * r_max );
    return;
}

////////////////////////////////////////////////////////////
// Location and momentum

void photon::init_loc(  )
{
    const double r_arg = pow( r_disk_min, -0.25 ) +
    	uni_rand( generator ) *
    	( pow( r_disk_max, -0.25 )
    	    - pow( r_disk_min, -0.25 ) );
    const double r   = pow( r_arg, -4 );
    const double phi = uni_rand( generator ) * 6.28319;
    x = { r * cos( phi ), r * sin( phi ), 0. };

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

    //////////////////////////////////////////////////
    // Calculate the cross section
    const double & eta = p[ 0 ]; // photon energy
    // Klein-Nishina factor from Mathematica... messy...
    double kn_factor( 0. );
    if( eta < tiny )
	kn_factor = 1.;
    else
	kn_factor = ((2*eta*(2 + eta*(1 + eta)*(8 + eta)))
	    /pow(1 + 2*eta,2) + (-2 + (-2 + eta)*eta)
	    *log(1 + 2*eta))/pow(eta,3) / ( 8. / 3. );
    const double sigma = kn_factor * 6.65246e-25;
    // Cross section in cm**2 ( = kn * sigma_T )
    //////////////////////////////////////////////////

    double d_tau( d_tau_fiducial ), dx( 0. );
    double tau_gone( 0. );
    double n_e0 = prof->n_e( x );
    
    while( tau_gone < tau )
    {
	const double r = radius_c(  );
	if( r > r_max || r < r_eh )
	{
	    continue_walking = false;
	    break;
	}
	else if( r < r_min )
	{
	    dx = r_min / 10.;
	    for( int i = 0; i < 3; ++ i )
		x[ i ] += dx * p[ i + 1 ] / eta;
	    continue;
	}
	
	// Predictor
        auto x_pre = x;
	dx = d_tau / ( n_e0 * sigma );
	for( int i = 0; i < 3; ++ i )
	    x_pre[ i ] += dx * p[ i + 1 ] / eta;
	const double n_e1 = prof->n_e( x_pre );
	// Corrector
	n_e0 = 0.5 * ( n_e0 + n_e1 );
	dx = d_tau / ( n_e0 * sigma );
	for( int i = 0; i < 3; ++ i )
	    x[ i ] += dx * p[ i + 1 ] / eta;
	// Save for the next step
	n_e0 = n_e1;
	
	tau_gone += d_tau;
	if( tau - tau_gone < d_tau )
	    d_tau = fabs( tau - tau_gone ) * ( 1 + tiny );
	// ( 1 + tiny ): the loop won't stuck
    }
    return;
}

void photon::locate_bin( const double & eta )
{
    auto p = bin_map.lower_bound( eta );
    if( p == bin_map.end(  ) )
	return;
    ++ res_e[ p->second ];
    return;
}

void photon::iterate_photon(  )
{
#pragma omp critical
    std::cout << "There are " << n_repeat << " photons "
	      << "being simulated on thread "
	      << omp_get_thread_num(  ) << std::endl;

    res_e.resize( eta_upper.size(  ), 0 );
    res_scat.resize( scat_max + 1,    0 );
    
    for( uint i = 0; i < n_repeat; ++ i )
    {
	reset(  );
	int j( 0 );
	for( j = 0; j < scat_max; ++ j )
	{
	    const double tau = exp_rand( generator );
	    step_walk( tau );
	    if( ! continue_walking )
	    	break;

	    if( is_simple_scat_model )
		elec.scatter_ph_simple
		    ( this->p, prof->theta( x ) );
	    else
		elec.scatter_ph
		    ( this->p, prof->theta( x ) );
	}
	++ res_scat[ j ];
	locate_bin( p[ 0 ] );
    }
    
    return;
}

const std::vector<photon::uint> & photon::get_res_e(  )
{
    return res_e;
}
const std::vector<photon::uint> & photon::get_res_scat(  )
{
    return res_scat;
}

const std::vector<double> & photon::get_eta_upper(  )
{
    return eta_upper;
}
