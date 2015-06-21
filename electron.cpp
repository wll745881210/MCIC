#include "electron.h"

#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
// Constructor and destructor

electron::electron(  ) : uni_rand( 0, 1 )
{
    generator.seed( rand_seed::get_seed( "electron" ) );
    return;
}

electron::~electron(  )
{
    return;
}

////////////////////////////////////////////////////////////
// Scatter a photon with back-and-forth Lorentz transform.

void electron::lorentz_trans
( const double & beta, const double & gamma,
  const double & mu  , const double & phi,
  std::array<double, 4> & vec )
{
    std::array< std::array< double, 4 >, 4 > Lambda = {};

    const double cmu = sqrt( 1 - mu * mu ); // sin( theta )
    const double b_x = beta * cmu * cos( phi );
    const double b_y = beta * cmu * sin( phi );
    const double b_z = beta * mu;
    const double b2  = beta * beta;
    const double gm1 = gamma - 1.;
    // This mu is NOT scattering mu!
    
    Lambda[ 0 ][ 0 ] = gamma;
    Lambda[ 0 ][ 1 ] = -gamma * b_x;
    Lambda[ 0 ][ 2 ] = -gamma * b_y;
    Lambda[ 0 ][ 3 ] = -gamma * b_z;
    Lambda[ 1 ][ 0 ] = Lambda[ 0 ][ 1 ];
    Lambda[ 1 ][ 1 ] = 1 + gm1 * b_x * b_x / b2;
    Lambda[ 1 ][ 2 ] = gm1 * b_x * b_y / b2;
    Lambda[ 1 ][ 3 ] = gm1 * b_x * b_z / b2;
    Lambda[ 2 ][ 0 ] = Lambda[ 0 ][ 2 ];
    Lambda[ 2 ][ 1 ] = Lambda[ 1 ][ 2 ];
    Lambda[ 2 ][ 2 ] = 1 + gm1 * b_y * b_y / b2;
    Lambda[ 2 ][ 3 ] = gm1 * b_y * b_z / b2;
    Lambda[ 3 ][ 0 ] = Lambda[ 0 ][ 3 ];
    Lambda[ 3 ][ 1 ] = Lambda[ 1 ][ 3 ];
    Lambda[ 3 ][ 2 ] = Lambda[ 2 ][ 3 ];
    Lambda[ 3 ][ 3 ] = 1 + gm1 * b_z * b_z / b2;

    std::array< double, 4 > res = {};
    for( unsigned i = 0; i < 4; ++ i )
	for( unsigned j = 0; j < 4; ++ j )
	    res[ i ] += Lambda[ i ][ j ] * vec[ j ];
    vec.swap( res );
    return;
}

void electron::rotate_back_ph
( std::array< double, 4 > & p_old,
  std::array< double, 4 > & p_new )
{
    const double p_norm =
	sqrt( pow( p_old[ 1 ], 2 ) + pow( p_old[ 2 ], 2 )
	    + pow( p_old[ 3 ], 2 ) );
    const double mu   = p_old[ 3 ] / p_norm;
    const double cmu  = sqrt( 1. - mu * mu );
    const double phi  = atan2( p_old[ 2 ], p_old[ 1 ] );
    const double cphi = cos( phi );
    const double sphi = sin( phi );
    
    std::array< std::array< double, 3 >, 3 > R = {};
    R[ 0 ][ 0 ] = cmu * cphi;
    R[ 0 ][ 1 ] = -sphi;
    R[ 0 ][ 2 ] = -cphi * mu;
    R[ 1 ][ 0 ] = sphi * cmu;
    R[ 1 ][ 1 ] = cphi;
    R[ 1 ][ 2 ] = -sphi * mu;
    R[ 2 ][ 0 ] = mu;
    R[ 2 ][ 1 ] = 0.;
    R[ 2 ][ 2 ] = cmu;
    
    for( int i = 0; i < 3; ++ i )
    {
	p_old[ i + 1 ] = 0.;
	for( int j = 0; j < 3; ++ j )
	    p_old[ i + 1 ] += R[ i ][ j ] * p_new[ j + 1 ];
    }
    
    return;
}

void electron::scatter_ph
( std::array< double, 4 > & p_ph, const double & theta )
{
    auto     ptr_gamma = rand_gamma::get_instance(  );
    const double gamma = ptr_gamma->get_rand_gamma( theta );
    const double beta  = sqrt( 1. - 1. / pow( gamma, 2 ) );
    const double mu_e  = uni_rand( generator ) * 2 - 1;
    const double phi_e = uni_rand( generator ) * 6.283185;
    
    lorentz_trans( beta, gamma, mu_e, phi_e, p_ph );

    auto    ptr_knscat = rand_knscat::get_instance(  );
    const double & eta = p_ph[ 0 ];
    const double mu_s  = ptr_knscat->get_rand_mu( eta );
    const double cmu_s = sqrt( 1. - mu_s * mu_s );
    const double phi_s = uni_rand( generator ) * 6.283185;

    const double f_compton = 1 / ( 1 + eta * ( 1 - mu_s ) );
    const double eta_new = p_ph[ 0 ] * f_compton;
    std::array<double, 4> p_ph_new;
    // Note the strange coordinate system!
    // Now the scattering axis is x-axis!	
    p_ph_new[ 0 ] = eta_new;
    p_ph_new[ 1 ] = eta_new * mu_s;
    p_ph_new[ 2 ] = eta_new * cmu_s * cos( phi_s );
    p_ph_new[ 3 ] = eta_new * cmu_s * sin( phi_s );
    
    p_ph[ 0 ] = eta_new;
    rotate_back_ph( p_ph, p_ph_new );

    lorentz_trans( -beta, gamma, mu_e, phi_e, p_ph );
        
    return;
}
