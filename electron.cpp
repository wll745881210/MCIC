#include "electron.h"

#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
// Static variables and functors

std::default_random_engine electron::generator;

////////////////////////////////////////////////////////////
// Constructor and destructor

electron::electron(  ) : uni_rand( 0, 1 )
{
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
  const double & mu  , std::array<double, 4> & vec )
{
    std::array< std::array< double, 4 >, 4 > Lambda = {};

    const double mu2 = mu * mu;
    const double cmu = sqrt( 1 - mu2 );
    
    Lambda[ 0 ][ 0 ] = gamma;
    Lambda[ 0 ][ 1 ] = -beta * gamma * mu;
    Lambda[ 0 ][ 2 ] = -beta * gamma * cmu;
    Lambda[ 1 ][ 0 ] = Lambda[ 0 ][ 1 ];
    Lambda[ 1 ][ 1 ] = 1 + ( gamma - 1 ) * mu2;
    Lambda[ 1 ][ 2 ] = ( gamma - 1 ) * mu * cmu;
    Lambda[ 2 ][ 0 ] = Lambda[ 0 ][ 2 ];
    Lambda[ 2 ][ 1 ] = Lambda[ 1 ][ 2 ];
    Lambda[ 2 ][ 2 ] = 1 + ( gamma - 1 ) * ( 1 - mu2 );
    Lambda[ 3 ][ 3 ] = 1.;

    std::array< double, 4 > res = {};
    for( unsigned i = 0; i < 4; ++ i )
	for( unsigned j = 0; j < 4; ++ j )
	    res[ i ] += Lambda[ i ][ j ] * vec[ j ];
    vec.swap( res );
    return;
}

double electron::get_mu_2d( const std::array< double, 4 > v )
{
    double norm( 0. );
    norm = sqrt( v[ 1 ] * v[ 1 ] + v[ 2 ] * v[ 2 ] );
    return v[ 1 ] / norm;
}

void electron::rotate_back_mu_2d
( std::array< double, 4 > & v, const double & mu )
{
    const double cmu = sqrt( 1. - mu * mu );
    const double x1  = v[ 1 ];
    const double x2  = v[ 2 ];
    v[ 1 ] = x1 * mu - x2 * cmu;
    v[ 2 ] = x2 * mu + x1 * cmu;
    return;
}

void electron::scatter_ph
( std::array< double, 4 > & p_ph, const double & theta )
{
    auto     ptr_gamma = rand_gamma::get_instance(  );
    const double gamma = ptr_gamma->get_rand_gamma( theta );
    const double beta  = sqrt( 1. - 1. / pow( gamma, 2 ) );
    const double mu_e  = uni_rand( generator ) * 2 - 1;
    
    lorentz_trans( beta, gamma, mu_e, p_ph );
    const double mu_ph = get_mu_2d( p_ph );

    auto ptr_knscat    = rand_knscat::get_instance(  );
    const double & eta = p_ph[ 0 ];
    const double mu_s  = ptr_knscat->get_rand_mu( eta );

    const double f_compton = 1 / ( 1 + eta * ( 1 - mu_s ) );
    p_ph[ 0 ] *= f_compton;
    p_ph[ 1 ]  = p_ph[ 0 ] * mu_s;
    p_ph[ 2 ]  = p_ph[ 0 ] * sqrt( 1. - mu_s * mu_s );
    rotate_back_mu_2d( p_ph, mu_ph );

    lorentz_trans( -beta, gamma, mu_e, p_ph );
    
    return;
}
