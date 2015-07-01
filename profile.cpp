#include "profile.h"

#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
// Static variables

profile * profile::singleton( nullptr );

////////////////////////////////////////////////////////////
// Initializers

profile::profile(  )
{
    return;
}

profile::~profile(  )
{
    return;
}

profile * profile::get_instance(  )
{
    if( singleton == nullptr )
	singleton = new profile;
    return singleton;
}

void profile::del_instance(  )
{
    if( singleton != nullptr )
	delete singleton;
    return;
}

void profile::init( input & args )
{
    // All in CGS...
    args.find_key( "r_fid",       r_fid,       1. );
    args.find_key( "n_fid",       n_fid,       1. );
    args.find_key( "theta_fid",   theta_fid,   1. );
    args.find_key( "n_pow_inner", n_pow_inner, 0. );
    args.find_key( "n_pow_outer", n_pow_outer, 0. );
    args.find_key( "t_pow_inner", t_pow_inner, 0. );
    args.find_key( "t_pow_outer", t_pow_outer, 0. );

    return;
}

////////////////////////////////////////////////////////////
// Get the values.
// If you want to use it for further purposes, modify here.

double profile::radius( const std::array< double, 3 > & x )
{
    double r2( 0. );
    for( int i = 0; i < 3; ++ i )
	r2 += x[ i ] * x[ i ];
    return sqrt( r2 );
}

double profile::n_e
( const std::array< double, 3 > & x )
{
    const double r = radius( x );
    const double idx
	= ( r > r_fid ? n_pow_outer : n_pow_inner );
    return n_fid * pow( r / r_fid, idx );
}

double profile::theta( const std::array< double, 3 > & x )
{
    const double r = radius( x );
    const double idx
	= ( r > r_fid ? t_pow_outer : t_pow_inner );
    return theta_fid * pow( r / r_fid, idx );
}
