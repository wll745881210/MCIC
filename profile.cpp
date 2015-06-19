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

////////////////////////////////////////////////////////////
// Get the values.
// If you want to 

double profile::radius( const std::array< double, 3 > & x )
{
    double r2( 0. );
    for( int i = 0; i < 3; ++ i )
	r2 += x[ i ] * x[ i ];
    return sqrt( r2 );
}

double profile::rho_ratio
( const std::array< double, 3 > & x )
{
    return 1.;
}

double profile::theta( const std::array< double, 3 > & x )
{
    return 0.5;
}
