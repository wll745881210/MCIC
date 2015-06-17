#include "rand_gamma.h"

#include <iostream>
#include <boost/math/special_functions/bessel.hpp>

double rand_gamma::pdf( const double & x )
{
    const double lnx = log( x );
    return norm_theta * ( 1 - theta * lnx ) *
	sqrt( pow( theta * lnx, 2 ) - 2 * theta * lnx );
}
