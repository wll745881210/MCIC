#include "rand_base.h"

#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////
// Initializers

rand_base::rand_base(  ) : cdf_rand( 0, 1 )
{
    return;
}

rand_base::~rand_base(  )
{
    for( auto & p : t_map )
	delete p.second;
    return;
}

////////////////////////////////////////////////////////////
// Distribution function is purely virtual for the base.
    
////////////////////////////////////////////////////////////
// Runge-Kutta integration--specifically designed

void rand_base::rk4
( double & c, const double & x, const double & dx )
{
    const double dc0 = pdf( x );
    const double dc1 = pdf( x + dx / 2 );
    const double dc3 = pdf( x + dx     );
    c += dx * ( dc0 + 4 * dc1 + dc3 ) / 6.;
    return;
}

//////////////////////////////
// This function is virtual--we shall do some
// "pre-/post-processings" when a derived random method is
// defined.
// Especially, we MUST initialize c_current.

std::map<double, double> *
rand_base::intg_single_t( const double & t )
{
    auto res = new std::map<double, double>;
    for( unsigned i = 0; i < x_vec.size(  ) - 1; ++ i )
    {
	const double & x = x_vec[ i ];
	
	const double dx  = x_vec[ i + 1 ] - x;
	res->insert( std::make_pair( c_current, x ) );
	rk4( c_current, x, dx );
    }

    return res;
}

void rand_base::prepare_intg(  )
{
    return;
}

void rand_base::integrate(  )
{
    prepare_intg(  );
    for( auto & t : t_vec )
    {
	auto p = intg_single_t( t );
    	t_map.insert( std::make_pair( t, p ) );
    }
    return;
}

////////////////////////////////////////////////////////////
// Interpolation

double rand_base::interp_single_t
( const double & c, std::map<double, double> * p )
{
    auto q = p->lower_bound( c );
    if( q == p->end(  ) )
    {
	-- q;
	return q->second;
    }
    else if( q == p->begin(  ) )
	return q->second;
    else
    {
	const double c1 = q->first;
	const double x1 = q->second;
	-- q;
	const double c0 = q->first;
	const double x0 = q->second;
	return x0 + ( x1 - x0 )  * ( c - c0 ) / ( c1 - c0 );
    }
}

double rand_base::interp_cdf
( const double & c, const double & t )
{
    auto p = t_map.lower_bound( t );
    if( p == t_map.end(  ) )
    {
	-- p;
	return interp_single_t( c, p->second );
    }
    else if( p == t_map.begin(  ) )
	return interp_single_t( c, p->second );
    else
    {
	const double t1 = p->first;
	const double x1 = interp_single_t( c, p->second );
	-- p;
	const double t0 = p->first;
	const double x0 = interp_single_t( c, p->second );

	return x0 + ( x1 - x0 ) * ( t - t0 ) / ( t1 - t0 );
    }
}

double rand_base::get_rand( const double & t )
{
    const double c = cdf_rand( generator );
    return interp_cdf( c, t );
}

    
