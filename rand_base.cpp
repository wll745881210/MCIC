#include "rand_base.h"

#include <cmath>

////////////////////////////////////////////////////////////
// Static variables
rand_base * rand_base::singleton( nullptr );

////////////////////////////////////////////////////////////
// Initializers

rand_gamma::rand_gamma(  ) : t_rand( 0, 1 ), d_x( -1 )
{
    return;
}

rand_gamma::~rand_gamma(  )
{
    for( auto p = theta_map.begin(  );
	 p != theta_map.end(  ); ++ p )
	delete ( p->second );

    return;
}

rand_gamma * rand_gamma::get_instance(  )
{
    if( singleton == nullptr )
	singleton = new rand_gamma;
    return singleton;
}

void rand_gamma::del_instance(  )
{
    if( singleton != nullptr )
	delete singleton;
    return;
}

void rand_gamma::set_intg_pts( const int & n_x )
{
    this->n_x = n_x;
    this->d_x = 1. / n_x;	// Final point is ( 0, 0 )
    return;
}

void rand_gamma::set_theta
( const double & theta0, const double & theta1,
  const int    & n_theta )
{
    this->ln_theta0  = log( theta0 );
    this->ln_theta1  = log( theta1 );
    this->n_theta    = n_theta;
    this->d_ln_theta = ( ln_theta1 - ln_theta0 )
	             / ( n_theta - 1 );
    return;
}


////////////////////////////////////////////////////////////
// Distribution function

    
////////////////////////////////////////////////////////////
// Runge-Kutta integration--specifically designed

void rand_gamma::rk4( double & p, const double & x )
{
    const double dp0 = pdf( x );
    const double dp1 = pdf( x - d_x / 2 );
    const double dp3 = pdf( x - d_x     );
    p -= d_x * ( dp0 + 4 * dp1 + dp3 ) / 6.;
    return;
}

void rand_gamma::intg_single_theta(  )
{
    double p( 1 );
    auto res = new std::map<double, double>;
    for( int i = 0; i < n_x; ++ i ) // Better than while
    {
	const double x = 1. - d_x * i;
	res->insert( std::make_pair( p, x ) );
	rk4( p, x );
	if( p < 0 )
	    break;
    }
    res->insert( std::make_pair( 0., 0. ) );

    theta_map.insert( std::make_pair( ln_theta, res ) );

    return;
}

void rand_gamma::integrate(  )
{
    for( int i = 0; i < n_theta; ++ i )
    {
	ln_theta = ln_theta0 + d_ln_theta * i;
	theta    = exp( ln_theta );
	norm_theta = exp( - 1. / theta )
	    / boost::math::cyl_bessel_k( 2, 1. / theta );
    	intg_single_theta(  );
    }
    return;
}

////////////////////////////////////////////////////////////
// Interpolation

double rand_gamma::interp_x
( const double & y, std::map<double, double> * p )
{
    auto q = p->lower_bound( y );
    if( q == p->end(  ) )
	return 1.;
    else if( q == p->begin(  ) )
	return 0.;
    else
    {
	const double y1 = q->first;
	const double x1 = q->second;
	-- q;
	const double y0 = q->first;
	const double x0 = q->second;
	return x0 + ( x1 - x0 ) / ( y1 - y0 ) * ( y - y0 );
    }
}

double rand_gamma::cdf
( const double & y, const double & lnt )
{
    auto p = theta_map.lower_bound( lnt );
    if( p == theta_map.end(  ) )
	throw "Electron temperature too high";
    else if( p == theta_map.begin(  ) )
	return interp_x( y, p->second );
    else
    {
	const double lnt1 = p->first;
	const double x1   = interp_x( y, p->second );
	-- p;
	const double lnt0 = p->first;
	const double x0   = interp_x( y, p->second );
	if( x0 > 1 || x1 > 1 )
	    throw "Interp error";

	return x0 + ( x1 - x0 ) *
	    ( lnt - lnt0 ) / ( lnt1 - lnt0 );
    }
}

double rand_gamma::get_rand_gamma( const double & theta )
{
    const double t = t_rand( generator );
    const double x = cdf( t, log( theta ) );
    return 1. - theta * log( x );
}

    
