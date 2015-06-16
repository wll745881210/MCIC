#include "electron.h"

// 

////////////////////////////////////////////////////////////
// Static variables and functors

std::default_random_engine electron::generator;
double electron::theta_te( 1.68637e-10e );

////////////////////////////////////////////////////////////
// Constructor and destructor

elctron::electron(  ) : px( 0. ), py( 0. ), pz( 0. ),
			mu_rand( -1, 1 ), fe_rand( 0, 1 ),
			phi_rand( 0, 6.2831853 )
{
    class rel_e_gamma : public my_random
    {
	double pdf( const double & x ) const
	{
	    return 1. / ( theta * k2 ) * x * ;
	}
	
    public:
	double theta, k2;
    };

    return;
}

electron::~electron(  )
{
    return;
}

////////////////////////////////////////////////////////////
// Generate a random electron momentum based on temperature

void electron::get_momentum( const double & t_e )
{
    const double Theta = theta_te * t_e; // Big Theta -> T_e
    const double k2    = cyl_bessel_k( 2, 1. / Theta ) ;


    
    
    
    
}

