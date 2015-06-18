#ifndef RAND_GAMMA_H_
#define RAND_GAMMA_H_

#include "rand_base.h"
#include "input.h"

#include <random>
#include <map>

////////////////////////////////////////////////////////////
// Here the distribution function is integrated to give the
// cumulative distribution, by which a uniform distribution
// is mapped onto desired distribution. The best way of
// doing this is to get the distribution function for
// x := exp( - ( gamma - 1 ) / theta ).
// In turn, the parameter t for class rand_base is
// t := log( theta ).
// The normalized PDF for x would then be
// f_x = exp( - 1 / theta ) / K_2( 1 / theta )
//     * ( 1 - theta * ln( x ) )
//     * ( theta**2 ln**2( x ) - 2 * theta * ln( x ) )^(1/2)
// and P = \int_0^x f_x' dx'.
// Then, integration is conducted from ( P = 1, x = 1 ) to
// ( P = 0, x = 0 ).
// BTW, when theta < 1e-2, we shall reduce our PDF to the NR
// Maxwellian case, but let's keep it as it is at this
// moment.

class rand_gamma : public rand_base
{
    ////////// Initializers //////////
private:
    static rand_gamma * singleton;
public:
    static rand_gamma * get_instance(  );
    static void         del_instance(  );
    void init( input & args );

    ////////// PDF related //////////
private:			// Data
    double theta, ln_theta, norm_theta;
private:			// Function
    double pdf( const double & x );

    ////////// Integration related //////////
private:
    std::map<double, double> *
    intg_single_t( const double & t );
    void prepare_intg(  );
    
public:				// Function
    double get_rand_gamma( const double & theta );
};

#endif
