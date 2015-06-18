#ifndef RAND_PLANCK_H_
#define RAND_PLANCK_H_

#include "rand_base.h"
#include "input.h"

////////////////////////////////////////////////////////////
// Planck distribution generator. Basically only using the
// one-dimensional feature of rand_base, since it is easy to
// reduce Planck distribution into a 1-D distribution.
// Now the random variable is:
// x := h \nu / ( k_B T )

class rand_planck : public rand_base
{
    ////////// Initializers //////////
private:
    static rand_planck * singleton;
public:
    static rand_planck * get_instance(  );
    static void          del_instance(  );
    void init( input & args );

    ////////// PDF related //////////
private:			// Function
    double pdf( const double & x );

    ////////// Integration related //////////
private:			// Data
    std::map<double, double> * res;
    static const double epsilon; // Tiny number
private:			// Function
    std::map<double, double> *
    intg_single_t( const double & t );
    void prepare_intg(  );    
public:				// Function
    double get_rand(  );    
};




#endif
