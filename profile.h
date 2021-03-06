#ifndef PROFILE_H_
#define PROFILE_H_

////////////////////////////////////////////////////////////
// Density and temperature profiles.
// Density values are normalized so that the value of
// density has unitary value at the center.
// Temperature values are normalized onto "theta", viz.
// theta := k_B T / ( m_e c^2 ).
// Originally rho_ratio and theta are very simple.
// If you want to use interpolation functions for them,
// go for it by modifying this file and "profile.cpp".

#include "input.h"

#include <array>

class profile
{
    ////////// Initializers //////////
private:
    static profile * singleton;
     profile(  );
    ~profile(  );
public:				// Function
    static profile * get_instance(  );
    static void      del_instance(  );
    void init( input & args );

    ////////// Profile values //////////
private:			// Data
    double r_fid;
    double n_fid,     n_pow_inner, n_pow_outer;
    double theta_fid, t_pow_inner, t_pow_outer;
private:			// Function
    double radius( const std::array< double, 3 > & x );
public:				// Function
    double n_e  ( const std::array< double, 3 > & x );
    double theta( const std::array< double, 3 > & x );
};

#endif
