#include "input.h"
#include "driver.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
// #include <omp.h>

////////////////////////////////////////////////////////////

int main( int argn, char * argv[  ] )
{
    try
    {
	std::string par_file_name;
	if( argn > 2 )
	    throw "Incorrect input parameter.\n"
		"See README for usage.";
	else if( argn < 2 )
	    par_file_name = "par.txt";
	else
	    par_file_name = argv[ 1 ];

	input args( par_file_name );
	driver( args );
    }
    catch( const char * err )
    {
	std::cerr << "\nError: " << err << std::endl;
	return -1;
    }

    return 0;
}

// int main(  )
// {
//     try
//     {
// 	auto test = rand_planck::get_instance(  );
// 	test->set_resolution( 1e2 );
// 	test->integrate(  );

// 	// std::vector<double> res;
// 	// for( int i = 0; i < 1e6; ++ i )
// 	//     res.push_back( test->get_rand_planck(  ) );

// 	// std::ofstream fout( "test.dat" );
// 	// for( unsigned i = 0; i < res.size(  ); ++ i )
// 	//     fout << res[ i ] << '\n';
// 	// rand_planck::del_instance(  );

// 	photon ph_test;
// 	ph_test.step_walk( 0.1 );
//     }
//     catch( const char * err )
//     {
// 	std::cerr << "\nError: " << err << std::endl;
// 	return -1;
//     }

//     return 0;
// }
