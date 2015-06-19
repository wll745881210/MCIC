#include "input.h"
#include "driver.h"
#include "electron.h"

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <fstream>

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
// 	// auto test = rand_planck::get_instance(  );
// 	// test->set_resolution( 1e2 );
// 	// test->integrate(  );

// 	// std::vector<double> res;
// 	// for( int i = 0; i < 1e6; ++ i )
// 	//     res.push_back( test->get_rand_planck(  ) );

// 	// std::ofstream fout( "test.dat" );
// 	// for( unsigned i = 0; i < res.size(  ); ++ i )
// 	//     fout << res[ i ] << '\n';
// 	// rand_planck::del_instance(  );

// electron test;
// std::array<double, 4> a = { 3.3, 2.2, 1.1, 0.5 };
// test.scatter_ph( a, 1 );

// }

//     catch( const char * err )
//     {
// 	std::cerr << "\nError: " << err << std::endl;
// 	return -1;
//     }

//     return 0;
// }
