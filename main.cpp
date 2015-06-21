#include "input.h"
#include "driver.h"
#include "electron.h"
#include "rand_gamma.h"

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
