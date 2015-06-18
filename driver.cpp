#include "driver.h"
#include "photon.h"
#include "rand_gamma.h"
#include "rand_knscat.h"
#include "rand_planck.h"

#include <vector>
#include <string>
#include <fstream>
#include <omp.h>

////////////////////////////////////////////////////////////
// Driver function

void driver( input & args )
{
    std::cout << "Initializing... " << std::flush;
    args.read(  );
    int n_thread( 1 );
    args.find_key( "n_thread", n_thread, 1 );
    if( n_thread > 1 )
	omp_set_num_threads( n_thread );
    else
	n_thread = omp_get_num_threads(  );
    auto p_gamma  = rand_gamma ::get_instance(  );
    auto p_knscat = rand_knscat::get_instance(  );
    auto p_planck = rand_planck::get_instance(  );
    p_gamma ->init( args );
    p_knscat->init( args );
    p_planck->init( args );
    photon::init( args );
    std::vector<photon> photon_arr( n_thread );
    std::cout << "Done." << std::endl;

    std::cout << "Doing MC simulation... " << std::flush;
#pragma omp parallel for
    for( int i = 0; i < n_thread; ++ i )
	photon_arr[ i ].iterate_photon(  );
    std::cout << "Done." << std::endl;
    
    std::cout << "Dumping data... " << std::flush;
    std::string output_path;
    args.find_key( "output_path", output_path, "mc.dat" );
    std::ofstream( output_path.c_str(  ) ) fout;    
    for( int i = 0; i < n_thread; ++ i )
    {
	auto & res = photon_arr[ i ].get_res(  );
	for( unsigned j = 0; j < res.size(  ); ++ j )
	    fout << res[ j ] << '\n';
    }
    std::cout << "Done." << std::endl;
    return;
}
