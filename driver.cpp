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
    auto p_prof   = profile    ::get_instance(  );
    p_gamma ->init( args );
    p_knscat->init( args );
    p_planck->init( args );
    p_prof  ->init( args );
    
    photon::init( args );
    std::vector<photon> photon_arr( n_thread );
    std::cout << "Done.\n\n" << std::endl;

    std::cout << "Doing MC simulation... " << std::endl;
#pragma omp parallel for
    for( int i = 0; i < n_thread; ++ i )
	photon_arr[ i ].iterate_photon(  );
    std::cout << "Done.\n\n" << std::endl;

    std::cout << "Dumping data... " << std::flush;
    
    std::string output_path;
    args.find_key( "output_path", output_path, "mc.dat" );
    std::ofstream fout_e( output_path.c_str(  ) );
    auto output_path_scat = output_path + "_nscat";
    std::ofstream fout_s( output_path_scat.c_str(  ) );

    int scat_max( 0 );
    args.find_key( "scat_max", scat_max, 0 );
    auto & e_upper = photon::get_eta_upper(  );
    std::vector<double> res_s( scat_max + 1     );
    std::vector<double> res_e( e_upper.size(  ) );
	
    for( int i = 0; i < n_thread; ++ i )
    {
	auto & res_i_e = photon_arr[ i ].get_res_e   (  );
	auto & res_i_s = photon_arr[ i ].get_res_scat(  );
	for( unsigned j = 0; j < res_e.size(  ); ++ j )
	    res_e[ j ] += res_i_e[ j ];
	for( unsigned j = 0; j < res_s.size(  ); ++ j )
	    res_s[ j ] += res_i_s[ j ];
    }
    for( unsigned i = 0; i < res_e.size(  ); ++ i )
	fout_e << e_upper[ i ] << '\t'
	       << res_e[ i ] << '\n';
    for( unsigned i = 0; i < res_s.size(  ); ++ i )
	fout_s << i << '\t' << res_s[ i ] << '\n';
    
    std::cout << "Done.\n" << std::endl;
    return;
}
