#include "input.h"
#include "photon.h"
#include "rand_gamma.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
// #include <omp.h>

////////////////////////////////////////////////////////////

// void driver( input & args )
// {
//     int n_thread( 1 );
//     args.find_key( "n_thread", n_thread, 1 );
//     if( n_thread > 1 )
// 	omp_set_num_threads( n_thread );

//     std::vector<photon> ph_list( n_thread );

    
//     double tau_max( 1.  );
//     int    itr_max( 20  );
//     args.find_key( "tau_max", tau_max, 1. );
//     args.find_key( "itr_max", itr_max, 20 );
    

//     double n_ph( 100 );
//     args.find_key( "num_photon", n_ph, 100 );
//     int n_walk( static_cast<int>( n_ph / n_thread ) );

//     photon::set_max_tau( tau_max );
//     photon::set_max_itr( itr_max );
//     photon::set_n_walk ( n_walk  );

//     std::cout << "Walking photons... " << std::flush;
// #pragma omp parallel for
//     for( unsigned i = 0; i < ph_list.size(  ); ++ i )
// 	ph_list[ i ].walk_photon(  );
//     std::cout << "Done." << std::endl;

//     std::vector<double> ph_bin( itr_max );
//     for( int i = 0; i < itr_max; ++ i )
// 	for( int j = 0; j < n_thread; ++ j )
// 	{
// 	    auto res = ph_list[ j ].result(  );
// 	    ph_bin[ i ] += res[ i ] / n_thread;
// 	}

//     std::string res_path;
//     args.find_key( "res_file", res_path, "esc_prob.txt" );
//     std::ofstream fout( res_path.c_str(  ) );
//     for( unsigned i = 0; i < ph_bin.size(  ); ++ i )
// 	fout << i << '\t' << ph_bin[ i ]<< '\n';

//     fout.close(  );
    
//     return;
// }

// int main( int argn, char * argv[  ] )
// {
//     try
//     {
// 	std::string par_file_name;
// 	if( argn > 2 )
// 	    throw "Incorrect input parameter.\n"
// 		"See README for usage.";
// 	else if( argn < 2 )
// 	    par_file_name = "par.txt";
// 	else
// 	    par_file_name = argv[ 1 ];

// 	input args( par_file_name );
// 	args.read(  );

// 	driver( args );
//     }
//     catch( const char * err )
//     {
// 	std::cerr << "\nError: " << err << std::endl;
// 	return -1;
//     }

//     return 0;
// }

int main(  )
{
    try
    {
	auto test = rand_gamma::get_instance(  );
	test->set_intg_pts( 1000 );
	test->set_theta( 1e-3, 1e0, 5 );
	test->integrate(  );

	std::vector<double> res;
	for( int i = 0; i < 1e6; ++ i )
	    res.push_back( test->get_rand_gamma( 1e0 ) );

	std::ofstream fout( "test.txt" );
	for( unsigned i = 0; i < res.size(  ); ++ i )
	    fout << res[ i ] << '\n';
	rand_gamma::del_instance(  );
    }
    catch( const char * err )
    {
	std::cerr << "\nError: " << err << std::endl;
	return -1;
    }

    return 0;
}
