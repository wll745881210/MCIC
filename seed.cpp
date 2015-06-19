#include "seed.h"

#include <random>
#include <chrono>
#include <omp.h>

const std::unordered_map<std::string, unsigned>
rand_seed::kind_operator =
{
    { "electron", 51151230 },
    { "photon",   13742137 },
    { "gamma",    47126899 },
    { "knscat",   62478482 },
    { "planck",   60453484 },
    { "default",  78772348 }
};

unsigned rand_seed::get_seed( const std::string kind )
{
    auto p = rand_seed::kind_operator.find( kind );
    const unsigned op = p->second;
    const unsigned th = omp_get_thread_num(  );

    typedef std::mt19937::result_type  seed_type;
    typename std::chrono::system_clock seed_clock;

    auto ch = static_cast<seed_type>
	( seed_clock.now( ).time_since_epoch( ).count( ) );

    ch ^= op;
    ch ^= th;
    return ch;
}

