#ifndef SEED_H_
#define SEED_H_

#include <unordered_map>
#include <string>

namespace rand_seed
{

extern const std::unordered_map< std::string, unsigned >
kind_operator;

unsigned get_seed( const std::string kind = "default" );

}

#endif
