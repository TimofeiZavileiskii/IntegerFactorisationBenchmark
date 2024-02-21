#pragma once
#include "gmp.h"
#include "gmpxx.h"
#include <vector>


void EcmJob(mpz_t& output, mpz_t& to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed);


void Ecm(mpz_t& output, mpz_t& to_factor, int thread_count);