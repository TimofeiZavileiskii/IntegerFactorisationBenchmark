#pragma once
#include "gmp.h"


void ThreadJobRho(mpz_t& output, mpz_t& to_factor, volatile bool& factored, int seed);


void PollardsRho(mpz_t& output, mpz_t& to_factor, int thread_count);