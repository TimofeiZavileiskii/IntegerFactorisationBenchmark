#pragma once
#include "gmp.h"
#include "gmpxx.h"
#include <set> 

void PollardsRho(mpz_t& output, mpz_t& to_factor, int thread_count);