#pragma once
#include "gmp.h"
#include <vector>
#include "FactorisationStats.h"

void EcmCuda(mpz_t output, mpz_t to_factor, int thread_count, const std::vector<long>& primes);
