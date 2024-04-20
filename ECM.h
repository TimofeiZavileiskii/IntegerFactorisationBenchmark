#pragma once
#include "gmp.h"
#include "gmpxx.h"
#include <vector>
#include "FactorisationStats.h"

enum EcmAlgorithm{
    Montgomery1,
    Montgomery2,
    Weierstrass1,
    Weierstrass2
};

enum BoundChoice{
    DEFAULT_BOUNDS,
    TABLE_BOUNDS,
    GPU_TABLE_BOUNDS
};


void ChooseBounds(long& bound, long& bound2, mpz_t to_factor, BoundChoice choice, int offset);


void Ecm(mpz_t output, mpz_t to_factor, int thread_count, EcmAlgorithm algorithm, const std::vector<long>& primes, StatsEcm* stats);

