#pragma once
#include "gmp.h"
#include "gmpxx.h"
#include <vector>

enum EcmAlgorithm{
    Montgomery1,
    Montgomery2,
    Weierstrass1,
    Weierstrass2
};


void Ecm(mpz_t output, mpz_t to_factor, int thread_count, EcmAlgorithm algorithm);
