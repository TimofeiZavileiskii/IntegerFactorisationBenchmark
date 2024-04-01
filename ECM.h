#pragma once
#include "gmp.h"
#include "gmpxx.h"
#include <vector>

void EcmTest(mpz_t& output, mpz_t& to_factor);

void Ecm(mpz_t& output, mpz_t& to_factor, int thread_count);