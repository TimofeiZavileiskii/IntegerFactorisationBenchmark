#pragma once
#include "gmp.h"
#include "gmpxx.h"
#include <vector>


void Ecm(mpz_t& output, mpz_t& to_factor, int thread_count);