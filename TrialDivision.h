#pragma once
#include "gmp.h"


void TrialDivision(mpz_t& output, mpz_t& to_factor, int thread_count);