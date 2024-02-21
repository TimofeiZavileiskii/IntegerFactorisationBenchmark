#pragma once
#include "gmp.h"


void ThreadJobTrial(mpz_t& output, mpz_t& to_factor, volatile bool& factored, mpz_t start, mpz_t end);


void TrialDivision(mpz_t& output, mpz_t& to_factor, int thread_count);