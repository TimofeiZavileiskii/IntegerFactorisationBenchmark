#pragma once

#include "gmp.h"
#include "gmpxx.h"

void generate_prime(mpz_t prime, mpz_t lower_bound, mpz_t upper_bound, gmp_randstate_t state);

mpz_class generate_rsa(int bit_size, gmp_randstate_t state);