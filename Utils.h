#pragma once
#include <vector>
#include "gmp.h"
#include "gmpxx.h"


bool inline IndexBitArray(long index, unsigned char* array);

void inline SetBitArray(long index, unsigned char* array);

void inline ZeroBitArray(long index, unsigned char* array);

void SieveOfEratosthenes(mpz_class upperbound, std::vector<mpz_class>& primes);