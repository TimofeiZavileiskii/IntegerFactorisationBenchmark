#pragma once
#include "gmp.h"
#include "Utils.h"


bool PollardsP1Job(mpz_t& output, mpz_t& to_factor, double bound_power);


void PollardsP1(mpz_t& output, mpz_t& to_factor, double bound_exp);