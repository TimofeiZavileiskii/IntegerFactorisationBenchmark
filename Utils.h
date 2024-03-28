#pragma once
#include <vector>
#include "gmp.h"
#include "gmpxx.h"
#include <cmath>

bool inline IndexBitArray(long index, unsigned char* array);

void inline SetBitArray(long index, unsigned char* array);

void inline ZeroBitArray(long index, unsigned char* array);

void SieveOfEratosthenes(mpz_class upperbound, std::vector<mpz_class>& primes);

struct MontgomeryParams{
    mpz_t mod, n_bar, r_mask, temp;
    uint bit_num;

    MontgomeryParams(){
        mpz_inits(mod, n_bar, r_mask, temp, NULL);
    }

    ~MontgomeryParams(){
        mpz_clears(mod, n_bar, r_mask, temp, NULL);
    }
};

inline void MontgomerySetup(mpz_t mod, MontgomeryParams& params){
    mpz_set(params.mod, mod);
    params.bit_num = (uint)log2(mpz_get_d(params.mod))+1;
    mpz_set_ui(params.r_mask, 1);
    mpz_mul_2exp(params.r_mask, params.r_mask,params. bit_num);
    mpz_invert(params.n_bar, params.r_mask, params.mod);
    mpz_mul(params.n_bar, params.n_bar, params.r_mask);
    mpz_sub_ui(params.n_bar, params.n_bar, 1);
    mpz_div(params.n_bar, params.n_bar, params.mod);
    mpz_sub_ui(params.r_mask, params.r_mask, 1);
}

inline void MontgomeryReduction(mpz_t a, MontgomeryParams& params){
    mpz_and(params.temp, a, params.r_mask);
    mpz_mul(params.temp, params.temp, params.n_bar);
    mpz_and(params.temp, params.temp, params.r_mask);
    mpz_mul(params.temp, params.temp, params.mod);
    mpz_add(params.temp, a, params.temp);
    mpz_div_2exp(a, params.temp, params.bit_num);
    if(mpz_cmp(a, params.mod) != -1){
        mpz_sub(a, a, params.mod);
    }
}