#pragma once
#include <iostream>
#include <vector>
#include "gmp.h"
#include "gmpxx.h"
#include <cmath>

bool inline IndexBitArray(long index, unsigned char* array);

void inline SetBitArray(long index, unsigned char* array);

void inline ZeroBitArray(long index, unsigned char* array);

void SieveOfEratosthenes(mpz_class upperbound, std::vector<mpz_class>& primes);

inline void print_mpz(std::string str, mpz_t n){
    std::cout << str;
    mpz_out_str(NULL, 10, n);
    std::cout << "\n";
}

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

inline void ModularSub(mpz_t c, mpz_t a, mpz_t b, mpz_t mod){
    mpz_sub(c, a, b);
    if(mpz_cmp_ui(c, 0) < 0){
        mpz_add(c, c, mod);
    }
}

inline void TransformIn(mpz_t a, MontgomeryParams& params){
    mpz_mul_2exp(a, a, params.bit_num);
    mpz_mod(a, a, params.mod);
}

inline void TransformOut(mpz_t a, MontgomeryParams& params){
    MontgomeryReduction(a, params);
}

inline void MontgomeryMul(mpz_t c, mpz_t a, mpz_t b, MontgomeryParams& params){
    mpz_mul(c, a, b);
    mpz_mod(c, c, params.mod);
    //MontgomeryReduction(c, params);
}