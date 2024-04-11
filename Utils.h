#pragma once
#include <iostream>
#include <vector>
#include "gmp.h"
#include "gmpxx.h"
#include <cmath>

# define m_arith 0 //montgomery arithemtic off = 0 on = 1

bool inline IndexBitArray(long index, unsigned char* array);

void inline SetBitArray(long index, unsigned char* array);

void inline ZeroBitArray(long index, unsigned char* array);

void SieveOfEratosthenes(long upperbound, std::vector<long>& primes);

inline void print_mpz(std::string str, mpz_t n){
    std::cout << str;
    mpz_out_str(NULL, 10, n);
    std::cout << "\n";
}

inline void initialise_rstate(gmp_randstate_t& rstate, long rseed){
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, rseed);
}

struct MontgomeryParams{
    mpz_t mod, n_bar, r_mask, temp;
    uint bit_num;

    inline MontgomeryParams(){
        mpz_inits(mod, n_bar, r_mask, temp, NULL);
    }

    inline ~MontgomeryParams(){
        mpz_clears(mod, n_bar, r_mask, temp, NULL);
    }
};

inline void MontgomerySetup(mpz_t mod, MontgomeryParams& params){
    mpz_set(params.mod, mod);
    params.bit_num = (uint)log2(mpz_get_d(params.mod))+2;
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

inline void ModularAdd(mpz_t c, mpz_t a, mpz_t b, mpz_t mod){
    mpz_add(c, a, b);
    if(mpz_cmp(c, mod) > -1){
        mpz_sub(c, c, mod);
    }
}

inline void TransformIn(mpz_t a, MontgomeryParams& params){
    #if m_arith==1
        mpz_mul_2exp(a, a, params.bit_num);
        mpz_mod(a, a, params.mod);
    #endif
}

inline void TransformOut(mpz_t a, MontgomeryParams& params){
    #if m_arith==1
        MontgomeryReduction(a, params);
    #endif
}

inline void MontgomeryMul(mpz_t c, mpz_t a, mpz_t b, MontgomeryParams& params){
    mpz_mul(c, a, b);
    #if m_arith==1
        MontgomeryReduction(c, params);
    #endif
    #if m_arith==0
        mpz_mod(c, c, params.mod);
    #endif
}

inline double log_base(double x, double base) {
    return log(x) / log(base);
}