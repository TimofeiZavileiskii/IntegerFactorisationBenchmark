#include "gmp.h"
#include "gmpxx.h"
#include <random>

#include <thread>
#include <iostream>
#include <stdlib.h>
#include "PollardsRho.h"
#include "Utils.h"

inline void hash_function(mpz_t x, mpz_t one, MontgomeryParams& params){
    MontgomeryMul(x, x, x, params);
    mpz_add_ui(x, x, 1);
}

void ThreadJobRho(mpz_t output, mpz_t to_factor, volatile bool& factored, int seed){
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, seed);

    //Initialise seeds
    mpz_t u, v, diff, factor, one, product, random_bound;
    MontgomeryParams params;
    //Initialise montgomery multiplication
    MontgomerySetup(to_factor, params);
    mpz_inits(u, v, diff, factor, one, product, random_bound, NULL);
    mpz_set_ui(one, 1);
    TransformIn(one, params);
    mpz_set(random_bound, to_factor);
    mpz_sub_ui(random_bound, random_bound, 1);
    mpz_urandomm(u, random_state, random_bound);
    mpz_set(v, u);

    int count = 0;
    int max_count = (int)log_base(mpz_get_d(to_factor), 2.0) * 5;

    while(!factored){
        //u = f(u)  v = f(f(v))
        hash_function(u, one, params);
        hash_function(v, one, params);
        hash_function(v, one, params);

        ModularSub(diff, u, v, params.mod);

        mpz_mul(product, product, diff);
        mpz_mod(product, product, to_factor);
        count++;
        //Check if factor is found
        if(mpz_cmp_ui(diff, 0) == 0){
            mpz_urandomm(u, random_state, random_bound);
            mpz_set(v, u);
        }

        if(count >= max_count){
            count = 0;

            mpz_gcd(factor, product, to_factor);
            if(mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, to_factor) < 0){
                mpz_set(output, factor);
                factored = true;
            }
            mpz_set_ui(product, 1);
        }
    }
    mpz_clears(u, v, diff, factor, one, product, random_bound, NULL);
}


void PollardsRho(mpz_t output, mpz_t to_factor, int thread_count){
    volatile bool factored = false;

    if(thread_count > 1){
        std::vector<std::thread> workers = std::vector<std::thread>(); 

        for(int i = 0; i < thread_count; i++){
            int random_seed = rand();
            workers.emplace_back(std::thread(ThreadJobRho, std::ref(output), std::ref(to_factor), std::ref(factored), random_seed));
        }

        for(std::thread &thread : workers){
            thread.join();
        }
    }
    else{
        int random_seed = rand();
        ThreadJobRho(output, to_factor, factored, random_seed);
    }
}