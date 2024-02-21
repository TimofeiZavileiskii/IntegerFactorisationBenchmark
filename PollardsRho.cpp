#include "gmp.h"
#include <iostream>
#include <stdlib.h>
#include <random>
#include <thread>
#include "PollardsRho.h"


void ThreadJobRho(mpz_t& output, mpz_t& to_factor, volatile bool& factored, int seed){
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, seed);

    mpz_t u, v, diff, factor, random_bound;
    mpz_inits(u, v, diff, factor, random_bound, NULL);
    mpz_set(random_bound, to_factor);
    mpz_sub_ui(random_bound, random_bound, 1);

    mpz_urandomm(u, random_state, random_bound);
    mpz_set(u, v);

    while(!factored){    
        mpz_mul(u, u, u);
        mpz_add_ui(u, u, 1);
        mpz_mod(u, u, to_factor);

        mpz_mul(v, v, v);
        mpz_add_ui(v, v, 1);
        mpz_mod(v, v, to_factor);
        mpz_mul(v, v, v);
        mpz_add_ui(v, v, 1);
        mpz_mod(v, v, to_factor);
        
        mpz_sub(diff, u, v);
        mpz_mod(diff, diff, to_factor);
        if(mpz_cmp_ui(diff, 0) == 0){
            mpz_urandomm(u, random_state, random_bound);
            mpz_set(u, v);
        }
        else{
            mpz_gcd(factor, diff, to_factor);
            if(mpz_cmp_ui(factor, 1) > 0){
                mpz_set(output, factor);
                factored = true;
            }
        }
    }
}


void PollardsRho(mpz_t& output, mpz_t& to_factor, int thread_count){
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