#include "gmp.h"
#include "TrialDivision.h"
#include <thread>
#include <vector>
#include <iostream>
#include <cmath>
#include <unistd.h>     

void ThreadJobTrial(mpz_t output, mpz_t to_factor, volatile bool& factored, volatile bool& can_continue, mpz_t temp_start, mpz_t temp_end){
    mpz_t try_mod, trial_divisor, end;

    mpz_inits(try_mod, trial_divisor, end, NULL);
    mpz_set(trial_divisor, temp_start);
    mpz_set(end, temp_end);
    can_continue = true;

    while(!factored && mpz_cmp(trial_divisor, end) < 0){
 
        mpz_mod(try_mod, to_factor, trial_divisor);
        if(mpz_cmp_ui(try_mod, 0) == 0){
            factored = true;
            mpz_set(output, trial_divisor);
        }
        mpz_add_ui(trial_divisor, trial_divisor, 1);
    }
    mpz_clears(try_mod, trial_divisor, end, NULL);
}

void TrialDivision(mpz_t output, mpz_t to_factor, int thread_count){
    if(thread_count == 0){
        thread_count = 8;
    }

    volatile bool factored = false;
    mpz_t lower_bound;
    mpz_init(lower_bound);
    double lower_bound_d = mpz_get_d(to_factor);
    lower_bound_d = std::sqrt(lower_bound_d) + 1;
    mpz_set_d(lower_bound, lower_bound_d);
    if(thread_count > 1){
        mpz_t range_size, range_start, range_end;
        mpz_inits(range_size, range_start, range_end, NULL);
        mpz_div_ui(range_size, lower_bound, thread_count);
        mpz_set_ui(range_start, 2);
        mpz_add(range_end, range_start, range_size);

        std::vector<std::thread> workers;
        for(int i = 0; i < thread_count; i++){
            volatile bool can_continue = false;
            workers.emplace_back(std::thread(ThreadJobTrial, std::ref(output), std::ref(to_factor), std::ref(factored), std::ref(can_continue), range_start, range_end));
            while(!can_continue){
            }

            mpz_set(range_start, range_end);
            if(i == thread_count-2){
                mpz_set(range_end, lower_bound);}
            else{
                mpz_add(range_end, range_start, range_size);}
        }
        for(int i = 0; i < thread_count; i++){
            workers[i].join();
        }
        mpz_clears(range_size, range_start, range_end, NULL);
    }
    else{
        mpz_t start;
        mpz_init(start);
        mpz_set_ui(start, 2);
        volatile bool can_continue = true;
        ThreadJobTrial(output, to_factor, factored, can_continue, start, lower_bound);
        mpz_clear(start);
    }
    mpz_clear(lower_bound);
}