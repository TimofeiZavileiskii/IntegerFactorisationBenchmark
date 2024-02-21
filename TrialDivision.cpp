#include "gmp.h"
#include "TrialDivision.h"
#include <thread>
#include <vector>
#include <iostream>
#include <cmath>

void ThreadJobTrial(mpz_t& output, mpz_t& to_factor, volatile bool& factored, mpz_t start, mpz_t end){
    mpz_t try_mod, trial_divisor;
    mpz_inits(try_mod, trial_divisor, NULL);
    mpz_set(trial_divisor, start);

    while(!factored && mpz_cmp(trial_divisor, end) < 0){
        mpz_mod(try_mod, to_factor, trial_divisor);
        if(mpz_cmp_ui(try_mod, 0) == 0){
            factored = true;
            mpz_set(output, trial_divisor);
        }
        mpz_add_ui(trial_divisor, trial_divisor, 1);
    }
}

void TrialDivision(mpz_t& output, mpz_t& to_factor, int thread_count){
    volatile bool factored = false;
    mpz_t lower_bound;
    mpz_init(lower_bound);

    double lower_bound_d = mpz_get_d(to_factor);
    lower_bound_d = std::sqrt(lower_bound_d);
    mpz_set_d(lower_bound, lower_bound_d);

    if(thread_count > 1){
        mpz_t range_size, range_start, range_end;
        mpz_inits(range_size, range_start, range_end, NULL);
        mpz_div_ui(range_size, lower_bound, thread_count);
        mpz_out_str(NULL, 10, range_size);
        mpz_set_ui(range_start, 2);
        mpz_add(range_end, range_start, range_size);

        std::vector<std::thread> workers;
        for(int i = 0; i < thread_count; i++){
            workers.emplace_back(std::thread(ThreadJobTrial, std::ref(output), std::ref(to_factor), std::ref(factored), range_start, range_end));
            mpz_set(range_start, range_end);
            if(i != thread_count-1){
                mpz_add(range_end, range_start, range_size);}
            else{
                mpz_set(range_end, lower_bound);}
        }

        for(int i = 0; i < thread_count; i++){
            workers[i].join();
        }
    }
    else{
        mpz_t start, end;
        mpz_inits(start, end, NULL);
        mpz_set_ui(start, 2);
        ThreadJobTrial(output, to_factor, factored, start, lower_bound);
    }
}