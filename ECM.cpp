#include "ECM.h"
#include "gmp.h"
#include "gmpxx.h"
#include <vector>
#include <iostream>
#include "Utils.h"
#include "EllipticCurves.h"
#include <cmath>
#include <random>
#include <thread>


void EcmJobWeisMonostage(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound, int& curves_tried){
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, starting_seed);
    
    WeistrassCurve curve;
    mpz_set(curve.n, to_factor);
    mpz_t f;
    mpz_inits(f, NULL);

    WeistrassPoint point;
    
    while(!factored)
    {
        curves_tried++;
        int outcome = curve.Initialise(random_state, point, to_factor, output);
        if(outcome == 0)
            continue;
        if(outcome == 1){
            factored = true;
            continue;
        }
        
        for(mpz_class& prime : primes){
            if(factored)
                break;

            long exponent = log_base(bound, prime.get_d());
            mpz_pow_ui(f, prime.get_mpz_t(), exponent);
            curve.MultPoints(point, f, output, factored);
        }
    }

    mpz_clears(f, NULL);
}

void EcmJobMontMonostage(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound, int& curves_tried)
    {
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, starting_seed);
    
    MontgomeryCurve curve(to_factor);
    mpz_t theta, random_bound, f, gcd;
    mpz_inits(theta, random_bound, f, gcd, NULL);

    MontgomeryPoint point;
    mpz_set(random_bound, to_factor);
    mpz_sub_ui(random_bound, random_bound, 7);
    while(!factored)
    {   
        curves_tried += 1;
        //Randomly select theta
        mpz_urandomm(theta, random_state, random_bound);
        mpz_add_ui(theta, theta, 6);

        //Generate x, z, C from theta
        bool inverse_not_exist = curve.ComputeC_Q(theta, point, output);
        if(inverse_not_exist){
            if(mpz_cmp(output, to_factor) == 0){
                continue;
            }
            else{
                return;
            }
        }

        for(mpz_class& prime : primes){
            if(factored)
                break;

            int exponent = log_base(bound, prime.get_d());
            for(int i = 0; i < exponent; i++){
                curve.MultPoints(point, prime.get_mpz_t());
            }
            //mpz_pow_ui(f, prime.get_mpz_t(), exponent);
        }
        mpz_gcd(gcd, point.z, to_factor);
        if(mpz_cmp_ui(point.z, 0) == 0){
            std::cout << "Z IS 0" << std::endl;
            continue;
        }
        if(mpz_cmp_ui(gcd, 1) > 0){
            mpz_set(output, gcd);
            factored = true;
        }
    }
    mpz_clears(theta, random_bound, f, gcd, NULL);
}


void EcmJobWeisDistage(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound, int& curves_tried){

}


void EcmJobMonstDistage(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound, int& curves_tried){

}


mpz_class ChooseBound(mpz_t to_factor){
    double to_factor_d = mpz_get_d(to_factor);
    double smallest_factor = sqrt(to_factor_d);
    double log_smallest_factor = log(smallest_factor);
    double bound = std::exp((sqrt(2.0))/2.0 * sqrt(log(log_smallest_factor)*log_smallest_factor));
    mpz_class B = bound;
    B -= B % 2;
    return B;
}


void Ecm(mpz_t output, mpz_t to_factor, int thread_count, EcmAlgorithm algorithm){
    if(thread_count == 0){
        thread_count = 8;
    }

    mpz_class B = ChooseBound(to_factor);

    std::vector<mpz_class> primes;
    SieveOfEratosthenes(B, primes);

    volatile bool factored = false;
    int sum = 0;
    if(thread_count == 1){
        long seed = rand();
        EcmJobMontMonostage(output, to_factor, factored, primes, seed, B.get_d(), sum);
    }
    else{
        std::vector<std::thread> workers = std::vector<std::thread>(); 
        int* curve_nums = new int[thread_count];
        for(int i = 0; i < thread_count; i++){
            long random_seed = rand();
            curve_nums[i] = 0;
            if(algorithm == Weierstrass1){
                workers.emplace_back(std::thread(EcmJobWeisMonostage, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d(), std::ref(curve_nums[i])));
            }
            else if(algorithm == Montgomery1){
                workers.emplace_back(std::thread(EcmJobMontMonostage, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d(), std::ref(curve_nums[i])));
            }
        }
        for(std::thread &thread : workers){
            thread.join();
        }

        for(int i = 0; i < thread_count; i++){
            sum += curve_nums[i];
        }
        delete[] curve_nums;
    }
    std::cout << "Total curve number: " << sum << "\n";
}

