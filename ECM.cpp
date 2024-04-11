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
    initialise_rstate(random_state, starting_seed);
    
    WeistrassCurve curve;
    WeistrassPoint point;
    mpz_set(curve.n, to_factor);
    mpz_t f;
    mpz_init(f);
    
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

    mpz_clear(f);
}


void EcmJobWeisDistage(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound1, double bound2, int& curves_tried){
    const int DIFF_COUNT = 80;
    
    gmp_randstate_t random_state;
    initialise_rstate(random_state, starting_seed);
    
    WeistrassCurve curve;
    WeistrassPoint original_point, point;
    mpz_set(curve.n, to_factor);
    mpz_t f;
    WeistrassPoint differences[DIFF_COUNT];
    mpz_init(f);
    
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
        
        int index;
        for(index = 0; primes[index].get_d() < bound1; index++){
            if(factored)
                break;

            long exponent = log_base(bound2, primes[index].get_d());
            mpz_pow_ui(f, primes[index].get_mpz_t(), exponent);
            curve.MultPoints(point, f, output, factored);
        }

        //Precompute the prime differences
        point.Copy(original_point);
        point.Copy(differences[0]);
        curve.AddPoints(differences[0], differences[0], output, factored);
        differences[0].Copy(point);
        for(int i = 1; i < DIFF_COUNT && !factored; i++){
            curve.AddPoints(point, differences[0], output, factored);
            point.Copy(differences[i]);
        }

        for(; index < primes.size() && !factored; index++){
            long diff = (primes[index].get_ui() - primes[index-1].get_ui()) >> 1; //Set prime to difference in place and get long and div by 2      
            curve.AddPoints(original_point, differences[diff], output, factored);
        }
    }

    mpz_clear(f);
}


void EcmJobMontMonostage(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound, int& curves_tried){
    gmp_randstate_t random_state;
    initialise_rstate(random_state, starting_seed);
    
    MontgomeryCurve curve(to_factor);
    mpz_t theta, random_bound, gcd;
    mpz_inits(theta, random_bound, gcd, NULL);

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
    mpz_clears(theta, random_bound, gcd, NULL);
}


void EcmJobMontDistage(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound1, double bound2, int& curves_tried){
    gmp_randstate_t random_state;
    initialise_rstate(random_state, starting_seed);

    
    MontgomeryCurve curve(to_factor);
    mpz_t theta, random_bound, f, gcd, product, multiple, point_r_xz, diff_x, sum_z;
    mpz_inits(theta, random_bound, f, gcd, multiple, point_r_xz, diff_x, sum_z, product, NULL);

    const int DIFF_SIZE = 5;

    MontgomeryPoint point, point_t, point_r, temp_p;
    MontgomeryPoint differences[DIFF_SIZE];
    mpz_t x_z_products[DIFF_SIZE];
    for(int i = 0; i < DIFF_SIZE; i++){
        mpz_init(x_z_products[i]);
    }

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
            if(mpz_cmp(output, to_factor) == 0){ continue; }
            else{ return; }
        }

        //stage 1
        int index;
        for(index = 0; primes[index].get_d() < bound1 && !factored; index++){
            int exponent = log_base(bound1, primes[index].get_d());
            mpz_pow_ui(f, primes[index].get_mpz_t(), exponent);
            curve.MultPoints(point, f);
        }

        mpz_gcd(gcd, point.z, to_factor);
        if(mpz_cmp_ui(point.z, 0) == 0){
          //  std::cout << "Z IS 0\n";
            continue;
        }
        if(mpz_cmp_ui(gcd, 1) > 0 && mpz_cmp(gcd, to_factor) < 0){
            mpz_set(output, gcd);
           // std::cout << "Factored in stage 1\n";
            factored = true;
            continue;
        }

        //initialise array of differences
        point.Copy(differences[0]);
        curve.DoublePoint(differences[0]);
        differences[0].Copy(differences[1]);
        curve.DoublePoint(differences[1]);
        for(int i = 2; i < DIFF_SIZE; i++){
            differences[i-1].Copy(differences[i]);
            curve.AddPoints(differences[0], differences[i], differences[i-2]);
        }
        for(int i = 0; i < DIFF_SIZE; i++){
            differences[i].ComputeXZ(x_z_products[i], to_factor);
        }

        //Investigate performance issues with stage 2 - rn it seems to work only by accident
        mpz_set_ui(product, 1);
        point.Copy(point_r);
        point.Copy(point_t);
        mpz_set_d(multiple, bound1-1);
        curve.MultPoints(point_r, multiple);
        mpz_sub_ui(multiple, multiple, 2*DIFF_SIZE);
        curve.MultPoints(point_t, multiple);

        //Perform stage 2
        for(int i = bound1-1; index < primes.size() && !factored; i += 2*DIFF_SIZE){
            point_r.ComputeXZ(point_r_xz, to_factor);
            for(; primes[index] <= (i + 2*DIFF_SIZE) && index < primes.size(); index++){ //Check that iteration is correct
                int dist = ((primes[index].get_ui() - i) >> 1) - 1;
                mpz_sub(diff_x, point_r.x, differences[dist].x);
                mpz_add(sum_z, point_r.z, differences[dist].z);
                mpz_mul(diff_x, diff_x, sum_z);
                mpz_mod(diff_x, diff_x, to_factor);
                mpz_sub(diff_x, diff_x, point_r_xz);
                mpz_add(diff_x, diff_x, x_z_products[dist]);
                mpz_mul(product, product, diff_x);
                mpz_mod(product, product, to_factor);
            }
            point_r.Copy(temp_p);
            curve.AddPoints(differences[DIFF_SIZE-1], point_r, point_t);
            temp_p.Copy(point_t);
        }

        mpz_gcd(gcd, product, to_factor);
        if(mpz_cmp_ui(point.z, 0) == 0){
           // std::cout << "Z IS 0\n";
            continue;
        }
        if(mpz_cmp_ui(gcd, 1) > 0 && mpz_cmp(gcd, to_factor) < 0){
            mpz_set(output, gcd);
            //std::cout << "Factored in stage 2\n";
            factored = true;
        }
    }

    for(int i = 0; i < DIFF_SIZE; i++){ 
        mpz_clear(x_z_products[i]);
    }
    mpz_clears(theta, random_bound, f, gcd, product, multiple, point_r_xz, diff_x, sum_z, NULL);
}


mpz_class ChooseBound(mpz_t to_factor){
    double to_factor_d = mpz_get_d(to_factor);
    double smallest_factor = sqrt(to_factor_d);
    double log_smallest_factor = log(smallest_factor);
    double bound = std::exp((sqrt(2.0))/2.0 * sqrt(log(log_smallest_factor)*log_smallest_factor));
    std::cout << "Bound is " << bound << "\n";
    mpz_class B = bound;
    B -= B % 2;
    return B;
}


void Ecm(mpz_t output, mpz_t to_factor, int thread_count, EcmAlgorithm algorithm){
    if(thread_count == 0){
        thread_count = 8;
    }

    mpz_class B;
    mpz_class B2;
    
    std::vector<mpz_class> primes;
    if(algorithm == Weierstrass1 || algorithm == Montgomery1){
        B = ChooseBound(to_factor);
        SieveOfEratosthenes(B, primes);
    }
    else{
        B = ChooseBound(to_factor);
        B2 = B*100;
        SieveOfEratosthenes(B2, primes);
        if(B < 260){
            B = 260;
            B2 = B * 100;
        }
    }
    

    volatile bool factored = false;
    int sum = 0;
    if(thread_count == 1){
        long seed = rand();
        switch (algorithm)
            {
            case Weierstrass1:
                EcmJobWeisMonostage(output, to_factor, factored, primes, seed, B.get_d(), sum);
                break;
            case Montgomery1:
                EcmJobMontMonostage(output, to_factor, factored, primes, seed, B.get_d(), sum);
                break;
            case Weierstrass2:
                EcmJobWeisDistage(output, to_factor, factored, primes, seed, B.get_d(), B2.get_d(), sum);
                break;
            case Montgomery2:
                EcmJobMontDistage(output, to_factor, factored, primes, seed, B.get_d(), B2.get_d(), sum);
                break;
            }
    }
    else{
        std::vector<std::thread> workers = std::vector<std::thread>(); 
        int* curve_nums = new int[thread_count];
        for(int i = 0; i < thread_count; i++){
            long random_seed = rand();
            curve_nums[i] = 0;
            switch (algorithm)
            {
            case Weierstrass1:
                workers.emplace_back(std::thread(EcmJobWeisMonostage, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d(), std::ref(curve_nums[i])));
                break;
            case Montgomery1:
                workers.emplace_back(std::thread(EcmJobMontMonostage, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d(), std::ref(curve_nums[i])));
                break;
            case Weierstrass2:
                workers.emplace_back(std::thread(EcmJobWeisDistage, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d(), B2.get_d(), std::ref(curve_nums[i])));
                break;
            case Montgomery2:
                workers.emplace_back(std::thread(EcmJobMontDistage, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d(), B2.get_d(), std::ref(curve_nums[i])));
                break;
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
