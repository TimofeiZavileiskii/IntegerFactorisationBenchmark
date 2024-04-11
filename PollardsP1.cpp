#include "gmp.h"
#include "gmpxx.h"
#include "Utils.h"
#include "PollardsP1.h"
#include <cmath>
#include <iostream>
/*
    a = 2
    B = int(math.pow(n, 1/3))
    primes = sieve_of_eratosthenes(B)

    for i, p in enumerate(primes):
        e = int(math.floor(math.log2(B)/math.log2(p)))
        f = pow(p, e)
        a = pow(a, f, n)
        if i % 10 == 0:
            q = math.gcd(a-1, n)
            if q > 1:
                return q
    return 1
*/


bool PollardsP1Job(mpz_t output, mpz_t to_factor, double bound_power){
    mpz_t base, base_minus_1, try_divisor;
    long bound, exponent;
    mpz_inits(base, exponent, base_minus_1, try_divisor, NULL);
    double bound_d = mpz_get_d(to_factor);
    double bound_d_log = std::log2(bound_d);
    
    bound_d = std::pow(bound_d, bound_power);
    bound = bound_d;
    mpz_set_ui(base, 2);
    
    std::vector<long> primes;
    SieveOfEratosthenes(bound, primes);
    unsigned int i = 0;
    while(i < primes.size()){
        long& prime = primes[i];

        double prime_d = prime;
        double exponent_d =  std::floor(bound_d_log/std::log2(prime_d));
        long exponent_l = (long)exponent_d;
        exponent = (pow(prime, exponent_l) + 0.5);
        mpz_powm_ui(base, base, exponent, to_factor);

        if(i % 10 == 0){
            mpz_sub_ui(base_minus_1, base, 1);
            mpz_gcd(try_divisor, base_minus_1, to_factor);
            if(mpz_cmp_ui(try_divisor, 1) > 0){
                mpz_set(output, try_divisor);
                return true;
            }
        }
        i++;
    }

    return false;
}


void PollardsP1(mpz_t output, mpz_t to_factor, int thread_count){
    bool factored = false;
    double bound = 100000;
    //Change the algorithm for the thread count
    factored = PollardsP1Job(output, to_factor, bound);
    
    if(!factored){
        mpz_set_ui(output, 1);
        std::cout << "FAILURE TO FACTORISE\n";
    }
}