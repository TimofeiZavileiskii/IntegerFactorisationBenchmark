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
    mpz_t base, bound, exponent, exponent2, base_minus_1, try_divisor;
    mpz_inits(base, bound, exponent, exponent2, base_minus_1, try_divisor, NULL);
    double bound_d = mpz_get_d(to_factor);
    double bound_d_log = std::log2(bound_d);
    
    bound_d = std::pow(bound_d, bound_power);
    mpz_set_ui(bound, bound_d);
    mpz_set_ui(base, 2);
    
    std::vector<mpz_class> primes;
    mpz_class bound_c = mpz_class(bound);
    SieveOfEratosthenes(bound_c, primes);
    unsigned int i = 0;
    while(i < primes.size()){
        mpz_class& prime = primes[i];

        double prime_d = prime.get_d();
        double exponent_d =  std::floor(bound_d_log/std::log2(prime_d));
        long exponent_l = (long)exponent_d;
        mpz_pow_ui(exponent2, prime.get_mpz_t(), exponent_l);
        mpz_powm(base, base, exponent2, to_factor);

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