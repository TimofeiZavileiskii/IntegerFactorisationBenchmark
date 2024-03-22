#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <cmath>
#include "gmp.h"
#include "gmpxx.h"


void generate_prime(mpz_t prime, mpz_t lower_bound, mpz_t upper_bound, gmp_randstate_t state){
    mpz_t gen_bound;
    mpz_init(gen_bound);
    mpz_sub(gen_bound, upper_bound, lower_bound);
    
    bool is_prime = false;
    while(!is_prime){
        mpz_urandomm(prime, state, gen_bound);
        mpz_add(prime, prime, lower_bound);
        int res = mpz_probab_prime_p(prime, 50);
        if(res > 0){
            is_prime = true;
        }
    }

    mpz_clear(gen_bound);
}


mpz_class generate_rsa(int bit_size, gmp_randstate_t state){
    mpz_t p;
    mpz_t q;
    mpz_t rsa;
    mpz_t upper_bound;
    mpz_t lower_bound;
    mpz_inits(p, q, rsa, upper_bound, lower_bound, NULL);

    mpz_set_ui(upper_bound, 2);
    mpz_pow_ui(upper_bound, upper_bound, (bit_size+1)/2);

    generate_prime(p, lower_bound, upper_bound, state);

    double p_d = mpz_get_d(p);
    int bit_size_q = bit_size - ((int)floor(log2(p_d)));
    mpz_set_ui(upper_bound, 2);
    mpz_pow_ui(upper_bound, upper_bound, bit_size_q);

    mpz_set_ui(lower_bound, 2);
    mpz_pow_ui(lower_bound, lower_bound, bit_size_q-1);

    generate_prime(q, lower_bound, upper_bound, state);
    mpz_mul(rsa, p, q);
    
    mpz_class output = mpz_class(rsa);
    mpz_clears(p, q, rsa, upper_bound, lower_bound, NULL);
    return output;
}


void generate_benchmark(int bit_bound, int seed){
    std::string filename = "rsa_numbers.csv"; 

    std::fstream benchmark;
    benchmark.open(filename, std::fstream::in | std::fstream::out | std::fstream::trunc);
    
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, seed);

    int numbers_per_bit = 10;
    for(int curr_size = 30; curr_size <= bit_bound; curr_size += 20){
        benchmark << curr_size << ",";

        for(int i = 0; i < numbers_per_bit; i++){
            benchmark << generate_rsa(curr_size, random_state).get_str();
            if(i != numbers_per_bit - 1){
                benchmark << ",";
            }
        }

        benchmark << "\n";
    }
    benchmark.close();
}


int main(int argc, char *argv[]){
    std::cout << "Program Start!\n";
    std::cout << "Enter Benchmark upper limit:\n";
    int bit_limit;
    std::cin >> bit_limit;
    std::cout << "Enter Seed:\n";
    int seed;
    std::cin >> seed;

    generate_benchmark(bit_limit, seed);
    
    std::cout << "End of program!" << std::endl;
    return 0;
}