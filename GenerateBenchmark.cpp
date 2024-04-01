#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <cmath>
#include "gmp.h"
#include "gmpxx.h"

#include "GeneratePrimes.h"


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