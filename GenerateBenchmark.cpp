#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <cmath>
#include "gmp.h"
#include "gmpxx.h"
#include <cxxopts.hpp>

#include "GeneratePrimes.h"


void generate_benchmark(int bit_bound, int numbers_per_size, int initial_size, int step, int seed){
    std::string filename = "rsa_numbers.csv"; 

    std::fstream benchmark;
    benchmark.open(filename, std::fstream::in | std::fstream::out | std::fstream::trunc);
    
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, seed);

    for(int curr_size = initial_size; curr_size <= bit_bound; curr_size += step){
        benchmark << curr_size << ",";

        for(int i = 0; i < numbers_per_size; i++){
            benchmark << generate_rsa(curr_size, random_state).get_str();
            if(i != numbers_per_size - 1){
                benchmark << ",";
            }
        }

        benchmark << "\n";
    }
    benchmark.close();
}

void record_benchmark_seed(unsigned int seed){
    std::ofstream seed_file;
    std::string filename = "seed.txt";
    seed_file.open(filename);
    seed_file << "Seed:" << seed << "\n";
    seed_file.close();
}

int main(int argc, char *argv[]){
    std::cout << "Program Start!\n";
    cxxopts::Options options("Factorisation Benchmark", "The framework for timing integer factorisation algorithms");
    options.add_options()
        ("m,max_size", "Up to which bit size generate benchmark", cxxopts::value<int>()->default_value("310"))
        ("n,numbers_per_size", "How many numbers per bit size generate", cxxopts::value<int>()->default_value("10"))
        ("s,seed", "Which random seed to use, uses system time if left empty", cxxopts::value<int>()->default_value("-1"))
        ("i,initial_size", "Integer size in bits to start from", cxxopts::value<int>()->default_value("30"))
        ("t,step", "What step size to make in the integer bit sizes", cxxopts::value<int>()->default_value("20"));

    auto parsed_arguments = options.parse(argc, argv);

    int bit_limit = parsed_arguments["max_size"].as<int>();
    int numbers_per_size = parsed_arguments["numbers_per_size"].as<int>();
    int seed = parsed_arguments["seed"].as<int>();
    int initial_size = parsed_arguments["initial_size"].as<int>();
    int step = parsed_arguments["step"].as<int>();

    long c_seed = time(NULL);
    srand(c_seed);

    if(seed == -1){
        seed = rand();
    }
    record_benchmark_seed(seed);

    generate_benchmark(bit_limit, numbers_per_size, initial_size, step, seed);
    
    std::cout << "End of program!" << std::endl;
    return 0;
}