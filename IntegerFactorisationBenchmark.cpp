#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>
#include <time.h>
#include "gmp.h"
#include "Utils.h"
#include "PollardsRhoCuda.h"
#include "PollardsRho.h"
#include "ECM.h"
#include "ECMCuda.h"
#include "PollardsP1.h"
#include "TrialDivisionCuda.h"
#include "TrialDivision.h"
#include "FactorisationStats.h"
#include <pari/pari.h>

#include <cxxopts.hpp>
#include <unistd.h>
#include <cstdlib>
#include <signal.h>
#include <set>


enum AlgorithmType{
    TRIAL_DIVISION,
    TRIAL_DIVISION_CUDA,
    POLLARDS_RHO,
    POLLARDS_RHO_CUDA,
    POLLARDS_P1,
    ECM_WEIERSTRASS,
    ECM_MONTGOMERY,
    ECM_WEIERSTRASS2,
    ECM_MONTGOMERY2,
    ECM_CUDA,
    PARI
};


AlgorithmType ParseAlgorithm(const std::string& algorithm_name){
    if(algorithm_name == "trial_division" || algorithm_name == "td"){
        return TRIAL_DIVISION;
    }
    else if(algorithm_name == "trial_division_cuda" || algorithm_name == "tdc"){
        return TRIAL_DIVISION_CUDA;
    }
    else if(algorithm_name == "pollards_rho" || algorithm_name == "pr"){
        return POLLARDS_RHO;
    }
    else if(algorithm_name == "pollards_rho_cuda" || algorithm_name == "prc"){
        return POLLARDS_RHO_CUDA;
    }
    else if(algorithm_name == "pollards_p1" || algorithm_name == "pp1"){
        return POLLARDS_P1;
    }
    else if(algorithm_name == "ecm_weierstrass" || algorithm_name == "ecmw"){
        return ECM_WEIERSTRASS;
    }
    else if(algorithm_name == "ecm_montgomery" || algorithm_name == "ecmm"){
        return ECM_MONTGOMERY;
    }
    else if(algorithm_name == "ecm_weierstrass_2" || algorithm_name == "ecmw2"){
        return ECM_WEIERSTRASS2;
    }
    else if(algorithm_name == "ecm_montgomery_2" || algorithm_name == "ecmm2"){
        return ECM_MONTGOMERY2;
    }
    else if(algorithm_name == "ecm_cuda" || algorithm_name == "ecmc"){
        return ECM_CUDA;
    }
    else if(algorithm_name == "pari"){
        return PARI;
    }
    else{
        throw std::logic_error("Algorithm with name " + algorithm_name + " not implemented\n");
    }
}


FactorisationStats* GetStatsObj(AlgorithmType type, int thread_count){
    std::string stats_filename = "AlgorithmStats.txt";

    if(type == ECM_WEIERSTRASS || type == ECM_WEIERSTRASS2 || type == ECM_MONTGOMERY || type == ECM_MONTGOMERY2){
        return new StatsEcm(stats_filename, thread_count);
    }
    else{
        return new FactorisationStats(stats_filename);
    }
}


void PariFactorise(mpz_t output, mpz_t to_factor){
    char *str_input = mpz_get_str(NULL, 10, to_factor);
    GEN pari_rsa = gp_read_str(str_input);
    GEN factors = Z_factor(pari_rsa);
    GEN factors_ints = gel(factors, 1);
    GEN factor1 = gel(factors_ints, 1);
    char *str_output = GENtostr(factor1);
    mpz_set_str(output, str_output, 10);
}


void run_algorithm(mpz_t output, mpz_t to_factor, int thread_count, AlgorithmType type, const std::vector<long>& primes, FactorisationStats* stats){
    switch(type){
        case TRIAL_DIVISION:
            TrialDivision(output, to_factor, thread_count);
            break;
        case TRIAL_DIVISION_CUDA:
            TrialDivisionCuda(output, to_factor, thread_count);
            break;
        case POLLARDS_RHO:
            PollardsRho(output, to_factor, thread_count);
            break;
        case POLLARDS_RHO_CUDA:
            PollardsRhoCuda(output, to_factor, thread_count);
            break;
        case POLLARDS_P1:
            PollardsP1(output, to_factor, thread_count);
            break;
        case ECM_WEIERSTRASS:
            Ecm(output, to_factor, thread_count, Weierstrass1, primes, (StatsEcm*)stats);
            break;
        case ECM_WEIERSTRASS2:
            Ecm(output, to_factor, thread_count, Weierstrass2, primes, (StatsEcm*)stats);
            break;
        case ECM_MONTGOMERY:
            Ecm(output, to_factor, thread_count, Montgomery1, primes, (StatsEcm*)stats);
            break;
        case ECM_MONTGOMERY2:
            Ecm(output, to_factor, thread_count, Montgomery2, primes, (StatsEcm*)stats);
            break;
        case ECM_CUDA:
            EcmCuda(output, to_factor, thread_count, primes);
            break;
        case PARI:
            PariFactorise(output, to_factor);
            break;
    }
}


void PariSetup(){
    ulong prime_limit = 1<<30;
    pari_init(5000000000, prime_limit);
}


float benchmark_number(mpz_t to_factor, int thread_count, AlgorithmType algorithm_type, const std::vector<long>& primes, FactorisationStats* stats){
    mpz_t factor, other_factor, remainder;
    mpz_inits(factor, other_factor, remainder, NULL);

    std::cout << "Number to factor: ";
    mpz_out_str(NULL, 10, to_factor);
    std::cout << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    run_algorithm(factor, to_factor, thread_count, algorithm_type, primes, stats);
    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<float> time = end_time - start_time;
    if(mpz_cmp_ui(factor, 1) == 0 || mpz_cmp_ui(factor, 0) == 0){
        std::cout << "Factorisation failed for main factor = ";
        mpz_out_str(NULL, 10, factor);
        std::cout << "\n";
        throw std::logic_error("Factorisation failed");
    }
    mpz_div(other_factor, to_factor, factor);
    mpz_mod(remainder, to_factor, factor);
    if(mpz_cmp_ui(other_factor, 1) == 0 || mpz_cmp_ui(other_factor, 0) == 0){
        std::cout << "Factorisation failed for other factor = ";
        mpz_out_str(NULL, 10, other_factor);
        std::cout << "\n";
        throw std::logic_error("Factorisation failed");
    }
    if(mpz_cmp_ui(remainder, 0) == 1){
        std::cout << "Factorisation failed for remainder = \n";
        mpz_out_str(NULL, 10, remainder);
        std::cout << "\n";
        throw std::logic_error("Factorisation failed");
    }
    std::cout << "Number factorised: ";
    mpz_out_str(NULL, 10, to_factor);
    std::cout << " = ";
    mpz_out_str(NULL, 10, factor);
    std::cout << " * ";
    mpz_out_str(NULL, 10, other_factor);
    std::cout << std::endl;

    std::cout << "Time spent: " << time.count() << " s" << std::endl;
    
    stats->OutputExtraStatsNumber();

    mpz_clears(factor, other_factor, remainder, NULL);

    return time.count();
}

void manual_input(AlgorithmType type, int thread_num, const std::vector<long>& primes, FactorisationStats* stats){
    mpz_t to_factor;
    mpz_inits(to_factor, NULL);
    bool to_continue = true;

    std::string response;

    while(to_continue){
        std::cout << "Enter number to factor:\n";
        mpz_inp_str(to_factor, NULL, 10);

        benchmark_number(to_factor, thread_num, type, primes, stats);

        std::cout << "Continue? (yes/y)?\n";
        std::cin >> response;

        to_continue = response == "yes" || response == "y";
    }
}

void write_entry(auto entry, std::string end_symbol, std::string filename){
    std::fstream benchmark_result;
    benchmark_result.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
    if (!benchmark_result ){
        benchmark_result.open(filename,  std::fstream::in | std::fstream::out | std::fstream::trunc);
    }
    benchmark_result << entry << end_symbol;
    benchmark_result.close();
}

void write_random_seed(unsigned int rseed, std::string& filename){
    std::fstream benchmark_result;
    benchmark_result.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
    if (!benchmark_result){
        benchmark_result.open(filename,  std::fstream::in | std::fstream::out | std::fstream::trunc);
    } 
    benchmark_result << "Seed," << rseed << "\n";
    benchmark_result.close();
}

void compute_stats(std::vector<float>& times, float& mean, float& variance, float& std){
    int n = times.size();
    mean = 0;
    variance = 0;
    std = 0;

    for(float time : times){
        mean += time;
    }
    mean /= (float)n;

    for(float time : times){
        variance += (time - mean) * (time - mean);
    }
    variance /= (float)(n - 1);
    std = sqrt(variance);
}


void signal_handler(int signum){
    std::cout << "\nInterrupt Caught -- adding stats to results and terminating\n";
    std::string results_path = "benchmark_results.csv";
    std::string amended_results_path = "amended_benchmark_results.csv";

    //Read entire file
    std::ifstream benchmark(results_path);
    std::vector<std::string> benchmark_lines;
    std::vector<float> times;
    if(!benchmark.is_open()) {
        std::cout << "Program terminate\nFailed to write stats" << std::endl;
        exit(signum);
    }
    
    while(!benchmark.eof()) {
        std::string buffer;
        getline(benchmark, buffer);
        benchmark_lines.push_back(buffer + '\n');
    }
    benchmark.close();

    //Add stats to the last line
    int end_index = benchmark_lines.size() - 1;
    std::string last_line = benchmark_lines[end_index];
    last_line.pop_back();
    last_line.pop_back();
    std::stringstream ss(last_line);
    std::string substr;
    bool first = true;
    while(ss.good()){
        getline(ss, substr, ',');
        if(first || substr.size() == 0){
            first = false;
            continue;
        }
        times.push_back(std::stof(substr));
    }
    float mean;
    float variance;
    float std;
    compute_stats(times, mean, variance, std);

    last_line = last_line + "," + std::to_string(mean) + "," + std::to_string(std);

    benchmark_lines[end_index] = last_line;

    //Write the end buffer
    std::ofstream amended_results(amended_results_path);
    for(std::string line: benchmark_lines) {
        amended_results << line;
    }
    amended_results.close();

    exit(signum);
}


void factorise_benchmark(AlgorithmType factorisation_algorithm, int thread_num, int to_factorise_count, int rseed, const std::vector<long>& primes, FactorisationStats* stats){
    std::vector<int> bit_sizes;
    
    std::string benchmark_name = "rsa_numbers.csv";
    std::string results_path = "benchmark_results.csv";
    std::ifstream benchmark(benchmark_name);

    write_random_seed(rseed, results_path);
    
    if(!benchmark.is_open()){
        throw std::runtime_error("Benchmark file failed to open");
    }

    std::string line;
    mpz_t rsa_num;
    mpz_init(rsa_num);
    int count = 0;
    int max_count = to_factorise_count;

    while(std::getline(benchmark, line))
    {
        std::vector<float> times_for_bitsize;

        std::stringstream ss(line);
        
        std::string substr;
        std::getline(ss, substr, ',');
        int bit_size = stoi(substr);
        bit_sizes.push_back(bit_size);

        write_entry(bit_size, ",", results_path);

        while (ss.good()) {
            std::string substr;
            std::getline(ss, substr, ',');

            mpz_set_str(rsa_num, substr.c_str(), 10);
            times_for_bitsize.push_back(benchmark_number(rsa_num, thread_num, factorisation_algorithm, primes, stats));
            count++;
            write_entry(times_for_bitsize[times_for_bitsize.size()-1], ",", results_path);
            if(count > max_count){
                break;
            }
        }

        float average, variance, std;
        compute_stats(times_for_bitsize, average, variance, std);
        std::cout << "------------------------- Bit Size: " << bit_size << " Average Time: " << average << std::endl;
        stats->OutputStatsAverage();

        write_entry(average, ",", results_path);
        write_entry(std, "\n", results_path);
        if(count > max_count){
            break;
        }
    }
    mpz_clear(rsa_num);
}

//Function to find differences between prime numbers up to a bound and print them to console
void study_prime_numbers(){
    std::vector<long> primes;
    long bound = 100000u;
    SieveOfEratosthenes(bound, primes);

    std::set<long> differences;
    
    bool first = true;
    long prev;
    long diff;
    long max_diff;
    for(long& p : primes){
        if(first){
            first = false;
            prev = p;
            continue;
        }
        diff = p - prev;
        differences.insert(diff);
        if(diff > max_diff){
            max_diff = diff;
        }
        prev = p;
    }

    std::cout << "Count is " << differences.size() << " with max diff: " << max_diff << "\n";
}

int main(int argc, char *argv[]){
    std::cout << "Program Start!\n";

    signal(SIGINT, signal_handler);

    cxxopts::Options options("Factorisation Benchmark", "The framework for timing integer factorisation algorithms");
    options.add_options()
        ("a,algorithm", "Factorisation algorithm to use", cxxopts::value<std::string>())
        ("m,manual", "To run the program in the manual input mode", cxxopts::value<bool>()->default_value("false"))
        ("c,count", "Cout of integers to factorise before closing the program", cxxopts::value<int>()->default_value("1000"))
        ("s,seed", "Random seed, chosen to be system time if left empty", cxxopts::value<int>()->default_value("-1"))
        ("t,thread_number", "Number of threads to be used for the algorithm", cxxopts::value<int>()->default_value("0"))
        ("p,prime_bound", "Up to which bound prime table should be generated", cxxopts::value<long>()->default_value("50000000")); //Set to 0, signifying default number of threads, determined individually for each algorithm 

    auto parsed_arguments = options.parse(argc, argv);

    std::string factorisation_algorithm = parsed_arguments["algorithm"].as<std::string>();
    bool manual_mode = parsed_arguments["manual"].as<bool>();
    int count = parsed_arguments["count"].as<int>();
    int thread_number = parsed_arguments["thread_number"].as<int>();
    int seed = parsed_arguments["seed"].as<int>();
    long prime_bound = parsed_arguments["prime_bound"].as<long>();

    if(seed == -1){
        seed = time(NULL);
    }
    srand(seed);

    AlgorithmType algorithm_type = ParseAlgorithm(factorisation_algorithm);
    
    if(algorithm_type == PARI){
        PariSetup();
    }

    std::vector<long> prime_table;

    std::cout << "Start computing the prime number table\n";
    auto start_time = std::chrono::high_resolution_clock::now();
    SieveOfEratosthenes(prime_bound, prime_table);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> time = end_time - start_time;
    std::cout << "Finished constructing the prime number table\nIt contains " << prime_table.size() << " primes, time taken: " << time.count() << "s\n";

    FactorisationStats* stats;
    stats = GetStatsObj(algorithm_type, thread_number);

    if(manual_mode){
        manual_input(algorithm_type, thread_number, prime_table, stats);
    }
    else{
        factorise_benchmark(algorithm_type, thread_number, count, seed, prime_table, stats);
    }

    std::cout << "End of program!\n";
    return 0;
}
