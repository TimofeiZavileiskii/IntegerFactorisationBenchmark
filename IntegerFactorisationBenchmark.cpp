#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include "gmp.h"
#include "Utils.h"
#include "PollardsRho.h"
#include "PollardsP1.h"
#include "TrialDivision.h"
#include "ECM.h"

void manual_input(){
    mpz_t to_factor, divisor, other_divisor, check;
    mpz_inits(to_factor, divisor, other_divisor, check, NULL);

    std::cout << "Enter number to factor:" << std::endl;
    mpz_inp_str(to_factor, NULL, 10);

    PollardsRho(divisor, to_factor, 1);

    mpz_div(other_divisor, to_factor, divisor);
    mpz_mod(check, to_factor, divisor);
    if(mpz_cmp_ui(check, 0) == 0){
        std::cout << "Success!" << std::endl;
        std::cout << "The factorisation of ";
        mpz_out_str(NULL, 10, to_factor);
        std::cout << ": ";
        mpz_out_str(NULL, 10, divisor);
        std::cout << " * ";
        mpz_out_str(NULL, 10, other_divisor);
        std::cout << std::endl;
    }
    else{
        std::cout << "Something went wrong\n";
        mpz_out_str(NULL, 10, divisor);
        std::cout << " does not divide ";
        mpz_out_str(NULL, 10, to_factor);
        std::cout << "\n";
    }
}

float benchmark_number(mpz_t& rsa_num){
    mpz_t factor, other_factor;
    mpz_inits(factor, other_factor, NULL);

    std::cout << "Number to factor: ";
    mpz_out_str(NULL, 10, rsa_num);
    std::cout << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    Ecm(factor, rsa_num, 8);
    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<float> time = end_time - start_time;

    mpz_div(other_factor, rsa_num, factor);
    std::cout << "Number factorised: ";
    mpz_out_str(NULL, 10, rsa_num);
    std::cout << " = ";
    mpz_out_str(NULL, 10, factor);
    std::cout << " * ";
    mpz_out_str(NULL, 10, other_factor);
    std::cout << std::endl;
    std::cout << "Time spent: " << time.count() << " s" << std::endl;

    return time.count();
}

void write_benchmark_result(std::vector<float>& benchmark_times, int bit_size, std::string& filename){
    std::fstream benchmark_result;
    benchmark_result.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
    if (!benchmark_result ){
        benchmark_result.open(filename,  std::fstream::in | std::fstream::out | std::fstream::trunc);
    } 
    benchmark_result << bit_size;
    for(float time : benchmark_times){
        benchmark_result << "," << time;
    }
    benchmark_result << "\n";
    benchmark_result.close();
}


void factorise_benchmark(){
    std::vector<int> bit_sizes;
    std::vector<float> times;
    std::string benchmark_name = "rsa_numbers.csv";

    std::ifstream benchmark(benchmark_name);
    
    if(!benchmark.is_open()){
        throw std::runtime_error("Benchmark file failed to open");
    }

    std::string line;
    int count = 0;
    int max_count = 60;
    while(std::getline(benchmark, line))
    {

        std::vector<float> times_for_bitsize;

        std::stringstream ss(line);
        
        std::string substr;
        std::getline(ss, substr, ',');
        int bit_size = stoi(substr);
        bit_sizes.push_back(bit_size);

        while (ss.good()) {
            std::string substr;
            std::getline(ss, substr, ',');
            mpz_t rsa_num;
            mpz_init(rsa_num);
            mpz_set_str(rsa_num, substr.c_str(), 10);
            times_for_bitsize.push_back(benchmark_number(rsa_num));
            count++;
            if(count > max_count){
                break;
            }
        }

        float average = 0;
        for(float time : times_for_bitsize){
            average += time;
        }
        average /= times_for_bitsize.size();
        times.push_back(average);

        std::cout << "------------------------- Bit Size: " << bit_size << " Average Time: " << average << std::endl;
        std::string results_path = "benchmark_results.csv";
        write_benchmark_result(times_for_bitsize, bit_size, results_path);
        if(count > max_count){
            break;
        }
    }
}

int main(){
    std::cout << "Program Start! " << std::endl;
    factorise_benchmark();
    
    std::cout << "End of program!";
    return 0;
}