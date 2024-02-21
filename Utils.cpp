#include "Utils.h"
#include <vector>
#include <iostream>
#include "gmp.h"
#include "gmpxx.h"


bool inline IndexBitArray(long index, unsigned char* array) {
    unsigned char bit = array[index >> 3];
    return (bool)((bit >> (index & 0x7)) & 1);
}

void inline SetBitArray(long index, unsigned char* array) {
    unsigned char& bit = array[index >> 3];
    bit |= 1 << (index & 0x7);
}

void inline ZeroBitArray(long index, unsigned char* array) {
    unsigned char& bit = array[index >> 3];
    bit &= 1 << (index & 0x7);
}

void SieveOfEratosthenes(mpz_class upperbound, std::vector<mpz_class>& primes) {
    const int max_array_size = 1 << 29;
    mpz_class nums_passed = 0;
    mpz_class max_array_size_comp = max_array_size;

    bool first_pass = true;

    while (nums_passed < upperbound) {
        //std::cout << "new array made!" << "\n";
        mpz_class nums_left = upperbound - nums_passed;
        int array_size = (nums_left < max_array_size) ? nums_left.get_si() : max_array_size;
        //std::cout << "Array size is " << array_size << "\n";
        mpz_class new_bound = nums_passed + array_size;
        
        //Set bits in the array to 1 for multiples of already identified primes
        unsigned char* sieve = new unsigned char[array_size]();
        for(mpz_class& prime : primes) {
            mpz_class array_index = (prime - (nums_passed % prime)) % prime;
            while (array_index < array_size) {
                SetBitArray(array_index.get_si(), sieve);
                array_index += prime;
            }
        }
        //std::cout << "Primes considered!" << "\n";
        //Sieve through the array as normal
        for (int i = first_pass ? 2 : 0; i < array_size; i++) {
            if (!IndexBitArray(i, sieve)) {
                mpz_class new_prime = nums_passed + i;
                primes.push_back(new_prime);
                mpz_class index_num = new_prime;

                while (index_num < new_bound) {
                    mpz_class access_array = index_num - nums_passed;
                    SetBitArray(access_array.get_si(), sieve);
                    index_num += new_prime;
                }
            }
        }

        //std::cout << "new bound is " << new_bound << "\n";
        nums_passed = new_bound;
        first_pass = false;
        delete[] sieve;
    }
}