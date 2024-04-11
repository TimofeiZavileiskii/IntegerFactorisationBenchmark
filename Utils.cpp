#include "Utils.h"
#include <vector>
#include <iostream>
#include "gmp.h"
#include "gmpxx.h"
#include "math.h"

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

void SieveOfEratosthenes(long upperbound, std::vector<long>& primes) {
    const int max_array_size = 1 << 29;
    long nums_passed = 0;
    long max_array_size_comp = max_array_size;

    bool first_pass = true;

    while (nums_passed < upperbound) {
        //std::cout << "new array made!" << "\n";
        long nums_left = upperbound - nums_passed;
        int array_size = (nums_left < max_array_size) ? nums_left : max_array_size;
        //std::cout << "Array size is " << array_size << "\n";
        long new_bound = nums_passed + array_size;
        
        //Set bits in the array to 1 for multiples of already identified primes
        unsigned char* sieve = new unsigned char[array_size]();
        for(long& prime : primes) {
            long array_index = (prime - (nums_passed % prime)) % prime;
            while (array_index < array_size) {
                SetBitArray(array_index, sieve);
                array_index += prime;
            }
        }
        //std::cout << "Primes considered!" << "\n";
        //Sieve through the array as normal
        for (int i = first_pass ? 2 : 0; i < array_size; i++) {
            if (!IndexBitArray(i, sieve)) {
                long new_prime = nums_passed + i;
                primes.push_back(new_prime);
                long index_num = new_prime;

                while (index_num < new_bound) {
                    long access_array = index_num - nums_passed;
                    SetBitArray(access_array, sieve);
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