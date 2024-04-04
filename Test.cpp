#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "Utils.h"
#include <gmp.h>
#include <gmpxx.h>
#include <random>
#include <iostream>

#include "GeneratePrimes.h"
#include "EllipticCurves.h"


void GetRandomState(gmp_randstate_t& random_state){
    long seed = rand();
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, seed);
}

TEST_CASE("Test Elliptic Curve Arithmetic"){
    srand(time(NULL));
    gmp_randstate_t random_state;
    GetRandomState(random_state);

    const int REPEAT_TEST = 10000;

    mpz_t mod, theta, target, mult1, mult2, diff, random_bound, gcd;
    mpz_inits(mod, theta, target, mult1, mult2, diff, random_bound, gcd, NULL);

    MontgomeryPoint point_start;
    MontgomeryPoint point1;
    MontgomeryPoint point2;
    MontgomeryPoint point_diff;
    MontgomeryPoint point_target;

    for(int i = 0; i < REPEAT_TEST; i++)
    {
        int bit_size = (rand() % 292) + 8;
        mpz_class mod_c = generate_rsa(bit_size, random_state);

        mpz_set(mod, mod_c.get_mpz_t());
        MontgomeryCurve curve(mod);

        mpz_set(random_bound, mod);
        mpz_sub_ui(random_bound, random_bound, 7);

        mpz_urandomm(theta, random_state, random_bound);
        mpz_add_ui(theta, theta, 6);
        bool inverse_not_exist = curve.ComputeC_Q(theta, point_start, random_bound); //random bound as dummy
     
        if(inverse_not_exist){  continue;   }

        point_start.Copy(point1);
        point_start.Copy(point2);
        point_start.Copy(point_diff);
        point_start.Copy(point_target);

        mpz_set(random_bound, mod);
        mpz_sub_ui(random_bound, random_bound, 5);
        mpz_urandomm(target, random_state, random_bound);
        mpz_add_ui(target, target, 3);
        mpz_sub_ui(random_bound, target, 2);
        mpz_urandomm(mult1, random_state, target);
        mpz_add_ui(mult1, mult1, 1);
        mpz_sub(mult2, target, mult1);

        mpz_sub(diff, mult2, mult1);
        if(mpz_cmp_ui(diff, 0) == 0){ continue; }
        mpz_abs(diff, diff);

        curve.MultPoints(point1, mult1);
        curve.MultPoints(point2, mult2);
        curve.MultPoints(point_diff, diff);
        curve.MultPoints(point_target, target);
        curve.AddPoints(point1, point2, point_diff);

        int inv_out1 = point_target.Reduce(mod);
        int inv_out2 = point2.Reduce(mod);

        int x_res = mpz_cmp(point_target.x, point2.x);
        int z_res = mpz_cmp(point_target.z, point2.z);

        if(inv_out1 == 0 && inv_out2 == 0){
            CHECK(x_res == 0);
            CHECK(z_res == 0);
        }
    }

    mpz_clears(mod, theta, target, mult1, mult2, diff, random_bound, gcd, NULL);
}

TEST_CASE("Prime Generation"){
    srand(time(NULL));

    const int REPEAT_TEST = 5;
    
    mpz_t prime, lower_bound, upper_bound;
    mpz_inits(prime, lower_bound, upper_bound, NULL);
    mpz_set_ui(lower_bound, 1);
    mpz_set_ui(upper_bound, 1);

    gmp_randstate_t random_state;
    GetRandomState(random_state);

    for(int i = 0; i < REPEAT_TEST; i++){
        int bit_size = rand() % 300;
        mpz_mul_2exp(upper_bound, upper_bound, bit_size);
        generate_prime(prime, lower_bound, upper_bound, random_state);
        CHECK(mpz_probab_prime_p(prime, 50) > 0);
    }

    mpz_clears(prime, lower_bound, upper_bound, NULL);
}

TEST_CASE("RSA Generation"){
    srand(time(NULL));

    const int REPEAT_TEST = 5;

    gmp_randstate_t random_state;
    GetRandomState(random_state);

    for(int i = 0; i < REPEAT_TEST; i++){
        int bit_size = (rand() % 293) + 7;
        mpz_class rsa_int = generate_rsa(bit_size, random_state).get_d();
        double rsa_integer_size = log2(generate_rsa(bit_size, random_state).get_d());
        
        CHECK(bit_size-2 < rsa_integer_size);
        CHECK(rsa_integer_size < bit_size+2);
    }
}

void CreateTestNums(mpz_t n1, mpz_t n2, mpz_t mod, MontgomeryParams& params, gmp_randstate_t random_state){
    mpz_t rem;
    mpz_init(rem);
    bool not_divisible_2 = false;
    while(!not_divisible_2){
        int bit_size = (rand() % 292) + 8;
        mpz_class mod_c = generate_rsa(bit_size, random_state);
        mpz_set(mod, mod_c.get_mpz_t());
        mpz_mod_ui(rem, mod, 2);
        if(mpz_cmp_ui(rem, 1) == 0){
            not_divisible_2 = true; 
        }
    }
    
    MontgomerySetup(mod, params);

    mpz_urandomm(n1, random_state, mod);
    mpz_urandomm(n2, random_state, mod);
    mpz_clear(rem);
}

TEST_CASE("Test Montgomery Arithmetic"){
    srand(time(NULL));
    gmp_randstate_t random_state;
    GetRandomState(random_state);

    mpz_t n1, n2, mon_n1, mon_n2, mod;
    mpz_inits(n1, n2, mon_n1, mon_n2, mod, NULL);
    
    const int REPEAT_TEST = 10000;

    //Verify Transformation In and Out
    for(int i = 0; i < REPEAT_TEST; i++){
        MontgomeryParams params;
        CreateTestNums(n1, n2, mod, params, random_state);
        mpz_set(mon_n1, n1);
        TransformIn(mon_n1, params);
        TransformOut(mon_n1, params);
        int comp_res = mpz_cmp(n1, mon_n1);
        CHECK(comp_res == 0);
    }

    //Verify Montgomery Multiplication
    for(int i = 0; i < REPEAT_TEST; i++){
        MontgomeryParams params;
        CreateTestNums(n1, n2, mod, params, random_state);

        mpz_set(mon_n1, n1);
        mpz_set(mon_n2, n2);

        mpz_mul(n1, n1, n2);
        mpz_mod(n1, n1, mod);

        TransformIn(mon_n1, params);
        TransformIn(mon_n2, params);
        MontgomeryMul(mon_n1, mon_n1, mon_n2, params);

        TransformOut(mon_n1, params);

        int comp_res = mpz_cmp(n1, mon_n1);
        CHECK(comp_res == 0);
    }

    //Verify Addition and Substraction
    for(int i = 0; i < REPEAT_TEST; i++){
        MontgomeryParams params;
        CreateTestNums(n1, n2, mod, params, random_state);

        mpz_set(mon_n1, n1);
        mpz_set(mon_n2, n2);

        TransformIn(mon_n1, params);
        TransformIn(mon_n2, params);
        if(rand() % 2 == 0){
            mpz_add(mon_n1, mon_n1, mon_n2);
            mpz_add(n1, n1, n2);
        }
        else{
            mpz_sub(mon_n1, mon_n1, mon_n2);
            mpz_sub(n1, n1, n2);
        }

        TransformOut(mon_n1, params);
        mpz_mod(mon_n1, mon_n1, mod);
        mpz_mod(n1, n1, mod);
        int comp_res = mpz_cmp(n1, mon_n1);
        CHECK(comp_res == 0);
    }

    //Verify GCD
    for(int i = 0; i < REPEAT_TEST; i++){
        MontgomeryParams params;
        CreateTestNums(n1, n2, mod, params, random_state);

        mpz_set(mon_n1, n1);

        mpz_gcd(n1, n1, mod);

        TransformIn(mon_n1, params);

        mpz_gcd(mon_n1, mon_n1, mod);

        int comp_res = mpz_cmp(n1, mon_n1);
        CHECK(comp_res == 0);
    }


    mpz_clears(n1, n2, mon_n1, mon_n2, mod, NULL);
}