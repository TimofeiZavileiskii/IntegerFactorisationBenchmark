#include "ECM.h"
#include "gmp.h"
#include "gmpxx.h"
#include <vector>
#include <iostream>
#include "Utils.h"
#include <cmath>
#include <random>
#include <thread>


struct Point{
    mpz_t x;
    mpz_t z;

    inline Point(){
        mpz_inits(x, z, NULL);
    }

    inline ~Point(){
        mpz_clears(x, z, NULL);
    }

    inline void Copy(Point& point){
        mpz_set(point.x, x);
        mpz_set(point.z, z);
    }
};

class Curve{
    private:
    Point p_copy;
    Point p_double;
    mpz_t x_sqr, z_sqr, x_z, x_z_sqr, x2z1;

    public:
    mpz_t c, mod, four; //Four is in montgomery space
    MontgomeryParams params;

    Curve(mpz_t n){
        mpz_inits(c, mod, x_sqr, z_sqr, x_z, x2z1, four, NULL);
        mpz_set(mod, n);
        MontgomerySetup(mod, params);
        mpz_set_ui(four, 4);
        TransformIn(four, params);
    }

    ~Curve(){
        mpz_clears(c, mod, x_sqr, z_sqr, x_z, x2z1, four, NULL);
    }

    void DoublePoint(Point& point){

        MontgomeryMul(x_sqr, point.x, point.x, params);
        MontgomeryMul(z_sqr, point.z, point.z, params);
        MontgomeryMul(x_z, point.x, point.z, params);
        MontgomeryMul(point.z, x_z, c, params);
        mpz_add(point.z, point.z, x_sqr);
        mpz_add(point.z, point.z, z_sqr);

        MontgomeryMul(point.z, point.z, x_z, params);
        MontgomeryMul(point.z, point.z, four, params);

        mpz_sub(x_sqr, x_sqr, z_sqr);
        MontgomeryMul(point.x, x_sqr, x_sqr, params);
    }

    void AddPoints(Point& point, Point& to_add, Point& p_min_add){
        #define x1x2 x_sqr
        #define z1z2 z_sqr
        #define x1z2 x_z

        MontgomeryMul(x1x2, point.x, to_add.x, params);
        MontgomeryMul(z1z2, point.z, to_add.z, params);
        MontgomeryMul(x1z2, point.x, to_add.z, params);
        MontgomeryMul(x2z1, to_add.x, point.z, params);

        ModularSub(to_add.x, x1x2, z1z2, mod);
        MontgomeryMul(to_add.x, to_add.x, to_add.x, params);
        MontgomeryMul(to_add.x, to_add.x, p_min_add.z, params);

        ModularSub(to_add.z, x1z2, x2z1, mod);
        MontgomeryMul(to_add.z, to_add.z, to_add.z, params);
        MontgomeryMul(to_add.z, to_add.z, p_min_add.x, params);
        #undef x1x2
        #undef z1z2
        #undef x1z2
    }

    void MultPoints(Point& point, mpz_t multiple, mpz_t factor, volatile bool& is_factored){
        if(mpz_cmp_ui(multiple, 1) == 0){
            return;
        }
        if(mpz_cmp_ui(multiple, 0) == 0){
            mpz_set_ui(point.x, 0);
            mpz_set_ui(point.z, 0);
            return;
        }
        if(mpz_cmp_ui(multiple, 2) == 0){
            DoublePoint(point);
            return;
        }
        
        point.Copy(p_copy);
        point.Copy(p_double);
        DoublePoint(p_double);
        
        int bit_count = (int)log2(mpz_get_d(multiple))+1;

        for(int i = bit_count-2; (i >= 0) && (!is_factored); i--){
            if(mpz_tstbit(multiple, i) == 1){
                AddPoints(p_double, p_copy, point);
                DoublePoint(p_double);
            }
            else{
                AddPoints(p_copy, p_double, point);
                DoublePoint(p_copy);
            }
        }

        if(mpz_tstbit(multiple, 0) == 1){
            AddPoints(p_copy, p_double, point);
            p_double.Copy(point);
        }
        else{
            DoublePoint(p_copy);
            p_copy.Copy(point);
        }
    }
};


inline bool ComputeC_Q(mpz_t theta, mpz_t u, mpz_t v, mpz_t diff_cb, mpz_t u3_v, mpz_t denom, mpz_t mod, Curve& curve, Point& point, mpz_t factor){
    mpz_mul(u, theta, theta);
    mpz_sub_ui(u, u, 5);
    mpz_mod(u, u, mod);
    mpz_mul_2exp(v, theta, 2);
    mpz_mod(v, v, mod);

    ModularSub(diff_cb, v, u, mod);
    mpz_powm_ui(diff_cb, diff_cb, 3, mod);
    mpz_mul_ui(u3_v, u, 3);
    mpz_add_ui(u3_v, u3_v, 3);
    mpz_mul(diff_cb, diff_cb, u3_v);
    mpz_mod(diff_cb, diff_cb, mod);


    mpz_powm_ui(denom, u, 3, mod);
    mpz_mul_ui(denom, denom, 3);
    mpz_mul(denom, denom, v);
    mpz_mod(denom, denom, mod);

    int inverse_exists = mpz_invert(u3_v, denom, mod);
    if(inverse_exists == 0){
        mpz_gcd(factor, denom, mod);
        return true;
    }
    
    mpz_mul(diff_cb, u3_v, diff_cb);
    mpz_sub_ui(diff_cb, diff_cb, 2);
    mpz_mod(curve.c, diff_cb, mod);

    mpz_powm_ui(point.x, u, 3, mod);
    mpz_powm_ui(point.z, v, 3, mod);

    return false;
}

inline double log_base(double x, double base) {
    return log(x) / log(base);
}

void EcmJob(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound)
    {
    mpz_class n = mpz_class(to_factor);
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, starting_seed);
    
    Curve curve(to_factor);
    mpz_t theta, u, v, random_bound, diff_cb, u3_v, denom, f, gcd;
    mpz_inits(theta, u, v, random_bound, diff_cb, u3_v, denom, f, gcd, NULL);

    Point point;
    mpz_set(random_bound, to_factor);
    mpz_sub_ui(random_bound, random_bound, 7);
    while(!factored)
    {
        //Find first point on the curve
        //Randomly select a, x and y
        mpz_urandomm(theta, random_state, to_factor);
        mpz_add_ui(theta, theta, 6);

        //Compute b value for the curve based on the random coordinates
        bool inverse_not_exist = ComputeC_Q(theta, u, v, diff_cb, u3_v, denom, to_factor, curve, point, output);
        if(inverse_not_exist){
            if(mpz_cmp(output, to_factor) == 0){
                continue;
            }
            else{
                return;
            }
        }

        for(mpz_class& prime : primes){
            if(factored)
                break;

            long exponent = log_base(bound, prime.get_d());
            mpz_pow_ui(f, prime.get_mpz_t(), exponent);
            curve.MultPoints(point, f, output, factored);
            mpz_gcd(gcd, point.z, to_factor);
            if(mpz_cmp_ui(point.z, 0) == 0){
                std::cout << "Z IS 0" << std::endl;
                throw 0;
            }
            if(mpz_cmp_ui(gcd, 1) > 0){
                mpz_set(output, gcd);
                factored = true;
            }
        }
    }
    mpz_clears(theta, u, v, random_bound, diff_cb, u3_v, denom, f, gcd, NULL);
}

void Ecm(mpz_t& output, mpz_t& to_factor, int thread_count){
    double to_factor_d = mpz_get_d(to_factor);
    double smallest_factor = sqrt(to_factor_d);
    double log_smallest_factor = log(smallest_factor);
    double bound = std::exp(sqrt(2.0)/2.0 * sqrt(log(log_smallest_factor)*log_smallest_factor));
    mpz_class B = bound;
    std::vector<mpz_class> primes;
    SieveOfEratosthenes(B, primes);

    volatile bool factored = false;

    if(thread_count == 1){
        long seed = rand();
        EcmJob(output, to_factor, factored, primes, seed, B.get_d());
    }
    else{
        std::vector<std::thread> workers = std::vector<std::thread>(); 

        for(int i = 0; i < thread_count; i++){
            long random_seed = rand();
            workers.emplace_back(std::thread(EcmJob, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d()));
        }
        for(std::thread &thread : workers){
            thread.join();
        }
    }
}
    
    

    