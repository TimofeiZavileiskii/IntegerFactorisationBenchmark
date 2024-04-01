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
    mpz_t x_sub_z;
    mpz_t x_add_z;

    inline Point(){
        mpz_inits(x, z, x_sub_z, x_add_z, NULL);
    }

    inline ~Point(){
        mpz_clears(x, z, x_sub_z, x_add_z, NULL);
    }

    inline void ComputeDiff(mpz_t mod){
        ModularSub(x_sub_z, x, z, mod);
        ModularAdd(x_add_z, x, z, mod);
    }

    inline void Copy(Point& point){
        mpz_set(point.x, x);
        mpz_set(point.z, z);
        mpz_set(point.x_sub_z, x_sub_z);
        mpz_set(point.x_add_z, x_add_z);
    }
};

class Curve{
    private:
    Point p_copy;
    Point p_double;

    mpz_t a_2_over_4, cross_1, sqr1, sqr2;

    public:
    mpz_t c; //Four is in montgomery space
    MontgomeryParams params;

    Curve(mpz_t n){
        mpz_inits(c, a_2_over_4, cross_1, sqr1, sqr2, NULL);
        MontgomerySetup(n, params);
    }

    ~Curve(){
        mpz_clears(c, a_2_over_4, cross_1, sqr1, sqr2, NULL);
    }

    void inline SetC(mpz_t c){
        mpz_set(this->c, c);
        mpz_add_ui(a_2_over_4, c, 2);
        mpz_set_ui(cross_1, 4);
        mpz_invert(cross_1, cross_1, params.mod);
        mpz_mul(a_2_over_4, a_2_over_4, cross_1);
        mpz_mod(a_2_over_4, a_2_over_4, params.mod);

        TransformIn(this->c, params);
        TransformIn(a_2_over_4, params);
    }

    void DoublePoint(Point& point){
        MontgomeryMul(sqr1, point.x_add_z, point.x_add_z, params);
        MontgomeryMul(sqr2, point.x_sub_z, point.x_sub_z, params);
        ModularSub(cross_1, sqr1, sqr2, params.mod);
        MontgomeryMul(point.x, sqr1, sqr2, params);
        MontgomeryMul(point.z, a_2_over_4, cross_1, params);
        ModularAdd(point.z, point.z, sqr2, params.mod);
        MontgomeryMul(point.z, point.z, cross_1, params);
        point.ComputeDiff(params.mod);
    }

    void AddPoints(Point& point, Point& to_add, Point& p_min_add){
        #define cross_2 sqr1
        MontgomeryMul(cross_1, point.x_sub_z, to_add.x_add_z, params);
        MontgomeryMul(cross_2, point.x_add_z, to_add.x_sub_z, params);
        ModularAdd(to_add.x, cross_1, cross_2, params.mod);
        MontgomeryMul(to_add.x, to_add.x, to_add.x, params);
        MontgomeryMul(to_add.x, to_add.x, p_min_add.z, params);

        mpz_sub(to_add.z, cross_1, cross_2);
        MontgomeryMul(to_add.z, to_add.z, to_add.z, params);
        MontgomeryMul(to_add.z, to_add.z, p_min_add.x, params);
        to_add.ComputeDiff(params.mod);
        #undef cross_2
    }

    void MultPoints(Point& point, mpz_t multiple, mpz_t factor, volatile bool& is_factored){
        if(mpz_cmp_ui(multiple, 0) == 0){
            mpz_set_ui(point.x, 0);
            mpz_set_ui(point.z, 0);
            return;
        }
        if(mpz_cmp_ui(multiple, 1) == 0){
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

        for(int i = bit_count-2; i > 0; i--){
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
            AddPoints(p_double, p_copy, point);
            p_copy.Copy(point);
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
    mpz_add(u3_v, u3_v, v);
    mpz_mul(diff_cb, diff_cb, u3_v);
    mpz_mod(diff_cb, diff_cb, mod);

    mpz_powm_ui(denom, u, 3, mod);
    mpz_mul_2exp(denom, denom, 2);
    mpz_mul(denom, denom, v);
    mpz_mod(denom, denom, mod);

    int inverse_exists = mpz_invert(u3_v, denom, mod);
    if(inverse_exists == 0){
        mpz_gcd(factor, denom, mod);
        return true;
    }
    
    mpz_mul(diff_cb, u3_v, diff_cb);
    mpz_sub_ui(diff_cb, diff_cb, 2);
    mpz_mod(diff_cb, diff_cb, mod);
    curve.SetC(diff_cb);

    mpz_powm_ui(point.x, u, 3, mod);
    TransformIn(point.x, curve.params);
    mpz_powm_ui(point.z, v, 3, mod);
    TransformIn(point.z, curve.params);
    point.ComputeDiff(mod);
    return false;
}

inline double log_base(double x, double base) {
    return log(x) / log(base);
}

void EcmJob(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, double bound, int& curves_tried)
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
        curves_tried += 1;
        //Randomly select theta
        mpz_urandomm(theta, random_state, random_bound);
        mpz_add_ui(theta, theta, 6);

        //Generate x, z, C from theta
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
        }
        mpz_gcd(gcd, point.z, to_factor);
        if(mpz_cmp_ui(point.z, 0) == 0){
            std::cout << "Z IS 0" << std::endl;
            continue;
        }
        if(mpz_cmp_ui(gcd, 1) > 0){
            mpz_set(output, gcd);
            factored = true;
        }
}
    mpz_clears(theta, u, v, random_bound, diff_cb, u3_v, denom, f, gcd, NULL);
}

void Ecm(mpz_t& output, mpz_t& to_factor, int thread_count){
    double to_factor_d = mpz_get_d(to_factor);
    double smallest_factor = sqrt(to_factor_d);
    double log_smallest_factor = log(smallest_factor);
    double bound = std::exp((sqrt(2.0))/2.0 * sqrt(log(log_smallest_factor)*log_smallest_factor));
    mpz_class B = bound;
    if(B%2 == 1){
        B -= 1;
    }
    std::vector<mpz_class> primes;
    SieveOfEratosthenes(B, primes);

    volatile bool factored = false;
    int sum = 0;
    if(thread_count == 1){
        long seed = rand();
        EcmJob(output, to_factor, factored, primes, seed, B.get_d(), sum);
    }
    else{
        std::vector<std::thread> workers = std::vector<std::thread>(); 
        int* curve_nums = new int[thread_count];
        for(int i = 0; i < thread_count; i++){
            long random_seed = rand();
            curve_nums[i] = 0;
            workers.emplace_back(std::thread(EcmJob, output, to_factor, std::ref(factored), std::ref(primes), random_seed, B.get_d(), std::ref(curve_nums[i])));
        }
        for(std::thread &thread : workers){
            thread.join();
        }

        for(int i = 0; i < thread_count; i++){
            sum += curve_nums[i];
        }
        delete[] curve_nums;
    }
    std::cout << "Total curve number: " << sum << "\n";
}
    
void EcmTest(mpz_t& output, mpz_t& to_factor){
    mpz_class n = mpz_class(to_factor);
    long starting_seed = rand();
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, starting_seed);
    
    Curve curve(to_factor);
    mpz_t theta, u, v, random_bound, diff_cb, u3_v, denom, f, gcd;
    mpz_inits(theta, u, v, random_bound, diff_cb, u3_v, denom, f, gcd, NULL);

    Point point;
    mpz_set(random_bound, to_factor);
    mpz_sub_ui(random_bound, random_bound, 7);

    mpz_urandomm(theta, random_state, random_bound);
    mpz_add_ui(theta, theta, 6);

    bool factored = false;
    mpz_t multiple, inv;
    mpz_inits(multiple, inv, NULL);
    Point p_8;
    Point p_2_copy;
    Point p_6;
    Point p_2;

    //Generate x, z, C from theta
    bool inverse_not_exist = ComputeC_Q(theta, u, v, diff_cb, u3_v, denom, to_factor, curve, point, output);
    point.Copy(p_8);
    curve.DoublePoint(p_8);
    curve.DoublePoint(p_8);
    curve.DoublePoint(p_8);

    point.Copy(p_6);
    point.Copy(p_2);
    curve.DoublePoint(p_2);
    p_2.Copy(p_2_copy);
    curve.DoublePoint(p_2);
    mpz_set_ui(multiple, 6);
    curve.MultPoints(p_6, multiple, output, factored);
    mpz_set_ui(multiple, 4);
    curve.AddPoints(p_6, p_2_copy, p_2);
    
    //curve.DoublePoint(p_2_5);
    mpz_invert(inv, p_8.z, to_factor);
    mpz_mul(p_8.z, p_8.z, inv);
    mpz_mod(p_8.z, p_8.z, to_factor);
    mpz_mul(p_8.x, p_8.x, inv);
    mpz_mod(p_8.x, p_8.x, to_factor);

    mpz_invert(inv, p_2_copy.z, to_factor);
    mpz_mul(p_2_copy.z, p_2_copy.z, inv);
    mpz_mod(p_2_copy.z, p_2_copy.z, to_factor);
    mpz_mul(p_2_copy.x, p_2_copy.x, inv);
    mpz_mod(p_2_copy.x, p_2_copy.x, to_factor);

    std::cout << "Point values:\nP_3 x:";
    mpz_out_str(NULL, 10, p_8.x);
    std::cout << " z: ";
    mpz_out_str(NULL, 10, p_8.z);
    std::cout << "\nM_3 x: ";
    mpz_out_str(NULL, 10, p_2_copy.x);
    std::cout << " z: ";
    mpz_out_str(NULL, 10, p_2_copy.z);
    std::cout << "\n";

    throw 0;
}