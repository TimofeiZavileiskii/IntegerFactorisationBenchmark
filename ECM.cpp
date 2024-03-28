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
    mpz_t y;
    bool is_inf;

    Point(){
        mpz_inits(x, y, NULL);
        is_inf = false;
    }

    ~Point(){
        mpz_clears(x, y, NULL);
    }
};

class Curve{
    private:
    Point p_copy;
    mpz_t neg_y, s, to_invert, y2, x_inv;
    
    public:
    mpz_t a, b, n;

    Curve(){
        mpz_inits(a, b, n, neg_y, s, to_invert, y2, x_inv, NULL);
    }

    ~Curve(){
        mpz_clears(a, b, n, to_invert, y2, x_inv, NULL);
    }

    void AddPoints(Point& point, Point& to_add, mpz_t factor, volatile bool& is_factored){
        if(to_add.is_inf)
            return;
        
        if(point.is_inf)
        {
            mpz_set(point.x, to_add.x);
            mpz_set(point.y, to_add.y);
            return;
        }

        mpz_sub(neg_y, n, to_add.y);
        bool x_equal = mpz_cmp(point.x, to_add.x) == 0;
        if(x_equal && (mpz_cmp(point.y, neg_y) == 0))
        {
            point.is_inf = true;
            return;
        }

        mpz_set_ui(to_invert, 1);
        if(x_equal){
            mpz_mul_ui(to_invert, point.y, 2);
            mpz_mod(to_invert, to_invert, n);
            int inverse_found = mpz_invert(y2, to_invert, n);
            if(inverse_found == 0){
                mpz_gcd(factor, n, to_invert);
                is_factored = true;
                return;
            }
            mpz_mul(s, point.x, point.x);
            mpz_mul_ui(s, s, 3);
            mpz_add(s, s, a);
            mpz_mul(s, s, y2);
            mpz_mod(s, s, n);
        }
        else{
            mpz_sub(to_invert, to_add.x, point.x);
            mpz_mod(to_invert, to_invert, n);
            int inverse_found = mpz_invert(x_inv, to_invert, n);
            if(inverse_found == 0){
                mpz_gcd(factor, n, to_invert);
                is_factored = true;
                return;
            }
            mpz_sub(s, to_add.y, point.y);
            mpz_mod(s, s, n);
            mpz_mul(s, s, x_inv);
            mpz_mod(s, s, n);
        }

        //to_invert = x3, new x coordinate
        mpz_mul(to_invert, s, s);
        mpz_sub(to_invert, to_invert, point.x);
        mpz_sub(to_invert, to_invert, to_add.x);
        mpz_mod(to_invert, to_invert, n);

        //neg_y = y3, new y coordinate
        mpz_sub(neg_y, point.x, to_invert);
        mpz_mul(neg_y, neg_y, s);
        mpz_sub(neg_y, neg_y, point.y);
        mpz_mod(point.y, neg_y, n);
        mpz_set(point.x, to_invert);
        point.is_inf = false;
    }

    void MultPointsRecursive(Point& point_copy, Point& point, mpz_t multiple, mpz_t factor, volatile bool& is_factored){
        if(mpz_cmp_ui(multiple, 0) == 0 || mpz_cmp_ui(multiple, 1) == 0)
            return;

        bool is_odd = mpz_odd_p(multiple) != 0;
        mpz_div_ui(multiple, multiple, 2);
        
        MultPointsRecursive(point_copy, point, multiple, factor, is_factored);
        if(is_factored)
            return;

        AddPoints(point, point, factor, is_factored);
        if(is_odd && !is_factored){
            AddPoints(point, point_copy, factor, is_factored);}
    }

    inline void MultPoints(Point& point, mpz_t multiple, mpz_t factor, volatile bool& is_factored){
        mpz_set(p_copy.x, point.x);
        mpz_set(p_copy.y, point.y);
        
        MultPointsRecursive(p_copy, point, multiple, factor, is_factored);
    }
};


inline double log_base(double x, double base) {
    return log(x) / log(base);
}

void EcmJob(mpz_t output, mpz_t to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed, 
    double bound)
    {
    mpz_class n = mpz_class(to_factor);
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, starting_seed);
    
    Curve curve;
    mpz_set(curve.n, to_factor);
    mpz_t x_cube, a_x, random_bound, g, pow_disc, b_sqr, disc, f;
    mpz_inits(x_cube, a_x, random_bound, g, pow_disc, b_sqr, disc, f, NULL);

    Point point = Point();
    
    while(!factored)
    {
        //Find first point on the curve
        //Randomly select a, x and y
        mpz_urandomm(curve.a, random_state, to_factor);
        mpz_urandomm(point.x, random_state, to_factor);
        mpz_urandomm(point.y, random_state, to_factor);
        
        //Compute b value for the curve based on the random coordinates
        mpz_mul(curve.b, point.y, point.y);
        mpz_mul(x_cube, point.x, point.x);
        mpz_mul(x_cube, x_cube, point.x);
        mpz_mul(a_x, point.x, curve.a);
        mpz_sub(curve.b, curve.b, x_cube);
        mpz_sub(curve.b, curve.b, a_x);
        mpz_mod(curve.b, curve.b, to_factor);

        //Check discriminant 
        mpz_pow_ui(pow_disc, curve.a, 3);
        mpz_mul(b_sqr, curve.b, curve.b);
        mpz_mul_ui(b_sqr, b_sqr, 27);
        mpz_mul_ui(pow_disc, pow_disc, 4);
        mpz_add(pow_disc, pow_disc, b_sqr);
        mpz_gcd(g, pow_disc, to_factor);

        if(mpz_cmp(g, to_factor) == 0)
            continue;
        
        if(mpz_cmp_ui(g, 1) > 0){
            mpz_set(output, g);
            factored = true;
            continue;
        }
        
        for(mpz_class& prime : primes){
            if(factored)
                break;

            long exponent = log_base(bound, prime.get_d());
            mpz_pow_ui(f, prime.get_mpz_t(), exponent);
            curve.MultPoints(point, f, output, factored);
        }
        
    }
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
    
    

    