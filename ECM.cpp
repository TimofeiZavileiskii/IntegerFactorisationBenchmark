#include "ECM.h"
#include "gmp.h"
#include "gmpxx.h"
#include <vector>
#include <iostream>
#include "Utils.h"
#include <cmath>
#include <random>
#include <thread>


class FactoredException: public std::exception
{
    public:
    mpz_class factor;

    FactoredException(mpz_class factor){
        this->factor = factor;
    }

    virtual const char* what() const throw()
    {
        return "Integer is factored";
    }
};


struct Point{
    public:
    mpz_class x;
    mpz_class y;
    bool is_inf;

    Point(mpz_class& x, mpz_class& y, bool is_inf){
        this->x = x;
        this->y = y;
        this->is_inf = is_inf;
    }
};

class CurveClass{
    public:

    mpz_class a;
    mpz_class b;
    mpz_class to_factor;

    CurveClass(mpz_class& a, mpz_class& b, mpz_class& to_factor){
        this->a = a;
        this->b = b;
        this->to_factor = to_factor;
    }

    Point AddPoints(Point& p1, Point& p2){
        if(p1.is_inf){
            return p2;
        }
        if(p2.is_inf){
            return p1;
        }

        if(p1.x == p2.x && p1.y == (to_factor - p2.y)){
            p1.is_inf = true;
            return p1;
        }

        mpz_class s = 1;
        if(p1.x == p2.x){
            mpz_class to_invert = (2 * p1.y) % to_factor;
            mpz_class y2;
            int inverse_found = mpz_invert(y2.get_mpz_t(), to_invert.get_mpz_t(), to_factor.get_mpz_t());
            if(inverse_found == 0){
                mpz_class factor;
                mpz_gcd(factor.get_mpz_t(), to_factor.get_mpz_t(), to_invert.get_mpz_t());
                FactoredException e(factor);
                throw e;
            }
            s = ((3 * (p1.x * p1.x) + a) * y2) % to_factor;
        }
        else{
            mpz_class to_invert = (p2.x - p1.x) % to_factor;
            mpz_class x_inv;
            int inverse_found = mpz_invert(x_inv.get_mpz_t(), to_invert.get_mpz_t(), to_factor.get_mpz_t());
            if(inverse_found == 0){
                mpz_class factor;
                mpz_gcd(factor.get_mpz_t(), to_factor.get_mpz_t(), to_invert.get_mpz_t());
                FactoredException e(factor);
                throw e;
            }
            s = (((p2.y - p1.y) % to_factor) * x_inv) % to_factor;
        }

        mpz_class x3 = (s*s - p1.x - p2.x) % to_factor;
        mpz_class y3 = (s * (p1.x - x3) - p1.y) % to_factor;
        return Point(x3, y3, false);
    }

    Point SubPoints(Point& p1, Point& p2){
        p2.y = to_factor - p2.y;
        return AddPoints(p1, p2);
    }

    Point MultPoints(Point p1, mpz_class mult){
        if(mult == 1)
            return p1;

        Point product = MultPoints(p1, mult / 2);
        Point square = AddPoints(product, product);
        if((mult & 1) == 1)
            return AddPoints(p1, square);
        else
            return square;
    }
};


void EcmJob(mpz_t& output, mpz_t& to_factor, volatile bool& factored, std::vector<mpz_class>& primes, long starting_seed){
    double to_factor_d = mpz_get_d(to_factor);
    to_factor_d = std::pow(to_factor_d, 1.0/3.0);
    double bound_log = std::log2(to_factor_d);
    mpz_class B = to_factor_d; 

    mpz_class n = mpz_class(to_factor);
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, starting_seed);
    
    while(!factored)
    {
        mpz_t a_s, x1_s, y1_s, random_bound;
        mpz_inits(a_s, x1_s, y1_s, random_bound, NULL);
        mpz_set_ui(random_bound, 50000);
        mpz_urandomm(a_s, random_state, random_bound);
        mpz_urandomm(x1_s, random_state, random_bound);
        mpz_urandomm(y1_s, random_state, random_bound);

        mpz_class a = mpz_class(a_s) % n;
        mpz_class x1 = mpz_class(x1_s) % n;
        mpz_class y1 = mpz_class(y1_s) % n;
        
        mpz_class b = (y1*y1 - x1*x1*x1 - a*x1) % n;

        mpz_t g, pow_disc;
        mpz_inits(g, pow_disc, NULL);
        mpz_class to_factor_c = mpz_class(to_factor);
        mpz_powm_ui(pow_disc, a.get_mpz_t(), 3, to_factor);
        
        mpz_class pow_disc_c = mpz_class(pow_disc);
        mpz_class disc = 4 * pow_disc_c + 27 * (b*b);
        mpz_gcd(g, disc.get_mpz_t(), to_factor);

        if(mpz_cmp(g, to_factor) == 0)
            continue;
        
        if(mpz_cmp_ui(g, 1) > 0){
            mpz_set(output, g);
            factored = true;
        }

        CurveClass curve = CurveClass(a, b, n);
        Point point = Point(x1, y1, false);
        try{
            for(mpz_class& prime : primes){
                double prime_log;
                double p_d = prime.get_d();
                long e = (long)(bound_log/std::log2(p_d));
                mpz_t f_s;
                mpz_init(f_s);
                mpz_pow_ui(f_s, prime.get_mpz_t(), e);
                mpz_class f = mpz_class(f_s);
                point = curve.MultPoints(point, f);
                if(factored){
                    break;
                }
            }
        }
        catch (FactoredException e){
            mpz_set(output, e.factor.get_mpz_t());
            factored = true;
        }
    }
}


void Ecm(mpz_t& output, mpz_t& to_factor, int thread_count){
    double to_factor_d = mpz_get_d(to_factor);
    to_factor_d = std::pow(to_factor_d, 1.0/5.0);
    double bound_log = std::log2(to_factor_d);
    mpz_class B = to_factor_d;
    std::cout << "Selected bound is ";
    mpz_out_str(NULL, 10, B.get_mpz_t());
    std::cout << std::endl;
    std::vector<mpz_class> primes;
    SieveOfEratosthenes(B, primes);

    volatile bool factored = false;

    if(thread_count == 1){
        long seed = rand();
        EcmJob(output, to_factor, factored, primes, seed);
    }
    else{
        std::vector<std::thread> workers = std::vector<std::thread>(); 

        for(int i = 0; i < thread_count; i++){
            long random_seed = rand();
            workers.emplace_back(std::thread(EcmJob, std::ref(output), std::ref(to_factor), std::ref(factored), std::ref(primes), random_seed));
        }
        for(std::thread &thread : workers){
            thread.join();
        }
    }
}
    
    

    