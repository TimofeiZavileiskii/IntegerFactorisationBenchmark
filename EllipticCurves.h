#include "gmp.h"
#include "gmpxx.h"
#include "Utils.h"


struct MontgomeryPoint{
    mpz_t x;
    mpz_t z;
    mpz_t x_sub_z;
    mpz_t x_add_z;

    inline MontgomeryPoint(){
        mpz_inits(x, z, x_sub_z, x_add_z, NULL);
    }

    inline ~MontgomeryPoint(){
        mpz_clears(x, z, x_sub_z, x_add_z, NULL);
    }

    inline void ComputeDiff(mpz_t mod){
        ModularSub(x_sub_z, x, z, mod);
        ModularAdd(x_add_z, x, z, mod);
    }

    inline void ComputeXZ(mpz_t xz, mpz_t mod){
        mpz_mul(xz, x, z);
        mpz_mod(xz, xz, mod);
    }

    inline void Copy(MontgomeryPoint& point){
        mpz_set(point.x, x);
        mpz_set(point.z, z);
        mpz_set(point.x_sub_z, x_sub_z);
        mpz_set(point.x_add_z, x_add_z);
    }

    inline int Reduce(mpz_t mod){
        mpz_t inv;
        mpz_init(inv);
        
        int inv_out = mpz_invert(inv, z, mod);
        mpz_mul(z, z, inv);
        mpz_mod(z, z, mod);

        mpz_mul(x, x, inv);
        mpz_mod(x, x, mod);

        mpz_clear(inv);
        return inv_out;
    }
};


inline bool test_bit(long n, int bit){
    return (n & (1 << (bit))) != 0;    
}


class MontgomeryCurve{
    private:
    MontgomeryPoint p_copy;
    MontgomeryPoint p_double;
    mpz_t cross_1, sqr1, sqr2;
    
    public:
    //Four is in montgomery space
    mpz_t a_2_over_4, c;
    MontgomeryParams params;

    inline MontgomeryCurve(mpz_t n){
        mpz_inits(c, a_2_over_4, cross_1, sqr1, sqr2, NULL);
        MontgomerySetup(n, params);
    }

    inline ~MontgomeryCurve(){
        mpz_clears(c, a_2_over_4, cross_1, sqr1, sqr2, NULL);
    }

    inline void SetC(mpz_t c){
        mpz_set(this->c, c);
        mpz_add_ui(a_2_over_4, c, 2);
        mpz_set_ui(cross_1, 4);
        mpz_invert(cross_1, cross_1, params.mod);
        mpz_mul(a_2_over_4, a_2_over_4, cross_1);
        mpz_mod(a_2_over_4, a_2_over_4, params.mod);
        
        TransformIn(this->c, params);
        TransformIn(a_2_over_4, params);
    }

    inline bool ComputeC_Q(mpz_t theta, MontgomeryPoint& point, mpz_t factor){
            #define diff_cb cross_1
            #define u3_v sqr1
            #define denom sqr2
            #define u point.x
            #define v point.z

            mpz_mul(u, theta, theta);
            mpz_sub_ui(u, u, 5);
            mpz_mod(u, u, params.mod);
            mpz_mul_2exp(v, theta, 2);
            mpz_mod(v, v, params.mod);

            ModularSub(diff_cb, v, u, params.mod);
            mpz_powm_ui(diff_cb, diff_cb, 3, params.mod);
            mpz_mul_ui(u3_v, u, 3);
            mpz_add(u3_v, u3_v, v);
            mpz_mul(diff_cb, diff_cb, u3_v);
            mpz_mod(diff_cb, diff_cb, params.mod);

            mpz_powm_ui(denom, u, 3, params.mod);
            mpz_mul_2exp(denom, denom, 2);
            mpz_mul(denom, denom, v);
            mpz_mod(denom, denom, params.mod);

            int inverse_exists = mpz_invert(u3_v, denom, params.mod);
            if(inverse_exists == 0){
                mpz_gcd(factor, denom, params.mod);
                return true;
            }
            
            mpz_mul(diff_cb, u3_v, diff_cb);
            mpz_sub_ui(diff_cb, diff_cb, 2);
            mpz_mod(diff_cb, diff_cb, params.mod);
            SetC(diff_cb);

            mpz_powm_ui(point.x, u, 3, params.mod);
            TransformIn(point.x, params);
            mpz_powm_ui(point.z, v, 3, params.mod);
            TransformIn(point.z, params);
            point.ComputeDiff(params.mod);
            return false;

            #undef diff_cb
            #undef u3_v
            #undef denom
            #undef u
            #undef v
    }

    inline void DoublePoint(MontgomeryPoint& point){
        MontgomeryMul(sqr1, point.x_add_z, point.x_add_z, params);
        MontgomeryMul(sqr2, point.x_sub_z, point.x_sub_z, params);
        ModularSub(cross_1, sqr1, sqr2, params.mod);
        MontgomeryMul(point.x, sqr1, sqr2, params);
        MontgomeryMul(point.z, a_2_over_4, cross_1, params);
        ModularAdd(point.z, point.z, sqr2, params.mod);
        MontgomeryMul(point.z, point.z, cross_1, params);
        point.ComputeDiff(params.mod);
    }

    inline void AddPoints(MontgomeryPoint& point, MontgomeryPoint& to_add, MontgomeryPoint& p_min_add){
        if((mpz_cmp_ui(point.x, 0) == 0) && (mpz_cmp_ui(point.z, 0) == 0)){
            return;
        }
        if((mpz_cmp_ui(to_add.x, 0) == 0) && (mpz_cmp_ui(to_add.z, 0) == 0)){
            point.Copy(to_add);
            return;
        }

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

    inline void MultPoints(MontgomeryPoint& point, long multiple){
        if(multiple == 0){
            mpz_set_ui(point.x, 0);
            mpz_set_ui(point.z, 0);
            return;
        }
        if(multiple == 1){
            return;
        }
        if(multiple == 2){
            DoublePoint(point);
            return;
        }
        
        point.Copy(p_copy);
        point.Copy(p_double);
        DoublePoint(p_double);
        
        int bit_count = (int)log2((double)multiple)+1;

        for(int i = bit_count-2; i > 0; i--){
            if(test_bit(multiple, i)){
                AddPoints(p_double, p_copy, point);
                DoublePoint(p_double);
            }
            else{
                AddPoints(p_copy, p_double, point);
                DoublePoint(p_copy);
            }
        }
        if(test_bit(multiple, 0)){
            AddPoints(p_double, p_copy, point);
            p_copy.Copy(point);
        }
        else{
            DoublePoint(p_copy);
            p_copy.Copy(point);
        }
    }
};


struct WeistrassPoint{
    mpz_t x;
    mpz_t y;
    bool is_inf;

    inline WeistrassPoint(){
        mpz_inits(x, y, NULL);
        is_inf = false;
    }

    inline ~WeistrassPoint(){
        mpz_clears(x, y, NULL);
    }

    inline void Copy(WeistrassPoint& point){
        mpz_set(point.x, x);
        mpz_set(point.y, y);
        point.is_inf = is_inf;
    }
};

class WeistrassCurve{
    private:
    WeistrassPoint p_copy;
    mpz_t neg_y, s, to_invert, y2, x_inv;
    
    public:
    mpz_t a, b, n;

    inline WeistrassCurve(mpz_t mod){
        mpz_inits(a, b, n, neg_y, s, to_invert, y2, x_inv, NULL);
        mpz_set(n, mod);
    }

    inline ~WeistrassCurve(){
        mpz_clears(a, b, n, to_invert, y2, x_inv, NULL);
    }

    inline int Initialise(gmp_randstate_t rstate, WeistrassPoint& point, mpz_t n, mpz_t output){
        //Find first point on the curve
        //Randomly select a, x and y
        mpz_urandomm(a, rstate, n);
        mpz_urandomm(point.x, rstate, n);
        mpz_urandomm(point.y, rstate, n);

        #define x_cube neg_y
        #define a_x s
        #define pow_disc y2
        #define b_sqr to_invert
        #define g x_inv
        
        //Compute b value for the curve based on the random coordinates
        mpz_mul(b, point.y, point.y);
        mpz_mul(x_cube, point.x, point.x);
        mpz_mul(x_cube, x_cube, point.x);
        mpz_mul(a_x, point.x, a);
        mpz_sub(b, b, x_cube);
        mpz_sub(b, b, a_x);
        mpz_mod(b, b, n);

        //Check discriminant 
        mpz_pow_ui(pow_disc, a, 3);
        mpz_mul(b_sqr, b, b);
        mpz_mul_ui(b_sqr, b_sqr, 27);
        mpz_mul_ui(pow_disc, pow_disc, 4);
        mpz_add(pow_disc, pow_disc, b_sqr);
        mpz_gcd(g, pow_disc, n);

        if(mpz_cmp(g, n) == 0){
            return 0;
        }
        if(mpz_cmp_ui(g, 1) > 0){
            mpz_set(output, g);
            return 1;
        }

        return 2;
        #undef x_cube
        #undef a_x
        #undef pow_disc
        #undef b_sqr
        #undef g
    }

    inline void AddPoints(WeistrassPoint& point, WeistrassPoint& to_add, mpz_t factor, volatile bool& is_factored){
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

    inline void MultPointsRecursive(WeistrassPoint& point_copy, WeistrassPoint& point, long multiple, mpz_t factor, volatile bool& is_factored){
        if(multiple == 0 || multiple == 1)
            return;

        bool is_odd = (multiple & 1) == 1;
        multiple = multiple >> 1;
        
        MultPointsRecursive(point_copy, point, multiple, factor, is_factored);
        if(is_factored)
            return;

        AddPoints(point, point, factor, is_factored);
        if(is_odd && !is_factored){
            AddPoints(point, point_copy, factor, is_factored);}
    }

    inline void MultPoints(WeistrassPoint& point, long multiple, mpz_t factor, volatile bool& is_factored){
        mpz_set(p_copy.x, point.x);
        mpz_set(p_copy.y, point.y);
        
        MultPointsRecursive(p_copy, point, multiple, factor, is_factored);
    }
};