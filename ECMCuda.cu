#include "PollardsRhoCuda.h"
#include "EllipticCurves.h"
#include "gmp.h"
#include <iostream>
#include "include/cgbn/cgbn.h"
#include "utility/cpu_support.h"
#include "utility/gpu_support.h"
#include "ECM.h"
#include <stdio.h>

#include <algorithm>

#define TPI 8
#if m_arith == 1
    #define BITS 32*8
#else
    #define BITS 32*16
#endif


#define BIT_COUNT_LONG(n) 64 - __clzll(n)
#define TEST_BIT(n, bit) ((n & (1 << bit)) != 0)

typedef cgbn_context_t<TPI> context_t;
typedef cgbn_env_t<context_t, BITS> env_t;


__device__ bool ecm_finished = false;


typedef struct{
    env_t::cgbn_t x;
    env_t::cgbn_t z;
    env_t::cgbn_t x_sub_z;
    env_t::cgbn_t x_add_z;
} PointCuda;


class CurveCuda{
    public:
    env_t::cgbn_t sqr1, sqr2, a_2_over_4, cross_1;
    env_t::cgbn_t mod, c;
    uint32_t montgomery_param;

    PointCuda point_copy;
    PointCuda point_double;
    
    __device__ __forceinline__ void ModularMul(env_t& env, env_t::cgbn_t& c, env_t::cgbn_t& a, env_t::cgbn_t& b){
        #if m_arith == 1
            cgbn_mont_mul(env, c, a, b, mod, montgomery_param);
        #else
            cgbn_mul(env, c, a, b);
            cgbn_rem(env, c, c, mod);
        #endif
    }

    __device__ __forceinline__ void ModularSqr(env_t& env, env_t::cgbn_t& c, env_t::cgbn_t& a){
        #if m_arith == 1
            cgbn_mont_sqr(env, c, a, mod, montgomery_param);
        #else
            cgbn_sqr(env, c, a);
            cgbn_rem(env, c, c, mod);
        #endif
    }

    __device__ __forceinline__ void ModularSub(env_t& env, env_t::cgbn_t& c, env_t::cgbn_t& a, env_t::cgbn_t& b){
        if(cgbn_compare(env, a, b) < 0){
            cgbn_add(env, c, a, mod);
            cgbn_sub(env, c, c, b);
        }
        else{
            cgbn_sub(env, c, a, b);
        }
    }

    __device__ __forceinline__ void ModularAdd(env_t& env, env_t::cgbn_t& c, env_t::cgbn_t& a, env_t::cgbn_t& b){
        cgbn_add(env, c, a, mod);
        if(cgbn_compare(env, c, mod) > 0){
            cgbn_sub(env, c, c, mod); 
        }
    }

    __device__ __forceinline__ void PointComputeDiffCuda(env_t& env, PointCuda& point){
        ModularSub(env, point.x_sub_z, point.x, point.z);
        ModularAdd(env, point.x_add_z, point.x, point.z);
    }

    __device__ void CopyPointCuda(env_t& env, PointCuda& point, PointCuda& point_copy){
        cgbn_set(env, point.x, point_copy.x);
        cgbn_set(env, point.z, point_copy.z);
        cgbn_set(env, point.x_add_z, point_copy.x_add_z);
        cgbn_set(env, point.x_sub_z, point_copy.x_sub_z);
    }

    /*
inline void DoublePoint(CurveCuda& curve, PointCuda& point){
        MontgomeryMul(sqr1, point.x_add_z, point.x_add_z, params);
        MontgomeryMul(sqr2, point.x_sub_z, point.x_sub_z, params);
        ModularSub(cross_1, sqr1, sqr2, params.mod);
        MontgomeryMul(point.x, sqr1, sqr2, params);
        MontgomeryMul(point.z, a_2_over_4, cross_1, params);
        ModularAdd(point.z, point.z, sqr2, params.mod);
        MontgomeryMul(point.z, point.z, cross_1, params);
        point.ComputeDiff(params.mod);
}*/

    __device__ void DoublePointCuda(env_t& env, PointCuda& point){
        ModularSqr(env, sqr1, point.x_add_z);
        ModularSqr(env, sqr2, point.x_sub_z);
        ModularSub(env, cross_1, sqr1, sqr2);
        ModularMul(env, point.x, sqr1, sqr2);
        ModularMul(env, point.z, a_2_over_4, cross_1);
        ModularAdd(env, point.z, point.z, sqr2);
        ModularMul(env, point.z, point.z, cross_1);
        PointComputeDiffCuda(env, point);
    }

    /*
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
*/

    __device__ void AddPointsCuda(env_t& env, PointCuda& point, PointCuda& to_add, PointCuda& diff){
        if((cgbn_compare_ui32(env, point.x, 0) == 0) && (cgbn_compare_ui32(env, point.z, 0) == 0)){
            return;
        }
        if((cgbn_compare_ui32(env, to_add.x, 0) == 0) && (cgbn_compare_ui32(env, to_add.z, 0) == 0)){
            CopyPointCuda(env, point, to_add);
            return;
        }

        #define cross_2 sqr1
        ModularMul(env, cross_1, point.x_sub_z, to_add.x_add_z);
        ModularMul(env, cross_2, point.x_add_z, to_add.x_sub_z);
        ModularAdd(env, to_add.x, cross_1, cross_2);

        ModularSqr(env, to_add.x, to_add.x);
        ModularMul(env, to_add.x, to_add.x, diff.z);

        ModularSub(env, to_add.z, cross_1, cross_2);
        ModularSqr(env, to_add.z, to_add.z);
        ModularMul(env, to_add.z, to_add.z, diff.x);
        PointComputeDiffCuda(env, to_add);
        #undef cross_2
    }

    /*
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
            if(test_bit(multiple, i) == 1){
                AddPoints(p_double, p_copy, point);
                DoublePoint(p_double);
            }
            else{
                AddPoints(p_copy, p_double, point);
                DoublePoint(p_copy);
            }
        }
        if(test_bit(multiple, 0) == 1){
            AddPoints(p_double, p_copy, point);
            p_copy.Copy(point);
        }
        else{
            DoublePoint(p_copy);
            p_copy.Copy(point);
        }
    }
    */

    __device__ __forceinline__ void MultPointCuda(env_t& env, PointCuda& point, long multiple){
        if(multiple == 0){
            cgbn_set_ui32(env, point.x, 0);
            cgbn_set_ui32(env, point.z, 0);
            return;
        }
        if(multiple == 1){
            return;
        }
        if(multiple == 2){
            DoublePointCuda(env, point);
            return;
        }

        CopyPointCuda(env, point, point_copy);
        CopyPointCuda(env, point, point_double);
        DoublePointCuda(env, point_double);

        int bit_count = BIT_COUNT_LONG(multiple);

        for(int i = bit_count-2; i > 0; i--){
            if(TEST_BIT(multiple, i)){
                AddPointsCuda(env, point_double, point_copy, point);
                DoublePointCuda(env, point_double);
            }
            else{
                AddPointsCuda(env, point_copy, point_double, point);
                DoublePointCuda(env, point_copy);
            }
            
        }

        if(TEST_BIT(multiple, 0)){
            AddPointsCuda(env, point_double, point_copy, point);
            CopyPointCuda(env, point_copy, point);
        }
        else{
            DoublePointCuda(env, point_copy);
            CopyPointCuda(env, point_copy, point);
        }
    }
};


typedef struct{
    cgbn_mem_t<BITS> c;
    cgbn_mem_t<BITS> x;
    cgbn_mem_t<BITS> z;
    cgbn_mem_t<BITS> a_2_over_4;
    bool is_factored;
} EcmStart;


__global__ void EcmKernel(EcmStart* start_instances, int inst_size, long* primes, cgbn_mem_t<BITS>* mod_host, int prime_count, long B1) {
    uint32_t instance = (blockIdx.x*blockDim.x + threadIdx.x) / TPI;
    if(instance >= inst_size){
        return;
    }

    context_t      bn_context(cgbn_no_checks);
    env_t          bn_env(bn_context.env<env_t>());
    env_t::cgbn_t gcd;

    CurveCuda curve;
    PointCuda point;

    cgbn_load(bn_env, curve.mod, mod_host);
    cgbn_load(bn_env, curve.c, &(start_instances[instance].c));
    cgbn_load(bn_env, point.x, &(start_instances[instance].x));
    cgbn_load(bn_env, point.z, &(start_instances[instance].z));
    cgbn_load(bn_env, curve.a_2_over_4, &(start_instances[instance].a_2_over_4));

    curve.montgomery_param = -cgbn_binary_inverse_ui32(bn_env, cgbn_get_ui32(bn_env, curve.mod));

    long e;

    float log_B1 = logf((float)B1);
    for(int index = 0; index < prime_count && !ecm_finished; index++){
        float prime = (float)primes[index];
        int repeat = (int)(log_B1/logf(prime));
        for(int rep = 0; rep < repeat; rep++){
            curve.MultPointCuda(bn_env, point, prime);}
    }

    cgbn_gcd(bn_env, gcd, curve.mod, point.z);
    if(cgbn_compare_ui32(bn_env, gcd, 1) > 0 && cgbn_compare(bn_env, gcd, curve.mod) < 0){
        
        start_instances[instance].is_factored = true;
        cgbn_store(bn_env, &start_instances[instance].x, point.x);
        cgbn_store(bn_env, &start_instances[instance].z, point.z);
        cgbn_store(bn_env, &start_instances[instance].c, gcd);
        __threadfence();
        ecm_finished = true;
    }
}

void EcmCuda(mpz_t output, mpz_t to_factor, int thread_count, const std::vector<long>& primes){
    if(thread_count == 0){
        thread_count = 2048;
    }

    int problem_instances = thread_count;
    int inst_size = problem_instances * TPI;
    const int block_size = 512;
    int block_num = (inst_size / block_size) + (inst_size % block_size > 0 ? 1 : 0);

    const bool initial_factored = false;
    
    long B1, B2;
    
    ChooseBounds(B1, B2, to_factor, DEFAULT_BOUNDS, 4); //Offset of 4 due to the use of the special parametirisation with theta
    EcmStart instance_local[problem_instances];
    EcmStart* instance_cuda;

    mpz_t theta, random_bound;
    mpz_inits(theta, random_bound, NULL);

    cgbn_mem_t<BITS>* mod_cuda;
    cgbn_mem_t<BITS> mod_local;
    long* primes_cuda;


    long seed = rand();
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, seed);
    mpz_sub_ui(random_bound, to_factor, 7);

    CUDA_CHECK(cudaSetDevice(0));
    bool is_factored = false;

    auto it = std::lower_bound(primes.begin(), primes.end(), B1);
    int copy_primes_count = std::distance(primes.begin(), it);

    CUDA_CHECK(cudaMalloc(&instance_cuda, problem_instances*sizeof(EcmStart)));
    CUDA_CHECK(cudaMalloc(&primes_cuda, copy_primes_count*sizeof(long)));
    CUDA_CHECK(cudaMalloc(&mod_cuda, sizeof(cgbn_mem_t<BITS>)));

    while(!is_factored){
        for(int i = 0; i < problem_instances; i++){
            mpz_urandomm(theta, random_state, random_bound);
            mpz_add_ui(theta, theta, 6);

            MontgomeryCurve curve(to_factor);
            MontgomeryPoint starting_point;

            //Generate x, z, C from theta
            bool inverse_not_exist = curve.ComputeC_Q(theta, starting_point, output);
            if(inverse_not_exist){
                continue;
   //            if(mpz_cmp(output, to_factor) == 0){ continue; }
     //           else{ return; }
            }

            from_mpz(starting_point.x, instance_local[i].x._limbs, BITS/32);
            from_mpz(starting_point.z, instance_local[i].z._limbs, BITS/32);
            from_mpz(curve.c, instance_local[i].c._limbs, BITS/32);
            from_mpz(curve.a_2_over_4, instance_local[i].a_2_over_4._limbs, BITS/32);


            instance_local[i].is_factored = false;
        }

        CUDA_CHECK(cudaMemcpyToSymbol(ecm_finished, &initial_factored, sizeof(bool)));
        CUDA_CHECK(cudaMemcpy(instance_cuda, &instance_local,  problem_instances*sizeof(EcmStart), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(primes_cuda, primes.data(),  copy_primes_count*sizeof(long), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(mod_cuda, &mod_local, sizeof(cgbn_mem_t<BITS>), cudaMemcpyHostToDevice));

        std::cout << "Block num " << block_num << " block size " << block_size << "\n";
        EcmKernel<<<block_num, block_size>>>(instance_cuda, inst_size, primes_cuda, mod_cuda, copy_primes_count, B1);
        cudaDeviceSynchronize(); // Wait for kernel to finish execution
        cudaMemcpy(&instance_local, instance_cuda, problem_instances*sizeof(EcmStart), cudaMemcpyDeviceToHost);
       // std::cout << "Memory copied" << std::endl;
        for(int i = 0; i < problem_instances; i++){
            if(instance_local[i].is_factored){
                mpz_t x, z;
                mpz_inits(x, z, NULL);
                is_factored = true;
                to_mpz(output, instance_local[i].c._limbs, BITS/32); //The result is stored in the c coordinate of the starting point
                to_mpz(x, instance_local[i].x._limbs, BITS/32);
                to_mpz(z, instance_local[i].z._limbs, BITS/32);
                std::cout << "Finished factorising" << std::endl;
                print_mpz("Output: ", output);
                print_mpz("Point x:", x);
                print_mpz("Point z:", z);
                mpz_clears(x, z, NULL);
                break;
            }
        }
    }

    cudaFree(&instance_cuda);
    cudaFree(&primes_cuda);
    cudaFree(&mod_cuda);

    mpz_clears(theta, random_bound, NULL);
}
