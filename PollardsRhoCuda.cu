#include "PollardsRhoCuda.h"
#include "gmp.h"
#include <iostream>
#include "include/cgbn/cgbn.h"
#include "utility/cpu_support.h"
#include "utility/gpu_support.h"

#define TPI 8
#define BITS 32*8

typedef cgbn_context_t<TPI> context_t;
typedef cgbn_env_t<context_t, BITS> env_t;


typedef struct{
    cgbn_mem_t<BITS> to_factor;
    cgbn_mem_t<BITS> start_point;
    cgbn_mem_t<BITS> a;
    bool is_factored;
} PollardStart;


__device__ bool factored = false;


__global__ void TryPollardRho(PollardStart* starting) {
    uint32_t instance = (blockIdx.x*blockDim.x + threadIdx.x) / TPI;

    context_t      bn_context(cgbn_no_checks);
    env_t          bn_env(bn_context.env<env_t>());
    
    env_t::cgbn_t to_factor, diff, point, point2, a, product;

    cgbn_load(bn_env, to_factor, &(starting[instance].to_factor));
    cgbn_load(bn_env, point, &(starting[instance].start_point));
    cgbn_set(bn_env, point2, point);
    cgbn_load(bn_env, a, &(starting->a));
    cgbn_set_ui32(bn_env, product, 1);

    uint32_t mont_c = -cgbn_binary_inverse_ui32(bn_env, cgbn_get_ui32(bn_env, to_factor));
    int count = 0;

    while(!factored){
        cgbn_mont_sqr(bn_env, point, point, to_factor, mont_c);
        cgbn_add(bn_env, point, point, a);
        cgbn_mont_sqr(bn_env, point2, point2, to_factor, mont_c);
        cgbn_add(bn_env, point2, point2, a);
        cgbn_mont_sqr(bn_env, point2, point2, to_factor, mont_c);
        cgbn_add(bn_env, point2, point2, a);


        if(cgbn_compare(bn_env, point, point2) < 0){
            cgbn_sub(bn_env, diff, point2, point);
        }
        else{
            cgbn_sub(bn_env, diff, point, point2);
        }
        
        if(cgbn_compare_ui32(bn_env, diff, 0) == 0){
            return;
        }

        cgbn_mul(bn_env, product, product, diff);
        cgbn_rem(bn_env, product, product, to_factor);

        count++;
        if(count == 100){
            count = 0;
            cgbn_gcd(bn_env, product, product, to_factor);//here product is the factor
            
            if(cgbn_compare_ui32(bn_env, product, 1) > 0 && cgbn_compare(bn_env, product, to_factor) < 0){
                factored = true;
                cgbn_store(bn_env, &(starting[instance].start_point), product);
                starting[instance].is_factored = true;
            }
            cgbn_set_ui32(bn_env, product, 1);
        }
    }
}


void PollardsRhoCuda(mpz_t output, mpz_t to_factor, int thread_count){
    if(thread_count == 0){
        thread_count = 2048;
    }

    int problem_instances = thread_count;
    int inst_size = problem_instances * TPI;
    const int block_size = 512;
    int block_num = inst_size / block_size;

    const bool initial_factored = false;
    
    PollardStart instance_local[problem_instances];
    PollardStart* instance_cuda;

    mpz_t temp_random;
    mpz_init(temp_random);

    long seed = rand();
    gmp_randstate_t random_state;
    gmp_randinit_default(random_state);
    gmp_randseed_ui(random_state, seed);

    CUDA_CHECK(cudaSetDevice(0));

    for(int i = 0; i < problem_instances; i++){
        PollardStart start;
        mpz_urandomm(temp_random, random_state, to_factor);
        from_mpz(temp_random, start.a._limbs, BITS/32);
        mpz_urandomm(temp_random, random_state, to_factor);
        from_mpz(temp_random, start.start_point._limbs, BITS/32);
        from_mpz(to_factor, start.to_factor._limbs, BITS/32);
        start.is_factored = false;
        instance_local[i] = start;
    }

    CUDA_CHECK(cudaMalloc(&instance_cuda, problem_instances * sizeof(PollardStart)));
    CUDA_CHECK(cudaMemcpyToSymbol(factored, &initial_factored, sizeof(bool)));
    CUDA_CHECK(cudaMemcpy(instance_cuda, &instance_local,  problem_instances * sizeof(PollardStart), cudaMemcpyHostToDevice));

    TryPollardRho<<<block_num, block_size>>>(instance_cuda);

    cudaMemcpy(&instance_local, instance_cuda, problem_instances*sizeof(PollardStart), cudaMemcpyDeviceToHost);
    cudaFree(&instance_cuda);

    for(int i = 0; i < problem_instances; i++){
        if(instance_local[i].is_factored){
            to_mpz(output, instance_local[i].start_point._limbs, BITS/32);
            break;
        }
    }

    mpz_clear(temp_random);
}
