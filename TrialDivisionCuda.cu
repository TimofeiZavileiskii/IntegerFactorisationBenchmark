#include "TrialDivisionCuda.h"
#include <gmp.h>
#include <iostream>
#include "include/cgbn/cgbn.h"
#include "utility/cpu_support.h"
#include "utility/gpu_support.h"

#define TPI 4
#define BITS 128

typedef cgbn_context_t<TPI> context_t;
typedef cgbn_env_t<context_t, BITS> env_t;


__device__ uint32_t range_size = 1;
__device__ bool factored = false;


typedef struct {
  cgbn_mem_t<BITS> to_factor;
  cgbn_mem_t<BITS> output;
} instance_t;


__global__ void try_division(instance_t* problem) {
    int32_t instance = (blockIdx.x*blockDim.x + threadIdx.x) / TPI;

    context_t      bn_context(cgbn_no_checks);
    env_t          bn_env(bn_context.env<env_t>());
    
    env_t::cgbn_t to_factor, range_start, try_mod, mod_out;

    cgbn_load(bn_env, to_factor, &(problem->to_factor));
    cgbn_set_ui32(bn_env, range_start, range_size);
    cgbn_mul_ui32(bn_env, range_start, range_start, instance);
    cgbn_add_ui32(bn_env, range_start, range_start, 2);

    if(instance == 0 && threadIdx.x % TPI == 0){
        printf(">> %i\n", range_size);
    }

    for(uint32_t i = 0; i < range_size && !factored; i++){
        cgbn_add_ui32(bn_env, try_mod, range_start, i);
        cgbn_rem(bn_env, mod_out, to_factor, try_mod);
        if(cgbn_compare_ui32(bn_env, mod_out, 0) == 0){
            env_t::cgbn_t other_factor;
            cgbn_div(bn_env, other_factor, to_factor, try_mod);
            
            if(cgbn_compare(bn_env, other_factor, try_mod) != -1){
                cgbn_store(bn_env, &(problem->output), try_mod);
            }
            else{
                cgbn_store(bn_env, &(problem->output), other_factor);
            }
            factored = true;
        }
    }
}


void TrialDivisionCuda(mpz_t& output, mpz_t& to_factor)
{
    const int num_kernels = 268435456;
    const int inst_size = num_kernels * TPI;
    const int block_size = 512;
    const int block_num = inst_size / block_size;
    
    const bool initial_factored = false;

    mpz_t bound, range_size_local;
    mpz_inits(bound, range_size_local, NULL);
    double bound_d = sqrt(mpz_get_d(to_factor));
    mpz_set_d(bound, bound_d);
    mpz_div_ui(range_size_local, bound, num_kernels);
    mpz_add_ui(range_size_local, range_size_local, 1);;
    
    uint32_t range_size_local_ui = mpz_get_ui(range_size_local);
    instance_t probelm_instance_local;
    instance_t* probelm_instance_cuda;

    CUDA_CHECK(cudaSetDevice(0));

    CUDA_CHECK(cudaMalloc(&probelm_instance_cuda, sizeof(instance_t)));

    from_mpz(to_factor, probelm_instance_local.to_factor._limbs, BITS/32);
    CUDA_CHECK(cudaMemcpyToSymbol(range_size, &range_size_local_ui, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(factored, &initial_factored, sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(probelm_instance_cuda, &probelm_instance_local, sizeof(instance_t), cudaMemcpyHostToDevice));

    try_division<<<block_num, block_size>>>(probelm_instance_cuda);

    cudaMemcpy(&probelm_instance_local, probelm_instance_cuda, sizeof(instance_t), cudaMemcpyDeviceToHost);
    cudaFree(&probelm_instance_cuda);
    to_mpz(output, probelm_instance_local.output._limbs, BITS/32);

    mpz_clears(bound, range_size_local, NULL);
}