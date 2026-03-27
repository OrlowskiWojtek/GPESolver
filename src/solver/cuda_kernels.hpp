#ifndef CUDA_KERNELS_HPP
#define CUDA_KERNELS_HPP

#include "cufft.h"
#include "solver/fft_context.hpp"

#ifdef __cplusplus
extern "C" {
#endif

// ============================== //
__global__
void kernel_fill_from_psi(
    cufftDoubleReal* rho, 
    const cufftDoubleComplex* psi, 
    int nx, int ny, int nz, 
    int full_nx, int full_ny, int full_nz, 
    double n_atoms);

void launch_kernel_fill_from_psi(
    cufftDoubleReal* rho, 
    const cufftDoubleComplex* psi, 
    int nx, int ny, int nz, 
    int full_nx, int full_ny, int full_nz, 
    double n_atoms);

// ============================== //

__global__
void kernel_multiply_dipole(
    complex_type* data, 
    const complex_type* kernel, 
    int N);

void launch_kernel_multiply_dipole(
    complex_type* data, 
    const complex_type* kernel, 
    int N);

// ============================== //

__global__
void kernel_kinetic(
    complex_type* rho, 
    const real_type* kernel, 
    int N);

void launch_kernel_kinetic(
    complex_type* rho, 
    const real_type* kernel, 
    int N);

// ============================== //

__global__
void kernel_half_potential_step(
    complex_type* psi,
    real_type* pote
);

void launch_kernel_half_potential_step(
    complex_type* psi,
    real_type* pote
);

#ifdef __cplusplus
}
#endif

#endif

