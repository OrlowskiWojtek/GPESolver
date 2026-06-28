#ifndef CUDA_KERNELS_HPP
#define CUDA_KERNELS_HPP

#include "cufft.h"
#include "solver/fft_solvers/fft_context.hpp"

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

// =============================== //
// TODO new -> to refactor
__global__ void kernel_copy_to_fi3d_gpu(
    const real_type* __restrict__ d_rho_r,
    double* __restrict__ fi3d_gpu,
    int nx, int ny, int nz,
    int full_nx, int full_ny, int full_nz,
    double norm_factor
);

void launch_kernel_copy_to_fi3d_gpu(
    real_type* d_rho_r,
    double* fi3d_gpu,
    int nx, int ny, int nz,
    int full_nx, int full_ny, int full_nz,
    double norm_factor
);

// ============================= //
__global__ 
void kernel_copy_with_norm(
    const complex_type* __restrict__ src,
    cuDoubleComplex* __restrict__ dst,
    double norm_factor,
    int N
);

void launch_kernel_copy_with_norm(
    const complex_type* src,
    cuDoubleComplex* dst,
    int N
);

#ifdef __cplusplus
}
#endif

#endif
