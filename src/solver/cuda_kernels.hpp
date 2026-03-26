#ifndef CUDA_KERNELS_HPP
#define CUDA_KERNELS_HPP

#include "cufft.h"

#ifdef __cplusplus
extern "C" {
#endif

//void kernel_fill_from_psi(
//    cufftDoubleReal* rho, 
//    const cufftDoubleComplex* psi, 
//    int nx, int ny, int nz, 
//    int full_nx, int full_ny, int full_nz,
//    double n_atoms);

__global__
void kernel_multiply_dipole(
    cufftDoubleComplex* data, 
    const cufftDoubleComplex* kernel, 
    int N);

__global__
void kernel_extract_result(
    const cufftDoubleReal* full_grid,
    cufftDoubleReal* output,
    int nx, int ny, int nz,
    int full_nx, int full_ny, int full_nz,
    double norm_factor);

void launch_kernel_multiply_dipole(
    cufftDoubleComplex* data, 
    const cufftDoubleComplex* kernel, 
    int N);

void launch_kernel_extract_result(
    const cufftDoubleReal* full_grid,
    cufftDoubleReal* output,
    int nx, int ny, int nz,
    double norm_factor);

#ifdef __cplusplus
}
#endif

#endif

