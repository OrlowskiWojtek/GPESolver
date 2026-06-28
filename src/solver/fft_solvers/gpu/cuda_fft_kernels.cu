#include "solver/fft_solvers/gpu/cuda_fft_kernels.hpp"
#include "cuComplex.h"

// ========Kernel multiply dipole======= //

__global__
void kernel_multiply_dipole(
    complex_type* data, 
    const complex_type* kernel, 
    int N) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        double ar = data[idx].x;
        double ai = data[idx].y;
        double br = kernel[idx].x;
        double bi = kernel[idx].y;
        
        data[idx].x = ar * br - ai * bi;
        data[idx].y = ar * bi + ai * br;
    }
}

void launch_kernel_multiply_dipole(
    complex_type* data, 
    const complex_type* kernel, 
    int N) 
{
    kernel_multiply_dipole<<<(N+255)/256, 256>>>(data, kernel, N);
}

// ========Kernel fill psi======= //

__global__
void kernel_fill_from_psi(
    real_type* rho, 
    const complex_type* psi, 
    int nx, int ny, int nz, 
    int full_nx, int full_ny, int full_nz,
    double n_atoms){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int rho_idx = (i * full_ny + j) * full_nz + k;
    int psi_idx = (i * ny + j) * nz + k;
    if (i >= full_nx || j >= full_ny || k >= full_nz) {
        return;
    }

    if (i < nx && j < ny && k < nz) {
        double psi_sq = psi[psi_idx].x * psi[psi_idx].x + 
                        psi[psi_idx].y * psi[psi_idx].y;

        rho[rho_idx] = psi_sq * n_atoms;
    } else {
        rho[rho_idx] = 0;
    }
}

void launch_kernel_fill_from_psi(
    real_type* rho, 
    const complex_type* psi, 
    int nx, int ny, int nz, 
    int full_nx, int full_ny, int full_nz,
    double n_atoms){

    int blockSize = 8;
    dim3 blocks((full_nx+blockSize-1)/blockSize, 
                (full_ny+blockSize-1)/blockSize, 
                (full_nz+blockSize-1)/blockSize);
    dim3 threads(blockSize, blockSize, blockSize);
    kernel_fill_from_psi<<<blocks, threads>>>(rho, psi, nx, ny, nz, 
                                              full_nx, full_ny, full_nz,
                                              n_atoms);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
}

// ========Kernel kinetic======= //

__global__
void kernel_kinetic(
    complex_type* rho, 
    const real_type* kernel, 
    int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < N) {   
        double s;
        double c;
        
        double psi_re = rho[idx].x;
        double psi_im = rho[idx].y;
        sincos(kernel[idx], &s, &c);

        rho[idx].x = psi_re * c - psi_im * s;
        rho[idx].y = psi_re * s + psi_im * c;
    }
}

void launch_kernel_kinetic(
    complex_type* rho, 
    const real_type* kernel, 
    int N)
{
    kernel_kinetic<<<(N+255)/256, 256>>>(rho, kernel, N);
}

// ===============Kernel copy to fi3d data ============ //

__global__ void kernel_copy_to_fi3d_gpu(
    const real_type* __restrict__ d_rho_r,
    double* __restrict__ fi3d_gpu,
    int nx, int ny, int nz,
    int full_nx, int full_ny, int full_nz,
    double norm_factor
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;


    if (i < nx && j < ny && k < nz) {
        size_t idx = (static_cast<size_t>(i) * ny + j) * nz + k;
        size_t idx_full = (static_cast<size_t>(i) * full_ny + j) * full_nz + k;
        fi3d_gpu[idx] = d_rho_r[idx_full] * norm_factor;
    }
}

void launch_kernel_copy_to_fi3d_gpu(
    real_type* d_rho_r,
    double* fi3d_gpu,
    int nx, int ny, int nz,
    int full_nx, int full_ny, int full_nz,
    double norm_factor
) {
    dim3 block(8, 8, 8);
    dim3 grid((nx + block.x - 1) / block.x,
              (ny + block.y - 1) / block.y,
              (nz + block.z - 1) / block.z);

    kernel_copy_to_fi3d_gpu<<<grid, block>>>(
        d_rho_r, fi3d_gpu, nx, ny, nz,
        full_nx, full_ny, full_nz, norm_factor
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
}

__global__ 
void kernel_copy_with_norm(
    const complex_type* __restrict__ src,
    cuDoubleComplex* __restrict__ dst,
    double norm_factor,
    int N
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < N) {
        dst[idx].x = src[idx].x * norm_factor;
        dst[idx].y = src[idx].y * norm_factor;
    }
}

void launch_kernel_copy_with_norm(
    const complex_type* src,
    cuDoubleComplex* dst,
    int N
) {
    int block = 256;
    int grid = (N + block - 1) / block;
    
    double norm_factor = 1.0 / static_cast<double>(N);
    
    kernel_copy_with_norm<<<grid, block>>>(src, dst, norm_factor, N);
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
}
