#include "cuda_kernels.hpp"

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
