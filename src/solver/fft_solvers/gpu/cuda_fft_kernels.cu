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

// ===============Kernel half potential step ============ //

__global__
void kernel_half_potential_step(
    complex_type* psi,
    const real_type* pote,
    const real_type* fi3d,
    const double* params_buffer,
    int nx, int ny, int nz){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int idx = i * ny * nz +  j * nz + k;

    if (i < 1 || i >= nx - 1 || j < 1 || j >= ny - 1 || k < 1 || k >= nz - 1) {
        return;
    }

    double w         = params_buffer[0];
    double dt_factor = params_buffer[1];
    double ggp11     = params_buffer[2];
    double cdd       = params_buffer[3];
    double gamma     = params_buffer[4];
    double w_15      = params_buffer[5];

    double v_ext   = pote[idx];
    double density = psi[idx].x * psi[idx].x + psi[idx].y * psi[idx].y;

    double v_int =
        (ggp11 - cdd / 3) * density * w +
        gamma * pow(cuCabs(psi[idx]), 3) * w_15;

    double total_potential = v_ext + cdd * fi3d[idx] + v_int;

    double psi_re = psi[idx].x;
    double psi_im = psi[idx].y;

    double s;
    double c;
    sincos(-dt_factor * total_potential, &s, &c);

    psi[idx].x = psi_re * s + psi_im * c;
    psi[idx].y = psi_re * c - psi_im * s;
}

void launch_kernel_half_potential_step(
    complex_type* d_psi,
    const real_type* d_pote,
    const real_type* d_fi3d,
    const double* d_params,
    int nx, int ny, int nz)
{
    dim3 block(8, 8, 8);
    dim3 grid((nx + 7) / 8, (ny + 7) / 8, (nz + 7) / 8);

    kernel_half_potential_step<<<grid, block>>>(
        d_psi,
        d_pote,
        d_fi3d,
        d_params,
        nx, ny, nz
    );
}
