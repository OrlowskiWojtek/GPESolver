#include "cuda_kernels.hpp"

//void kernel_fill_from_psi(
//    cufftDoubleReal* rho, 
//    const cufftDoubleComplex* psi, 
//    int nx, int ny, int nz, 
//    int full_nx, int full_ny, int full_nz,
//    double n_atoms){
//    int i = blockIdx.x * blockDim.x + threadIdx.x;
//    int j = blockIdx.y * blockDim.y + threadIdx.y;
//    int k = blockIdx.z * blockDim.z + threadIdx.z;
//
//    if (i < nx && j < ny && k < nz) {
//        int src_idx = (i * ny + j) * nz + k;
//        int dst_idx = (i * full_ny + j) * full_nz + k;
//
//        double psi_sq = psi[src_idx].x * psi[src_idx].x + 
//                        psi[src_idx].y * psi[src_idx].y;
//
//        rho[dst_idx] = psi_sq * n_atoms;
//    }
//}

__global__
void kernel_multiply_dipole(
    cufftDoubleComplex* data, 
    const cufftDoubleComplex* kernel, 
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
    cufftDoubleComplex* data, 
    const cufftDoubleComplex* kernel, 
    int N) 
{
    kernel_multiply_dipole<<<(N+255)/256, 256>>>(data, kernel, N);
}

__global__
void kernel_extract_result(
    const cufftDoubleReal* full_grid,
    cufftDoubleReal* output,
    int nx, int ny, int nz,
    double norm_factor)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < nx && j < ny && k < nz) {
        int src_idx = (i * ny + j) * nz + k;
        int dst_idx = (i * ny + j) * nz + k;
        output[dst_idx] = full_grid[src_idx] * norm_factor;
    }
}

void launch_kernel_extract_result(
    const cufftDoubleReal* full_grid,
    cufftDoubleReal* output,
    int nx, int ny, int nz,
    double norm_factor)
{
    kernel_extract_result<<<1, 256>>>(full_grid, output, nx, ny, nz, norm_factor);
}
