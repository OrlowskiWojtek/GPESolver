#ifndef BINARY_MIXTURES_KERNELS_CUH
#define BINARY_MIXTURES_KERNELS_CUH

#include <cuComplex.h>

#define CUDA_KERNEL_DECLARE(kernel_name, ...) \
    __global__ void kernel_name(__VA_ARGS__); \
    void launch_##kernel_name(__VA_ARGS__)

bool initialize_gauss();

void launch_kernel_calc_lhy(
    const cuDoubleComplex* data,
    double* __restrict__ d_norm,
    const double dxdydz,
    int N
);

void launch_kernel_imag_time_iteration(
    cuDoubleComplex *__restrict__ cpsi_a,
    cuDoubleComplex *__restrict__ cpsi_b,
    const double *__restrict__ pote,
    const double *__restrict__ fi3d_a,
    const double *__restrict__ fi3d_b,
    const double *__restrict__ flhy_a,
    const double *__restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double n_atoms_a_15,
    const double n_atoms_b_15,
    const double cdd_11,
    const double cdd_12,
    const double cdd_22,
    const double ggp_11,
    const double ggp_12,
    const double ggp_22,
    const int nx,
    const int ny,
    const int nz,
    const double dx,
    const double dy,
    const double dz,
    const double dt
);

__global__ 
void kernel_imag_time_iteration(
    cuDoubleComplex *__restrict__ cpsi_a,
    cuDoubleComplex *__restrict__ cpsi_b,
    const double *__restrict__ pote,
    const double *__restrict__ fi3d_a,
    const double *__restrict__ fi3d_b,
    const double *__restrict__ flhy_a,
    const double *__restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double n_atoms_a_15,
    const double n_atoms_b_15,
    const double cdd_11,
    const double cdd_12,
    const double cdd_22,
    const double ggp_11,
    const double ggp_12,
    const double ggp_22,
    const int nx,
    const int ny,
    const int nz,
    const double dx,
    const double dy,
    const double dz,
    const double dt
);

void launch_kernel_calc_lhy(
    const cuDoubleComplex* __restrict__ psi_a,
    const cuDoubleComplex* __restrict__ psi_b,
    double* __restrict__ flhy_a,
    double* __restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double g11,
    const double g12,
    const double g22,
    const double cdd11,
    const double cdd12,
    const double cdd22,
    const int nx,
    const int ny,
    const int nz
);

__global__
void kernel_calc_lhy(
    const cuDoubleComplex* __restrict__ psi_a,
    const cuDoubleComplex* __restrict__ psi_b,
    double* __restrict__ flhy_a,
    double* __restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double g11,
    const double g12,
    const double g22,
    const double cdd11,
    const double cdd12,
    const double cdd22,
    const int nx,
    const int ny,
    const int nz
);

#endif
