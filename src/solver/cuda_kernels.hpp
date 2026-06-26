#ifndef IMAG_TIME_KERNEL_CUH
#define IMAG_TIME_KERNEL_CUH

#include <cuComplex.h>
#include "context/context.hpp"

__global__ void imag_time_iteration_kernel(
    const cuDoubleComplex *__restrict__ cpsi, // funkcja falowa (tylko do odczytu)
    const double *__restrict__ pote,          // potencjał zewnętrzny (tylko do odczytu)
    const double *__restrict__ fi3d,          // potencjał dipolowy (tylko do odczytu)
    cuDoubleComplex *__restrict__ cpsii,      // wynik
    int nx,
    int ny,
    int nz,
    double dx,
    double dy,
    double dz,
    double m,
    double imag_dt,
    double cdd,
    double ggp11,
    double gamma,
    double w,   // n_atoms
    double w_15 // w_15
);

void launch_kernel_imag_time_iteration(const cuDoubleComplex *d_cpsi,
                                       const double *d_pote,
                                       const double *d_fi3d,
                                       cuDoubleComplex *d_cpsii,
                                       int nx,
                                       int ny,
                                       int nz,
                                       double dx,
                                       double dy,
                                       double dz,
                                       double m,
                                       double imag_dt,
                                       double cdd,
                                       double ggp11,
                                       double gamma,
                                       double w,
                                       double w_15);
// ================================ //
__global__ void kernel_normalize(
    cuDoubleComplex* __restrict__ data,
    int N,
    double norm_factor
);

void launch_kernel_normalize(
    cuDoubleComplex* data,
    int N,
    double norm_factor
);

__global__ void kernel_calc_norm_sq(
    const cuDoubleComplex* __restrict__ data,
    double* __restrict__ result,
    int N
);

double launch_kernel_calc_norm(
    const cuDoubleComplex* data,
    int N
);

// ================================== //

__global__
void kernel_calc_energies(
    const cuDoubleComplex* __restrict__ psi,
    const double* __restrict__ pote,
    const double* __restrict__ fi3d,
    double* __restrict__ e_kin,
    double* __restrict__ e_pot,
    double* __restrict__ e_int,
    double* __restrict__ e_ext,
    double* __restrict__ e_bmf,
    int nx, int ny, int nz,
    double dx, double dy, double dz,
    double m, double ggp11, double cdd, double gamma,
    double n_atoms, double w_15
);

void launch_kernel_calc_energies(
    const cuDoubleComplex* psi,
    const double* pote,
    const double* fi3d,
    energies_t& ene,
    int nx, int ny, int nz,
    double dx, double dy, double dz,
    double m, double ggp11, double cdd, double gamma,
    double n_atoms, double w_15
);

#endif
