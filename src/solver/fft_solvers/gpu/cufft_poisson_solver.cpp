#include "solver/fft_solvers/gpu/cufft_poisson_solver.hpp"
#include "output.hpp"
#include "parameters/parameters.hpp"
#include <cmath>

#include "solver/fft_solvers/gpu/cuda_fft_kernels.hpp"

CUFFTPoissonSolver::CUFFTPoissonSolver() {
}

void CUFFTPoissonSolver::prepare_containers() {
    int nx = 2 * p->nx;
    int ny = 2 * p->ny;
    int nz = 2 * p->nz;

    double dx = p->dx;
    double dy = p->dy;
    double dz = p->dz;

    double dkx = 2. * M_PI / (nx * dx);
    double dky = 2. * M_PI / (ny * dy);
    double dkz = 2. * M_PI / (nz * dz);

    int N = nx * ny * (nz / 2 + 1);

    h_Vdip_k = new complex_type[N];

    for (int i = 0; i < nx; ++i) {
        double kx = (i <= (nx / 2)) ? i * dkx : (i - nx) * dkx;
        for (int j = 0; j < ny; ++j) {
            double ky = (j <= (ny / 2)) ? j * dky : (j - ny) * dky;
            for (int k = 0; k < nz / 2 + 1; ++k) {
                size_t idx = (i * ny + j) * (nz / 2 + 1) + k;
                double kz  = (k <= (nz / 2)) ? k * dkz : (k - nz) * dkz;

                double k2 = kx * kx + ky * ky + kz * kz;

                if (k2 > 1e-30) {
                    real(h_Vdip_k[idx]) = (kz * kz / k2);
                    imag(h_Vdip_k[idx]) = 0.0;
                } else {
                    real(h_Vdip_k[idx]) = 0.0;
                    imag(h_Vdip_k[idx]) = 0.0;
                }
            }
        }
    }

    auto err = cudaMalloc(&d_Vdip_k, sizeof(complex_type) * N);
    if (err != cudaSuccess) {
        OutputFormatter::printError("Can't aloc Vdip memory");
        OutputFormatter::printError(cudaGetErrorString(err));
    }

    cudaMemcpy(d_Vdip_k, h_Vdip_k, N * sizeof(complex_type), cudaMemcpyHostToDevice);
}

void CUFFTPoissonSolver::execute() {
    int nx    = 2 * p->nx;
    int ny    = 2 * p->ny;
    int nz    = 2 * p->nz;
    int N     = nx * ny * nz;
    int N_out = nx * ny * (nz / 2 + 1);

    // Copy psi_data onto device
    cudaMemcpy(d_psi, psi->get_data(), psi->size() * sizeof(complex_type), cudaMemcpyHostToDevice);

    // get BEC density from psi
    launch_kernel_fill_from_psi(d_rho_r, d_psi, p->nx, p->ny, p->nz, nx, ny, nz, p->n_atoms);

    //  Forward FFT
    cufftExecD2Z(plan_fwd, d_rho_r, d_rho_k);

    // multiply by dipole kernel
    launch_kernel_multiply_dipole(d_rho_k, d_Vdip_k, N_out);

    // Backward FFT
    cufftExecZ2D(plan_bwd, d_rho_k, d_rho_r);
    cudaMemcpy(h_rho_r, d_rho_r, N * sizeof(real_type), cudaMemcpyDeviceToHost);

    double norm_factor = 1.0 / static_cast<double>(N);
    auto &rfi3d        = *fi3d;

    for (int i = 0; i < nx / 2; i++) {
        for (int j = 0; j < ny / 2; j++) {
            for (int k = 0; k < nz / 2; k++) {
                size_t idx     = (i * ny + j) * nz + k;
                rfi3d(i, j, k) = h_rho_r[idx] * norm_factor;
            }
        }
    }
}

CUFFTPoissonSolver::~CUFFTPoissonSolver() {
    cudaFree(d_rho_r);
    cudaFree(d_rho_k);
    cudaFree(d_Vdip_k);
    cudaFree(d_psi);

    cudaFreeHost(h_rho_r);

    delete[] h_Vdip_k;
}

void CUFFTPoissonSolver::prepare_transforms() {
    int nx    = 2 * p->nx;
    int ny    = 2 * p->ny;
    int nz    = 2 * p->nz;
    int N     = nx * ny * nz;
    int N_out = nx * ny * (nz / 2 + 1);

    auto err = cudaMalloc(&d_rho_r, sizeof(real_type) * N);
    if (err != cudaSuccess) {
        OutputFormatter::printError("Can't aloc d_rho_r memory");
        OutputFormatter::printError(cudaGetErrorString(err));
    }

    err = cudaMalloc(&d_rho_k, sizeof(complex_type) * N_out);
    if (err != cudaSuccess) {
        OutputFormatter::printError("Can't aloc d_rho_k memory");
        OutputFormatter::printError(cudaGetErrorString(err));
    }

    err = cudaMalloc(&d_psi, sizeof(complex_type) * psi->size());
    if (err != cudaSuccess) {
        OutputFormatter::printError("Can't aloc d_psi memory");
        OutputFormatter::printError(cudaGetErrorString(err));
    }

    cudaMallocHost(&h_rho_r, sizeof(real_type) * N);

    cufftPlan3d(&plan_fwd, nx, ny, nz, CUFFT_D2Z);
    cufftPlan3d(&plan_bwd, nx, ny, nz, CUFFT_Z2D);

    OutputFormatter::printInfo("Planned FFT with cuFFT");
}
