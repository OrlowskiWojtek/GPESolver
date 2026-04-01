#include "solver/fft_poisson_solver.hpp"
#include "output.hpp"
#include "parameters/parameters.hpp"
#include <cmath>

#ifdef USE_CUDA
#include "cuda_kernels.hpp"
#endif

PoissonSolver::PoissonSolver() {
}

void PoissonSolver::prepare_containers() {
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

#ifdef USE_CUDA
    auto err = cudaMalloc(&d_Vdip_k, sizeof(complex_type) * N);
    if (err != cudaSuccess) {
        OutputFormatter::printError("Can't aloc Vdip memory");
        OutputFormatter::printError(cudaGetErrorString(err));
    }

    cudaMemcpy(d_Vdip_k, h_Vdip_k, N * sizeof(complex_type), cudaMemcpyHostToDevice);
#endif
}

void PoissonSolver::execute() {
    int nx    = 2 * p->nx;
    int ny    = 2 * p->ny;
    int nz    = 2 * p->nz;
    int N     = nx * ny * nz;
    int N_out = nx * ny * (nz / 2 + 1);

#ifdef USE_CUDA
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
#else
    auto &rpsi = *psi;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                size_t idx = (i * ny + j) * nz + k;
                if (i < nx / 2 && j < ny / 2 && k < nz / 2) {
                    h_rho_r[idx] = std::norm(rpsi(i, j, k)) * p->n_atoms;
                } else {
                    h_rho_r[idx] = 0;
                }
            }
        }
    }

    fftw_execute(plan_fwd);

    for (int idx = 0; idx < N_out; ++idx) {
        double ar          = real(h_rho_k[idx]);
        double ai          = imag(h_rho_k[idx]);
        double br          = real(h_Vdip_k[idx]);
        double bi          = imag(h_Vdip_k[idx]);
        real(h_rho_k[idx]) = ar * br - ai * bi;
        imag(h_rho_k[idx]) = ar * bi + ai * br;
    }

    fftw_execute(plan_bwd);
#endif

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

PoissonSolver::~PoissonSolver() {
#ifdef USE_CUDA
    cudaFree(d_rho_r);
    cudaFree(d_rho_k);
    cudaFree(d_Vdip_k);
    cudaFree(d_psi);

    cudaFreeHost(h_rho_r);
#else
    fftw_free(h_rho_r);
    fftw_free(h_rho_k);
#endif

    delete[] h_Vdip_k;
}

void PoissonSolver::prepare_transforms() {
    int nx    = 2 * p->nx;
    int ny    = 2 * p->ny;
    int nz    = 2 * p->nz;
    int N     = nx * ny * nz;
    int N_out = nx * ny * (nz / 2 + 1);

#ifdef USE_CUDA
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
#else
    h_rho_r = (real_type *)fftw_malloc(sizeof(real_type) * N);
    h_rho_k = (complex_type *)fftw_malloc(sizeof(complex_type) * N_out);

    fftw_plan_with_nthreads(FFTW_N_THREADS);
    plan_fwd = fftw_plan_dft_r2c_3d(nx, ny, nz, h_rho_r, h_rho_k, FFTW_MEASURE);
    fftw_plan_with_nthreads(FFTW_N_THREADS);
    plan_bwd = fftw_plan_dft_c2r_3d(nx, ny, nz, h_rho_k, h_rho_r, FFTW_MEASURE);

    OutputFormatter::printInfo("Planned FFTW with " + std::to_string(FFTW_N_THREADS) +
                               " threads (Poisson).");
#endif
}
