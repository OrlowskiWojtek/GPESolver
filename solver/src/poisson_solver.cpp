#include "include/poisson_solver.hpp"
#include "include/output.hpp"
#include "include/params.hpp"
#include <cmath>
#include <fftw3.h>
#include <omp.h>
#include <stdexcept>

int PoissonSolver::FFTW_N_THREADS = 4;

PoissonSolver::PoissonSolver()
    : p(PhysicalParameters::getInstance()) {}

void PoissonSolver::prepare(StdMat3D<std::complex<double>> *cpsi, StdMat3D<double> *_fi3d) {
    prepare_transforms();

    psi  = cpsi;
    fi3d = _fi3d;

    int nx = 2 * p->nx;
    int ny = 2 * p->ny;
    int nz = 2 * p->nz;

    double dx = p->dx;
    double dy = p->dy;
    double dz = p->dz;

    double dkx = 2. * M_PI / (nx * dx);
    double dky = 2. * M_PI / (ny * dy);
    double dkz = 2. * M_PI / (nz * dz);

    for (int i = 0; i < nx; ++i) {
        double kx = (i <= (nx / 2)) ? i * dkx : (i - nx) * dkx;
        for (int j = 0; j < ny; ++j) {
            double ky = (j <= (ny / 2)) ? j * dky : (j - ny) * dky;
            for (int k = 0; k < nz; ++k) {
                double kz = (k <= (nz / 2)) ? k * dkz : (k - nz) * dkz;

                double k2  = kx * kx + ky * ky + kz * kz;
                size_t idx = (i * ny + j) * nz + k;

                if (k2 > 1e-30) {
                    Vdip_k[idx][0] = (kz * kz / k2);
                    Vdip_k[idx][1] = 0.0;
                } else {
                    Vdip_k[idx][0] = 0.0;
                    Vdip_k[idx][1] = 0.0;
                }
            }
        }
    }
}

void PoissonSolver::execute() {
    int nx = 2 * p->nx;
    int ny = 2 * p->ny;
    int nz = 2 * p->nz;
    int N  = nx * ny * nz;

    auto &rpsi = *psi;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                size_t idx = (i * ny + j) * nz + k;

                if (i < nx / 2 && j < ny / 2 && k < nz / 2) {
                    double val    = std::norm(rpsi(i, j, k));
                    rho_r[idx][0] = val * p->n_atoms;
                    rho_r[idx][1] = 0.;
                } else {
                    rho_r[idx][0] = 0.;
                    rho_r[idx][1] = 0.;
                }
            }
        }
    }

    fftw_execute(plan_fwd);

    for (int idx = 0; idx < N; ++idx) {
        double ar     = rho_k[idx][0];
        double ai     = rho_k[idx][1];
        double br     = Vdip_k[idx][0];
        double bi     = Vdip_k[idx][1];
        rho_k[idx][0] = ar * br - ai * bi;
        rho_k[idx][1] = ar * bi + ai * br;
    }

    fftw_execute(plan_bwd);

    double norm_factor = 1.0 / static_cast<double>(N);
    auto &rfi3d        = *fi3d;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                if (i >= nx / 2 || j >= ny / 2 || k >= nz / 2) {
                    continue;
                }

                size_t idx     = (i * ny + j) * nz + k;
                rfi3d(i, j, k) = rho_r[idx][0] * norm_factor;
            }
        }
    }
}

PoissonSolver::~PoissonSolver() {
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);

    fftw_free(rho_r);
    fftw_free(rho_k);
    fftw_free(Vdip_k);
}

void PoissonSolver::prepare_transforms() {
    int nx = 2 * p->nx;
    int ny = 2 * p->ny;
    int nz = 2 * p->nz;
    int N  = nx * ny * nz;

    rho_r  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    rho_k  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    Vdip_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    FFTW_N_THREADS = omp_get_max_threads();

    int res = fftw_init_threads();
    if (res == 0) {
        throw std::runtime_error("FFTW thread initialization failed!");
    }

    fftw_plan_with_nthreads(FFTW_N_THREADS);
    plan_fwd = fftw_plan_dft_3d(nx, ny, nz, rho_r, rho_k, FFTW_FORWARD, FFTW_MEASURE);
    plan_bwd = fftw_plan_dft_3d(nx, ny, nz, rho_k, rho_r, FFTW_BACKWARD, FFTW_MEASURE);

    OutputFormatter::printInfo("Planned FFTW with " + std::to_string(FFTW_N_THREADS) + " threads.");
}
