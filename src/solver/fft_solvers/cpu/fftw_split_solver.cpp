#include "solver/fft_solvers/cpu/fftw_split_solver.hpp"
#include "output.hpp"

FFTWRealTimeSplitSolver::FFTWRealTimeSplitSolver() {
}

void FFTWRealTimeSplitSolver::prepare_containers() {
    int nx = p->nx;
    int ny = p->ny;
    int nz = p->nz;

    double dx = p->dx;
    double dy = p->dy;
    double dz = p->dz;

    int N            = nx * ny * nz;
    double dt        = p->real_time_dt;
    h_kinetic_factor = new double[N];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double kx = (i <= nx / 2) ? (2.0 * M_PI * i / (nx * dx))
                                          : (2.0 * M_PI * (i - nx) / (nx * dx));
                double ky = (j <= ny / 2) ? (2.0 * M_PI * j / (ny * dy))
                                          : (2.0 * M_PI * (j - ny) / (ny * dy));
                double kz = (k <= nz / 2) ? (2.0 * M_PI * k / (nz * dz))
                                          : (2.0 * M_PI * (k - nz) / (nz * dz));

                double k_sq = kx * kx + ky * ky + kz * kz;

                double factor = -dt / (2.0 * p->m) * k_sq;

                h_kinetic_factor[i * ny * nz + j * nz + k] = factor;
            }
        }
    }
}

void FFTWRealTimeSplitSolver::execute() {
    int nx = p->nx;
    int ny = p->ny;
    int nz = p->nz;
    int N  = nx * ny * nz;

    auto &rpsi = *psi;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                size_t idx = (i * ny + j) * nz + k;

                real(h_rho_r[idx]) = rpsi(i, j, k).real();
                imag(h_rho_r[idx]) = rpsi(i, j, k).imag();
            }
        }
    }

    fftw_execute(plan_fwd);

    for (int i = 0; i < N; i++) {
        double factor = h_kinetic_factor[i];

        std::complex<double> evolution_operator = std::exp(std::complex<double>(0.0, factor));
        std::complex<double> psi_k(h_rho_k[i][0], h_rho_k[i][1]);
        std::complex<double> res = psi_k * evolution_operator;
        h_rho_k[i][0]            = res.real();
        h_rho_k[i][1]            = res.imag();
    }

    fftw_execute(plan_bwd);

    double norm_factor = 1.0 / static_cast<double>(N);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                size_t idx = (i * ny + j) * nz + k;
                rpsi(i, j, k) =
                    std::complex<double>(real(h_rho_r[idx]), imag(h_rho_r[idx])) * norm_factor;
            }
        }
    }
}

FFTWRealTimeSplitSolver::~FFTWRealTimeSplitSolver() {
    fftw_free(h_rho_r);
    fftw_free(h_rho_k);

    delete[] h_kinetic_factor;
}

void FFTWRealTimeSplitSolver::prepare_transforms() {
    int nx = p->nx;
    int ny = p->ny;
    int nz = p->nz;
    int N  = nx * ny * nz;

    h_rho_r = (complex_type *)fftw_malloc(sizeof(complex_type) * N);
    h_rho_k = (complex_type *)fftw_malloc(sizeof(complex_type) * N);

    fftw_plan_with_nthreads(FFTW_N_THREADS);
    plan_fwd = fftw_plan_dft_3d(nx, ny, nz, h_rho_r, h_rho_k, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan_with_nthreads(FFTW_N_THREADS);
    plan_bwd = fftw_plan_dft_3d(nx, ny, nz, h_rho_k, h_rho_r, FFTW_BACKWARD, FFTW_MEASURE);
    OutputFormatter::printInfo("Planned FFTW with " + std::to_string(FFTW_N_THREADS) +
                               " threads (RTSplitSolver).");
}
