#include "include/solver.hpp"
#include "include/numerical_params.hpp"
#include "include/params.hpp"
#include <fftw3.h>

GrossPitaevskiSolver::GrossPitaevskiSolver()
    : params(PhysicalParameters::getInstance()) {

    init_containers();
    save_x_cut_to_file();
}

void GrossPitaevskiSolver::solve() {
    calc_initial_state();
}

void GrossPitaevskiSolver::calc_initial_state() {
    for (size_t iter = 0; iter < NumericalParameters::iter_imag_evo; iter++) {
        imag_time_iter();
    }

    save_x_cut_to_file();
    save_xy_cut_to_file();
}

void GrossPitaevskiSolver::imag_time_iter() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    // todo: push into double psi_norm();
    double xnorma = 0;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                xnorma += std::norm(cpsi(i, j, k));
            }
        }
    }

    xnorma *= params->get_dxdydz();

    calc_fi3d();

    //! WARNING: assumption dx = dy
    // \todo hide into imaginary step
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double v = pote(i, j, k);
                std::complex<double> c1 =
                    -0.5 / (params->m * std::pow(params->dx, 2)) *
                        (cpsi(i - 1, j, k) + cpsi(i + 1, j, k) - 2. * cpsi(i, j, k)) -
                    0.5 / (params->m * std::pow(params->dy, 2)) *
                        (cpsi(i, j - 1, k) + cpsi(i, j + 1, k) - 2. * cpsi(i, j, k)) -
                    0.5 / (params->m * std::pow(params->dz, 2)) *
                        (cpsi(i, j, k - 1) + cpsi(i, j, k + 1) - 2. * cpsi(i, j, k)) +
                    cpsi(i, j, k) * (v + params->cdd * fi3d(i, j, k));
                cpsii(i, j, k) = cpsi(i, j, k) - NumericalParameters::imag_time_dt * c1;
            }
        }
    }

    double w = params->n_atoms / xnorma;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                cpsii(i, j, k) = cpsii(i, j, k) -
                                 NumericalParameters::imag_time_dt *
                                     ((params->ggp11-params->cdd/3) * std::norm(cpsi(i, j, k)) * cpsi(i, j, k) * w +
                                      params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
                                          cpsi(i, j, k) * std::pow(w, 1.5));
            }
        }
    }

    cpsi = cpsii;

    calc_norm();
}

void GrossPitaevskiSolver::init_containers() {
    init_potential();

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    cpsii.resize(nx, ny, nz);
    cpsi.resize(nx, ny, nz);
    cpsin.resize(nx, ny, nz);

    fi3do.resize(nx, ny, nz);
    fi3d.resize(nx, ny, nz);

    init_with_cos();
}

void GrossPitaevskiSolver::init_potential() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    pote.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                pote(i, j, k) = pote_value(i, j, k);
            }
        }
    }
}

double GrossPitaevskiSolver::pote_value(int ix, int iy, int iz) {
    double x = params->get_x(ix);
    double y = params->get_y(iy);
    double z = params->get_z(iz);

    //! \todo: change wrl and wzl
    double vx = -params->b * std::pow(x, 2) + params->aa * std::pow(x, 4);
    double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->wrl, 2);
    double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->wzl, 2);

    return vx + vy + vz;
}

double GrossPitaevskiSolver::pote_released_value(int ix, int iy, int iz) {
    double x = params->get_x(ix);
    double y = params->get_y(iy);
    double z = params->get_z(iz);

    //! \todo: change wrl and wzl
    double vx = params->aa * std::pow(x, 4);
    double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->wrl, 2);
    double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->wzl, 2);

    return vx + vy + vz;
}

void GrossPitaevskiSolver::init_with_cos() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double x = params->get_x(i);
                double y = params->get_y(j);
                double z = params->get_z(k);

                double rrr = (static_cast<int>(params->nx / 2) * params->dx);

                double cos_x = std::cos(2 * M_PI * x / rrr);
                double cos_y = std::cos(2 * M_PI * y / rrr);
                double cos_z = std::cos(1 * M_PI * z / rrr);

                double val                = cos_x * cos_y * cos_z;
                std::complex<double> cval = std::complex<double>(val, 0.);
                cpsi(i, j, k)             = cval;
                cpsii(i, j, k)            = cval;
                cpsin(i, j, k)            = cval;
            }
        }
    }

    calc_norm();
}

void GrossPitaevskiSolver::calc_fi3d() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    size_t N = params->nx * params->ny * params->nz;

    // FFTW memory (aligned)
    fftw_complex *psi_r  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *psi_k  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *Vdip_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    // === 1. n(r) = |ψ|² ===
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                size_t idx    = (i * ny + j) * nz + k;
                double val    = std::norm(cpsi(i, j, k)); // |ψ|²
                psi_r[idx][0] = val;
                psi_r[idx][1] = 0.;
            }
        }
    }

    // === 2. FFT[ n(r) ] ===
    fftw_plan plan_fwd = fftw_plan_dft_3d(nx, ny, nz, psi_r, psi_k, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);

    double dkx = 2. * M_PI / (params->nx * params->dx);
    double dky = 2. * M_PI / (params->ny * params->dy);
    double dkz = 2. * M_PI / (params->nz * params->dz);

    for (int i = 0; i < nx; ++i) {
        double kx = (i <= (nx / 2)) ? i * dkx : (i - nx) * dkx;
        for (int j = 0; j < ny; ++j) {
            double ky = (j <= (ny / 2)) ? j * dky : (j - ny) * dky;
            for (int k = 0; k < nz; ++k) {
                double kz = (k <= (nz / 2)) ? k * dkz : (k - nz) * dkz;

                double k2  = kx * kx + ky * ky + kz * kz;
                size_t idx = (i * ny + j) * nz + k;

                if (k2 > 1e-12) {
                    Vdip_k[idx][0] = (4. * M_PI / 3.) * (3.0 * kz * kz / k2 - 1.0);
                    Vdip_k[idx][1] = 0.0;
                } else {
                    Vdip_k[idx][0] = 0.0;
                    Vdip_k[idx][1] = 0.0;
                }
            }
        }
    }

    // === 4. Pomnóż w przestrzeni pędów: psi_k *= Vdip_k ===
    for (size_t idx = 0; idx < N; ++idx) {
        std::complex<double> a(psi_k[idx][0], psi_k[idx][1]);
        std::complex<double> b(Vdip_k[idx][0], Vdip_k[idx][1]);
        std::complex<double> c = a * b;
        psi_k[idx][0]          = c.real();
        psi_k[idx][1]          = c.imag();
    }

    // === 5. Odwróć FFT ===
    fftw_plan plan_bwd = fftw_plan_dft_3d(nx, ny, nz, psi_k, psi_r, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_bwd);

    // === 6. Wynik (z normalizacją) ===
    double norm_factor = 1.0 / static_cast<double>(nx * ny * nz);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                size_t idx    = (i * ny + j) * nz + k;
                fi3d(i, j, k) = psi_r[idx][0] * norm_factor;
            }
        }
    }

    // === 7. Zwolnij pamięć i plany ===
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);

    fftw_free(psi_r);
    fftw_free(psi_k);
    fftw_free(Vdip_k);
}

void GrossPitaevskiSolver::calc_norm() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    double xnorma = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                xnorma += std::norm(cpsi(i, j, k));
            }
        }
    }

    xnorma *= params->get_dxdydz();

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                cpsi(i, j, k) /= std::sqrt(xnorma);
            }
        }
    }
}

void GrossPitaevskiSolver::save_xy_cut_to_file() {
    std::ofstream file("cut.dat");

    int z_zero_idx = params->nz / 2 + 1;

    for (int i = 0; i < params->nx; i++) {
        for (int j = 0; j < params->ny; j++) {
            file << params->get_x(i) << "\t" << params->get_y(j) << "\t"
                 << std::norm(cpsi(i, j, z_zero_idx)) << "\t" << fi3d(i, j, z_zero_idx) << "\t"
                 << pote(i, j, z_zero_idx) << "\n";
        }
        file << std::endl;
    }

    file << std::flush;
    file.close();
}

void GrossPitaevskiSolver::save_x_cut_to_file() {
    std::ofstream file("xcut.dat");

    int z_zero_idx = params->nz / 2 + 1;
    int y_zero_idx = params->ny / 2 + 1;

    for (int i = 0; i < params->nx; i++) {
        file << params->get_x(i) << "\t" << std::abs(cpsi(i, y_zero_idx, z_zero_idx)) << "\t"
             << fi3d(i, y_zero_idx, z_zero_idx) << "\t" << pote(i, y_zero_idx, z_zero_idx)
             << std::endl;
    }

    file.close();
}
