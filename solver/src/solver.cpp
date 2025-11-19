#include "include/solver.hpp"
#include "include/numerical_params.hpp"
#include "include/params.hpp"
#include <fftw3.h>
#include <fstream>
#include <iostream>

GrossPitaevskiSolver::GrossPitaevskiSolver()
    : params(PhysicalParameters::getInstance())
    , poisson_solver(std::make_unique<PoissonSolver>()) {

    init_containers();
    poisson_solver->prepare(&cpsi, &fi3d);
    calc_norm();
}

void GrossPitaevskiSolver::solve() {
    calc_initial_state();

    // free_potential_well();
    // calc_evolution();
}

void GrossPitaevskiSolver::calc_initial_state() {
    for (size_t iter = 0; iter < NumericalParameters::iter_imag_evo; iter++) {
        imag_time_iter();

        if (iter % 100 == 0) {
            save_xy_cut_to_file(iter);
        }
    }
}

void GrossPitaevskiSolver::calc_evolution() {
    for (size_t iter = 0; iter < NumericalParameters::iter_real_evo; iter++) {
        real_time_iter();
        if (iter % 1000 == 0) {
            std::cout << "Saving to file" << iter << std::endl;
            save_xy_cut_to_file(iter);
        }
    }
}

void GrossPitaevskiSolver::imag_time_iter() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    // todo: push into double psi_norm();
    calc_norm();
    // xnorma = 0;
    // for (int i = 1; i < nx - 1; i++) {
    //     for (int j = 1; j < ny - 1; j++) {
    //         for (int k = 1; k < nz - 1; k++) {
    //             xnorma += std::norm(cpsi(i, j, k));
    //         }
    //     }
    // }

    // xnorma *= params->get_dxdydz();

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
    // hide into something
    double w = params->n_atoms / xnorma;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                cpsii(i, j, k) =
                    cpsii(i, j, k) - NumericalParameters::imag_time_dt *
                                         ((params->ggp11 - params->cdd / 3) *
                                              std::norm(cpsi(i, j, k)) * cpsi(i, j, k) * w +
                                          params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
                                              cpsi(i, j, k) * std::pow(w, 1.5));
            }
        }
    }

    cpsi = cpsii;
    calc_norm();
    normalize();
}

void GrossPitaevskiSolver::real_time_iter() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    cpsin = cpsi;
    cpsii = cpsi;

    // todo: push into double psi_norm();
    xnorma = 0;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                xnorma += std::norm(cpsi(i, j, k));
            }
        }
    }

    xnorma *= params->get_dxdydz();

    calc_fi3d();
    std::complex<double> dt(0, -NumericalParameters::real_time_dt);

    for (int iter = 0; iter < 2; iter++) {
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
                        cpsi(i, j, k) * (v + params->cdd * fi3d(i, j, k)); // need old fi3d here

                    std::complex<double> c2 =
                        -0.5 / (params->m * std::pow(params->dx, 2)) *
                            (cpsin(i - 1, j, k) + cpsin(i + 1, j, k) - 2. * cpsin(i, j, k)) -
                        0.5 / (params->m * std::pow(params->dy, 2)) *
                            (cpsin(i, j - 1, k) + cpsin(i, j + 1, k) - 2. * cpsin(i, j, k)) -
                        0.5 / (params->m * std::pow(params->dz, 2)) *
                            (cpsin(i, j, k - 1) + cpsin(i, j, k + 1) - 2. * cpsin(i, j, k)) +
                        cpsin(i, j, k) * (v + params->cdd * fi3d(i, j, k));

                    cpsii(i, j, k) = cpsi(i, j, k) + dt * (c1 + c2) / 2.;
                }
            }
        }

        double w = params->n_atoms / xnorma;
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int k = 1; k < nz - 1; k++) {
                    std::complex<double> c1 =
                        ((params->ggp11 - params->cdd / 3) * std::norm(cpsi(i, j, k)) *
                             cpsi(i, j, k) * w +
                         params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) * cpsi(i, j, k) *
                             std::pow(w, 1.5));

                    std::complex<double> c2 =
                        ((params->ggp11 - params->cdd / 3) * std::norm(cpsin(i, j, k)) *
                             cpsin(i, j, k) * w +
                         params->gamma * std::pow(std::abs(cpsin(i, j, k)), 3) * cpsin(i, j, k) *
                             std::pow(w, 1.5));

                    cpsii(i, j, k) = cpsii(i, j, k) + dt * (c1 + c2) / 2.;
                }
            }
        }

        cpsin = cpsii;
    }

    cpsi = cpsin;
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
    poisson_solver->execute(xnorma);
}

void GrossPitaevskiSolver::calc_norm() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    xnorma = 0.;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                xnorma += std::norm(cpsi(i, j, k));
            }
        }
    }

    xnorma *= params->get_dxdydz();
}

void GrossPitaevskiSolver::normalize() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                cpsi(i, j, k) /= std::sqrt(xnorma);
            }
        }
    }
}

void GrossPitaevskiSolver::save_xy_cut_to_file(int iter) {
    std::ofstream file("cut" + std::to_string(iter) + ".dat");

    int z_zero_idx = params->nz / 2 + 1;

    for (int i = 0; i < params->nx; i++) {
        for (int j = 0; j < params->ny; j++) {
            file << params->get_x(i) << "\t" << params->get_y(j) << "\t"
                 << std::norm(cpsi(i, j, z_zero_idx)) << "\t" << fi3d(i, j, z_zero_idx) << "\t"
                 << pote(i, j, z_zero_idx) << "\n";
        }
        // file << std::endl;
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

void GrossPitaevskiSolver::free_potential_well() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    pote.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                pote(i, j, k) = pote_released_value(i, j, k);
            }
        }
    }
}
