#include "include/solver.hpp"
#include "include/numerical_params.hpp"
#include "include/params.hpp"
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <omp.h>

GrossPitaevskiSolver::GrossPitaevskiSolver()
    : params(PhysicalParameters::getInstance())
    , poisson_solver(std::make_unique<PoissonSolver>())
    , file_manager(std::make_unique<FileManager>()) {

    std::ofstream file("energy.dat");
    file.close();

    try {
        file_manager->load_params();
    } catch (const std::exception &e) {
        std::cerr << "Error loading params: " << e.what() << std::endl;
        std::cerr << "Using default parameters." << std::endl;
    }

    init_containers();
    file_manager->set_data_pointer(&cpsi);

    try {
        file_manager->load_initial_state();
    } catch (const std::exception &e) {
        std::cerr << "Error loading initial state: " << e.what() << std::endl;
        std::cerr << "Initializing with cosine function." << std::endl;
        init_with_cos();
    }

    poisson_solver->prepare(&cpsi, &fi3d);
    calc_norm();
    normalize();
}

void GrossPitaevskiSolver::solve() {
    calc_initial_state();
    save_initial_state_to_file();

    // load_initial_state_from_file()

    free_potential_well();
    calc_evolution();
}

void GrossPitaevskiSolver::calc_initial_state() {
    for (size_t iter = 0; iter < NumericalParameters::iter_imag_evo; iter++) {
        imag_time_iter();
        calc_energy();

        // if (iter % 100 == 0) {
        //     save_xy_cut_to_file(iter);
        //     save_ene_to_file(iter);
        // }
    }
}

void GrossPitaevskiSolver::calc_evolution() {
    for (size_t iter = 0; iter < NumericalParameters::iter_real_evo; iter++) {
        real_time_iter();
        calc_energy();

        if (iter % 100 == 0) {
            save_xy_cut_to_file(iter);
            save_ene_to_file(iter);
        }
    }
}

void GrossPitaevskiSolver::imag_time_iter() {
    calc_fi3d();

    imag_iter_linear_step();
    imag_iter_nonlinear_step();

    calc_norm();
    normalize();
}

void GrossPitaevskiSolver::imag_iter_linear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

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
}

void GrossPitaevskiSolver::imag_iter_nonlinear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    double w = params->n_atoms;
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
}

void GrossPitaevskiSolver::real_time_iter() {
    cpsin = cpsi;
    cpsii = cpsi;

    // calc_norm();
    calc_fi3d();

    for (int iter = 0; iter < 2; iter++) {
        real_iter_linear_step();
        real_iter_nonlinear_step();
        cpsin = cpsii;
    }

    cpsi = cpsin;
    // calc_norm();
    // normalize();
}

void GrossPitaevskiSolver::real_iter_linear_step() {
    std::complex<double> dt(0, -NumericalParameters::real_time_dt);
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

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
}

void GrossPitaevskiSolver::real_iter_nonlinear_step() {
    std::complex<double> dt(0, -NumericalParameters::real_time_dt);
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    double w = params->n_atoms;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                std::complex<double> c1 = ((params->ggp11 - params->cdd / 3) *
                                               std::norm(cpsi(i, j, k)) * cpsi(i, j, k) * w +
                                           params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
                                               cpsi(i, j, k) * std::pow(w, 1.5));

                std::complex<double> c2 = ((params->ggp11 - params->cdd / 3) *
                                               std::norm(cpsin(i, j, k)) * cpsin(i, j, k) * w +
                                           params->gamma * std::pow(std::abs(cpsin(i, j, k)), 3) *
                                               cpsin(i, j, k) * std::pow(w, 1.5));

                cpsii(i, j, k) = cpsii(i, j, k) + dt * (c1 + c2) / 2.;
            }
        }
    }
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
    normalize();
}

void GrossPitaevskiSolver::calc_fi3d() {
    poisson_solver->execute();
}

void GrossPitaevskiSolver::calc_norm() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    xnorma = 0.0;
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

void GrossPitaevskiSolver::calc_energy() {
    e_kin = 0.;
    e_pot = 0.;
    e_int = 0.;
    e_ext = 0.;
    e_bmf = 0.;

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    std::complex<double> grad_psi_x;
    std::complex<double> grad_psi_y;
    std::complex<double> grad_psi_z;

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                // Kinetic energy
                grad_psi_x = (cpsi(i + 1, j, k) - cpsi(i - 1, j, k)) / (2 * params->dx);
                grad_psi_y = (cpsi(i, j + 1, k) - cpsi(i, j - 1, k)) / (2 * params->dy);
                grad_psi_z = (cpsi(i, j, k + 1) - cpsi(i, j, k - 1)) / (2 * params->dz);
                e_kin += (std::norm(grad_psi_x) + std::norm(grad_psi_y) + std::norm(grad_psi_z));

                // Potential energy
                e_pot += pote(i, j, k) * std::norm(cpsi(i, j, k));

                // Interaction energy
                e_int += (params->ggp11 - params->cdd / 3) * std::norm(cpsi(i, j, k)) *
                         std::norm(cpsi(i, j, k));

                // Dipole-dipole interaction energy
                e_ext += params->cdd * fi3d(i, j, k) * std::norm(cpsi(i, j, k));

                // beyond mean-field energy
                e_bmf += params->gamma * std::pow(std::abs(cpsi(i, j, k)), 5);
            }
        }
    }

    e_kin *= params->get_dxdydz() / (2 * params->m);
    e_pot *= params->get_dxdydz();
    e_int *= params->get_dxdydz() * params->n_atoms;
    e_ext *= params->get_dxdydz();
    e_bmf *= params->get_dxdydz() * std::pow(params->n_atoms, 1.5);

    e_total = e_kin + e_pot + e_int + e_ext + e_bmf;
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

void GrossPitaevskiSolver::save_ene_to_file(int iter) {
    std::ofstream file("energy.dat", std::ios::app);

    file << iter << "\t" << e_kin << "\t" << e_pot << "\t" << e_int << "\t" << e_ext << "\t"
         << e_bmf << "\t" << e_total << std::endl;

    file.close();
}

void GrossPitaevskiSolver::load_initial_state_from_file() {

    calc_norm();
    normalize();
}

void GrossPitaevskiSolver::save_initial_state_to_file() {
}
