#include "solver/solver.hpp"
#include "output.hpp"
#include "parameters/parameters.hpp"
#include "solver/fft_rt_split_solver.hpp"
#include "solver/numerical_params.hpp"
#include "units.hpp"
#include <chrono>

double NumericalParameters::real_time_dt = 1.00e10;
double NumericalParameters::imag_time_dt = 1.25e11;

GrossPitaevskiSolver::GrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : params(PhysicalParameters::getInstance())
    , poisson_solver(std::make_unique<PoissonSolver>())
    , rt_split_solver(std::make_unique<RealTimeSplitSolver>())
    , p_mediator(mediator)
    , p_sctx(SimulationContext::getInstance()) {
}

void GrossPitaevskiSolver::initialize() {
    init_containers();

    poisson_solver->prepare(&cpsi, &fi3d, &pote);
    rt_split_solver->prepare(&cpsi, &fi3d, &pote);

    calc_norm();
    normalize();
}

void GrossPitaevskiSolver::solve() {
    switch (params->calc_strategy.type) {
    case CalcStrategy::Type::IMAGINARY_TIME:
        calc_initial_state();

        p_mediator->save_initial_state(cpsi);
        break;
    case CalcStrategy::Type::REAL_TIME:
        free_potential_well();
        calc_evolution();
        break;
    case CalcStrategy::Type::FULL:
        calc_initial_state();
        p_mediator->save_initial_state(cpsi);
        free_potential_well();
        calc_evolution();
        break;
    case CalcStrategy::Type::SPEED_TEST:
        run_speed_test();
        break;
    }
}

void GrossPitaevskiSolver::calc_initial_state() {
    OutputFormatter::printInfo("Starting imaginary time evolution");

    for (size_t iter = 1; iter <= params->iter_imag; iter++) {
        imag_time_iter();
        calc_energy();
    }

    OutputFormatter::printInfo("Imaginary time evolution completed");

    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Minimized energy [meV]: ",
                                       UnitConverter::ene_au_to_meV(ene.e_total));
    OutputFormatter::printBorderLine();

    p_mediator->save_energies(enes);
}

void GrossPitaevskiSolver::calc_evolution() {
    OutputFormatter::printInfo("Starting real time evolution");

    for (size_t iter = 1; iter <= params->iter_real; iter++) {
        real_time_iter();
        calc_energy();

        if (iter % 1000 == 0) {
            p_mediator->save_checkpoint(cpsi);
        }
    }

    OutputFormatter::printInfo("Real time evolution completed");
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Finished on energy [meV]: ",
                                       UnitConverter::ene_au_to_meV(ene.e_total));
    OutputFormatter::printBorderLine();

    p_mediator->save_data(cpsi);
    p_mediator->save_energies(enes);
}

void GrossPitaevskiSolver::imag_time_iter() {
    calc_fi3d();
    imag_iter_linear_step();
    imag_iter_nonlinear_step();
    calc_norm();
    normalize();

    //{
    //    auto start = std::chrono::high_resolution_clock::now();
    //    calc_fi3d();
    //    auto end = std::chrono::high_resolution_clock::now();
    //    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //    std::cout << "calc_fi3d: " << duration.count() << " ms\n";
    //}

    //{
    //    auto start = std::chrono::high_resolution_clock::now();
    //    imag_iter_linear_step();
    //    auto end = std::chrono::high_resolution_clock::now();
    //    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //    std::cout << "imag_iter_linear_step: " << duration.count() << " ms\n";
    //}

    //{
    //    auto start = std::chrono::high_resolution_clock::now();
    //    imag_iter_nonlinear_step();
    //    auto end = std::chrono::high_resolution_clock::now();
    //    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //    std::cout << "imag_iter_nonlinear_step: " << duration.count() << " ms\n";
    //}

    //{
    //    auto start = std::chrono::high_resolution_clock::now();
    //    calc_norm();
    //    normalize();
    //    auto end = std::chrono::high_resolution_clock::now();
    //    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //    std::cout << "calc_norm: " << duration.count() << " ms\n";
    //}
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
                                              cpsi(i, j, k) * params->w_15);
            }
        }
    }

    cpsi = cpsii;
}

void GrossPitaevskiSolver::real_time_iter() {

    real_fft_potential_half_step();
    calc_fi3d();
    real_fft_kinetic_step();
    real_fft_potential_half_step();

    //    // Timing for real_fft_potential_half_step()
    //    auto start1 = std::chrono::high_resolution_clock::now();
    //    real_fft_potential_half_step();
    //    auto end1 = std::chrono::high_resolution_clock::now();
    //    std::chrono::duration<double, std::milli> elapsed1 = end1 - start1;
    //    std::cout << "real_fft_potential_half_step: " << elapsed1.count() << " ms\n";
    //
    //    // Timing for calc_fi3d()
    //    auto start2 = std::chrono::high_resolution_clock::now();
    //    calc_fi3d();
    //    auto end2 = std::chrono::high_resolution_clock::now();
    //    std::chrono::duration<double, std::milli> elapsed2 = end2 - start2;
    //    std::cout << "calc_fi3d: " << elapsed2.count() << " ms\n";
    //
    //    // Timing for real_fft_kinetic_step()
    //    auto start3 = std::chrono::high_resolution_clock::now();
    //    real_fft_kinetic_step();
    //    auto end3 = std::chrono::high_resolution_clock::now();
    //    std::chrono::duration<double, std::milli> elapsed3 = end3 - start3;
    //    std::cout << "real_fft_kinetic_step: " << elapsed3.count() << " ms\n";
    //
    //    // Timing for real_fft_potential_half_step()
    //    auto start4 = std::chrono::high_resolution_clock::now();
    //    real_fft_potential_half_step();
    //    auto end4 = std::chrono::high_resolution_clock::now();
    //    std::chrono::duration<double, std::milli> elapsed4 = end4 - start4;
    //    std::cout << "real_fft_potential_half_step: " << elapsed4.count() << " ms\n";
}

void GrossPitaevskiSolver::init_containers() {
    init_potential();

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    cpsii.resize(nx, ny, nz);
    cpsi.resize(nx, ny, nz);
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
    double x = p_sctx->get_x(ix);
    double y = p_sctx->get_y(iy);
    double z = p_sctx->get_z(iz);

    double vx = -params->b * std::pow(x, 2) + params->aa * std::pow(x, 4);

    double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
    double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

    return vx + vy + vz;
}

double GrossPitaevskiSolver::pote_released_value(int ix, int iy, int iz) {
    double x = p_sctx->get_x(ix);
    double y = p_sctx->get_y(iy);
    double z = p_sctx->get_z(iz);

    double vx = params->aa * std::pow(x, 4);
    double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
    double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

    return vx + vy + vz;
}

double GrossPitaevskiSolver::pote_offset_value(int ix, int iy, int iz) {
    double x = p_sctx->get_x(ix);
    double y = p_sctx->get_y(iy);
    double z = p_sctx->get_z(iz);

    double vx = 0;
    if (x < 0)
        vx = params->aa * std::pow(x, 4);

    double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
    double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

    return vx + vy + vz;
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
                // pote(i, j, k) = pote_released_value(i, j, k);
                pote(i, j, k) = pote_offset_value(i, j, k);
            }
        }
    }
}

void GrossPitaevskiSolver::calc_energy() {
    ene.e_kin = 0.;
    ene.e_pot = 0.;
    ene.e_int = 0.;
    ene.e_ext = 0.;
    ene.e_bmf = 0.;

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
                grad_psi_x = -(cpsi(i + 1, j, k) + cpsi(i - 1, j, k) - 2. * cpsi(i, j, k)) /
                             (std::pow(params->dx, 2));
                grad_psi_y = -(cpsi(i, j + 1, k) + cpsi(i, j - 1, k) - 2. * cpsi(i, j, k)) /
                             (std::pow(params->dy, 2));
                grad_psi_z = -(cpsi(i, j, k + 1) + cpsi(i, j, k - 1) - 2. * cpsi(i, j, k)) /
                             (std::pow(params->dz, 2));
                ene.e_kin +=
                    ((grad_psi_x + grad_psi_y + grad_psi_z) * std::conj(cpsi(i, j, k))).real();

                // Potential energy
                ene.e_pot += pote(i, j, k) * std::norm(cpsi(i, j, k));

                // Interaction energy
                ene.e_int +=
                    0.5 * params->ggp11 * std::norm(cpsi(i, j, k)) * std::norm(cpsi(i, j, k));

                // Dipole-dipole interaction energy
                ene.e_ext +=
                    0.5 * params->cdd * fi3d(i, j, k) * std::norm(cpsi(i, j, k)) * params->n_atoms;
                ene.e_ext -= params->cdd / 3. * std::pow(std::norm(cpsi(i, j, k)), 2) / 2. *
                             std::pow(params->n_atoms, 2);

                // beyond mean-field energy
                ene.e_bmf += 2. / 5. * params->gamma * std::pow(std::abs(cpsi(i, j, k)), 5);
            }
        }
    }

    ene.e_kin *= params->get_dxdydz() / (2 * params->m) * params->n_atoms;
    ene.e_pot *= params->get_dxdydz() * params->n_atoms;
    ene.e_int *= params->get_dxdydz() * std::pow(params->n_atoms, 2);
    ene.e_ext *= params->get_dxdydz();
    ene.e_bmf *= params->get_dxdydz() * std::pow(params->n_atoms, 2.5);

    ene.sum();

    enes.emplace_back(ene);
}

void GrossPitaevskiSolver::run_speed_test() {
    OutputFormatter::printInfo("Starting speed test...");

    params->set_default_values();
    params->init_parameters();
    params->print();
    init_containers();

    poisson_solver->prepare(&cpsi, &fi3d, &pote);
    rt_split_solver->prepare(&cpsi, &fi3d, &pote);

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t iter = 0; iter < NumericalParameters::iter_imag_speed_test; iter++) {
        imag_time_iter();
    }

    auto end                              = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    OutputFormatter::printInfo("Imaginary time speed test completed.");
    OutputFormatter::printInfo("Total time for " +
                               std::to_string(NumericalParameters::iter_imag_speed_test) +
                               " iterations: " + std::to_string(elapsed.count()) + " seconds.");
    OutputFormatter::printInfo(
        "Average time per iteration: " +
        std::to_string(elapsed.count() /
                       static_cast<double>(NumericalParameters::iter_imag_speed_test)) +
        " seconds.");

    start = std::chrono::high_resolution_clock::now();

    for (size_t iter = 0; iter < NumericalParameters::iter_real_speed_test; iter++) {
        real_time_iter();
    }

    end     = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    OutputFormatter::printInfo("Real time speed test completed.");
    OutputFormatter::printInfo("Total time for " +
                               std::to_string(NumericalParameters::iter_real_speed_test) +
                               " iterations: " + std::to_string(elapsed.count()) + " seconds.");
    OutputFormatter::printInfo(
        "Average time per iteration: " +
        std::to_string(elapsed.count() /
                       static_cast<double>(NumericalParameters::iter_real_speed_test)) +
        " seconds.");
}

void GrossPitaevskiSolver::real_fft_potential_half_step() {
    int nx           = params->nx;
    int ny           = params->ny;
    int nz           = params->nz;
    double w         = params->n_atoms;
    double dt_factor = NumericalParameters::real_time_dt / 2.;

    double psi_re;
    double psi_im;
    double s;
    double c;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double v_ext   = pote(i, j, k);
                double density = std::norm(cpsi(i, j, k));

                double v_int = (params->ggp11 - params->cdd / 3) * density * w +
                               params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) * params->w_15;

                double total_potential = v_ext + params->cdd * fi3d(i, j, k) + v_int;

                // cpsi(i, j, k) *= std::exp(std::complex<double>(0.0, -dt_factor *
                // total_potential));
                psi_re = cpsi(i, j, k).real();
                psi_im = cpsi(i, j, k).imag();
                sincos(-dt_factor * total_potential, &s, &c);

                cpsi(i, j, k).imag(psi_re * s + psi_im * c);
                cpsi(i, j, k).real(psi_re * c - psi_im * s);
            }
        }
    }
}

void GrossPitaevskiSolver::real_fft_kinetic_step() {
    rt_split_solver->execute();
}

void GrossPitaevskiSolver::load_buffer(const wavefunction_t &wvf) {
    cpsii = wvf;
    cpsi  = wvf;

    calc_norm();
    normalize();
}
