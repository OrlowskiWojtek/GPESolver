#include "solver/solver.hpp"
#include "output.hpp"
#include "parameters/parameters.hpp"
#include "fft_solvers/fft_export.hpp"
#include "units.hpp"
#include <chrono>

double AbstractGrossPitaevskiSolver::real_time_dt = 1.00e10;
double AbstractGrossPitaevskiSolver::imag_time_dt = 1.25e11;

AbstractGrossPitaevskiSolver::AbstractGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : params(PhysicalParameters::getInstance())
    , p_mediator(mediator)
    , p_sctx(SimulationContext::getInstance()) {

#ifdef USE_CUDA
    poisson_solver = std::make_unique<CUFFTPoissonSolver>();
    rt_split_solver = std::make_unique<CUFFTRealTimeSplitSolver>();
#endif
}

void AbstractGrossPitaevskiSolver::initialize() {
    init_containers();

    poisson_solver->prepare(&cpsi, &fi3d, &pote);
    rt_split_solver->prepare(&cpsi, &fi3d, &pote);

    calc_norm();
    normalize();
}

void AbstractGrossPitaevskiSolver::solve() {
    iter_time_ms = std::chrono::steady_clock::now();
    
    switch (params->calc_strategy.type) {
    case CalcStrategy::Type::IMAGINARY_TIME:
        calc_initial_state();
        p_mediator->save_initial_state(cpsi);
        break;
    case CalcStrategy::Type::REAL_TIME:
        //free_potential_well();
        p_mediator->save_checkpoint(cpsi);
        calc_cradle();
        calc_evolution();
        break;
    case CalcStrategy::Type::FULL:
        calc_initial_state();
        p_mediator->save_initial_state(cpsi);
        calc_cradle();
        // free_potential_well();
        calc_evolution();
        break;
    case CalcStrategy::Type::SPEED_TEST:
        // TODO remove run_speed_test();
        break;
    }
}

void AbstractGrossPitaevskiSolver::calc_initial_state() {
    OutputFormatter::printInfo("Starting imaginary time evolution");

    for (size_t iter = 0; iter <= params->iter_imag; iter++) {
        imag_time_iter();
        calc_energy();

        if(iter % 100 == 0){
            p_mediator->save_checkpoint(cpsi);
            summarize_imag_iter();
        }
    }

    OutputFormatter::printInfo("Imaginary time evolution completed");

    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Minimized energy [meV]: ",
                                       UnitConverter::ene_au_to_meV(ene.e_total));
    OutputFormatter::printBorderLine();

    p_mediator->save_energies(enes);
}

void AbstractGrossPitaevskiSolver::calc_evolution() {
    OutputFormatter::printInfo("Starting real time evolution");

    for (size_t iter = 1; iter <= params->iter_real; iter++) {
        real_time_iter();
        calc_energy();

        if (iter % 1000 == 0) {
            p_mediator->save_checkpoint(cpsi);
            summarize_real_iter();
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

void AbstractGrossPitaevskiSolver::imag_time_iter() {
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

void AbstractGrossPitaevskiSolver::imag_iter_linear_step() {
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
                cpsii(i, j, k) = cpsi(i, j, k) - imag_time_dt * c1;
            }
        }
    }
}

void AbstractGrossPitaevskiSolver::imag_iter_nonlinear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    double w = params->n_atoms;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                cpsii(i, j, k) =
                    cpsii(i, j, k) - imag_time_dt *
                                         ((params->ggp11 - params->cdd / 3) *
                                              std::norm(cpsi(i, j, k)) * cpsi(i, j, k) * w +
                                          params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
                                              cpsi(i, j, k) * params->w_15);
            }
        }
    }

    cpsi = cpsii;
}

void AbstractGrossPitaevskiSolver::real_time_iter() {

    real_fft_potential_half_step();
    calc_fi3d();
    real_fft_kinetic_step();
    real_fft_potential_half_step();
}

void AbstractGrossPitaevskiSolver::calc_cradle() {
    OutputFormatter::printInfo("Changing potential to move one droplet");

    p_mediator->request_cradle_potential();
    for(int iter = 0; iter < 2000; iter++){
        real_time_iter();
    }
    p_mediator->request_free_potential();

    OutputFormatter::printInfo("Cradle has been moved");
}

void AbstractGrossPitaevskiSolver::free_potential_well() {
    p_mediator->request_free_potential();
}

void AbstractGrossPitaevskiSolver::load_buffer(const wavefunction_t &wvf) {
    cpsii = wvf;
    cpsi  = wvf;

    calc_norm();
    normalize();
}

void AbstractGrossPitaevskiSolver::load_pote(const potential_t &pote_initialized) {
    pote = pote_initialized;
}

void AbstractGrossPitaevskiSolver::summarize_imag_iter(){
    auto now = std::chrono::steady_clock::now();
    
    int time_elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - iter_time_ms).count();
    OutputFormatter::printInfo("Time per 100 iterations: " + std::to_string(time_elapsed_ms) + " ms");

    iter_time_ms = now;
}

void AbstractGrossPitaevskiSolver::summarize_real_iter(){
    auto now = std::chrono::steady_clock::now();
    
    int time_elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - iter_time_ms).count();
    OutputFormatter::printInfo("Time per 1000 iterations: " + std::to_string(time_elapsed_ms) + " ms");

    iter_time_ms = now;
}
