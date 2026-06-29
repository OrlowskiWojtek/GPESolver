#include "solver/solver.hpp"
#include "fft_solvers/fft_export.hpp"
#include "output.hpp"
#include "parameters/parameters.hpp"
#include "units.hpp"
#include <chrono>

AbstractGrossPitaevskiSolver::AbstractGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : params(PhysicalParameters::getInstance())
    , p_mediator(mediator)
    , p_sctx(SimulationContext::getInstance())
    , buf_data(std::make_unique<CPUSolverData>()) {
}

void AbstractGrossPitaevskiSolver::initialize() {
    buf_data->allocate(params->nx, params->ny, params->nz);
    init_containers();

    prepare_fft();

    calc_norm();
    normalize();
}

void AbstractGrossPitaevskiSolver::solve() {
    iter_time_ms = std::chrono::steady_clock::now();

    switch (params->calc_strategy.type) {
    case CalcStrategy::Type::IMAGINARY_TIME:
        calc_initial_state();

        p_mediator->save_initial_state(buf_data->cpsi);
        break;
    case CalcStrategy::Type::REAL_TIME:
        // free_potential_well();
        calc_evolution();
        break;
    case CalcStrategy::Type::FULL:
        calc_initial_state();

        p_mediator->save_initial_state(buf_data->cpsi);
        free_potential_well();
        calc_evolution();
        break;
    }
}

void AbstractGrossPitaevskiSolver::calc_initial_state() {
    OutputFormatter::printInfo("Starting imaginary time evolution");

    for (size_t iter = 0; iter <= params->iter_imag; iter++) {
        imag_time_iter();

        if (iter % 100 == 0) {
            export_data();
            p_mediator->save_checkpoint(buf_data->cpsi);
            summarize_imag_iter();
            calc_energy();
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

    for (size_t iter = 0; iter <= params->iter_real; iter++) {
        real_time_iter();

        if (iter % 1000 == 0) {
            export_data();
            p_mediator->save_checkpoint(buf_data->cpsi);
            summarize_real_iter();
            calc_energy();
        }
    }

    OutputFormatter::printInfo("Real time evolution completed");
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Finished on energy [meV]: ",
                                       UnitConverter::ene_au_to_meV(ene.e_total));
    OutputFormatter::printBorderLine();

    p_mediator->save_data(buf_data->cpsi);
    p_mediator->save_energies(enes);
}

void AbstractGrossPitaevskiSolver::imag_time_iter() {
    calc_fi3d();
    imag_iter_linear_step();
    imag_iter_nonlinear_step();
    calc_norm();
    normalize();
}

void AbstractGrossPitaevskiSolver::real_time_iter() {

    real_fft_potential_half_step();
    calc_fi3d();
    real_fft_kinetic_step();
    real_fft_potential_half_step();
}

void AbstractGrossPitaevskiSolver::free_potential_well() {
    p_mediator->request_free_potential();
}

void AbstractGrossPitaevskiSolver::load_buffer(const wavefunction_t &wvf) {
    buf_data->cpsi  = wvf;
    buf_data->cpsii = wvf;

    import_data();
}

void AbstractGrossPitaevskiSolver::load_pote(const potential_t &pote_initialized) {
    buf_data->pote = pote_initialized;

    import_pote();
}

void AbstractGrossPitaevskiSolver::summarize_imag_iter() {
    auto now = std::chrono::steady_clock::now();

    int time_elapsed_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - iter_time_ms).count();
    OutputFormatter::printInfo("Time per 100 iterations: " + std::to_string(time_elapsed_ms) +
                               " ms");

    iter_time_ms = now;
}

void AbstractGrossPitaevskiSolver::summarize_real_iter() {
    auto now = std::chrono::steady_clock::now();

    int time_elapsed_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - iter_time_ms).count();
    OutputFormatter::printInfo("Time per 1000 iterations: " + std::to_string(time_elapsed_ms) +
                               " ms");

    iter_time_ms = now;
}
