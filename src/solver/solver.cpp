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
    iter_time_ms  = std::chrono::steady_clock::now();
    start_time_ms = std::chrono::steady_clock::now();

    switch (params->calc_strategy.type) {
    case CalcStrategy::Type::IMAGINARY_TIME:
        calc_initial_state();

        export_data();
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
            if(params->save_imag_time_data) { 
                export_data();
                p_mediator->save_checkpoint(buf_data->cpsi);
            }
            calc_energy();
            summarize_energies();
            summarize_imag_iter(iter);
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
        if (!params->const_edd)
            params->update_edd(iter);

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
    calc_norm();
    normalize();
}

void AbstractGrossPitaevskiSolver::load_pote(const potential_t &pote_initialized) {
    buf_data->pote = pote_initialized;

    import_pote();
}

void AbstractGrossPitaevskiSolver::summarize_imag_iter(int current_iter) {
    if(current_iter == 0){
        return;
    }

    auto now = std::chrono::steady_clock::now();
    double frc = static_cast<double>(current_iter) /  
                 static_cast<double>(params->iter_imag);

    int time_elapsed_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - iter_time_ms).count();

    int time_elapsed_ms_from_start =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_ms).count();
    int time_predicted_ms =  time_elapsed_ms_from_start * ( 1. / frc - 1.);

    OutputFormatter::printInfo("Time per 100 iterations: " +
                               std::to_string(time_elapsed_ms) +
                               " ms | finished " +
                               std::to_string(frc * 100.) +
                               " % | predicted time [s]: " +
                               std::to_string(time_predicted_ms / 1000));

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

void AbstractGrossPitaevskiSolver::summarize_energies() {
    if (enes.size() < 2) {
        return;
    }

    energies_t &enes_last = *(enes.end() - 2);
    energies_t &enes_now  = enes.back();
    
    if(std::isnan(enes_now.e_total)){
        assert(!"NaN value detected, aborting calculations");
    }

    if(std::isinf(enes_now.e_total)){
        assert(!"inf value detected, aborting calculations");
    }

    //! dip-dip interaction
    double ext_diff = std::abs(enes_last.e_ext - enes_now.e_ext) / std::abs(enes_last.e_ext);
    //! lhy correction
    double bmf_diff = std::abs(enes_last.e_bmf - enes_now.e_bmf) / std::abs(enes_last.e_bmf);
    //! kinetic energy
    double kin_diff = std::abs(enes_last.e_kin - enes_now.e_kin) / std::abs(enes_last.e_kin);
    //! external potential energy
    double pot_diff = std::abs(enes_last.e_pot - enes_now.e_pot) / std::abs(enes_last.e_pot);
    //! contact interaction energy
    double int_diff = std::abs(enes_last.e_int - enes_now.e_int) / std::abs(enes_last.e_int);

    double tot_diff = (ext_diff) + (bmf_diff) + (kin_diff) + (pot_diff) + (int_diff);

    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Energies");
    OutputFormatter::printScientificRow<double>(
        {"dip-dip", "contact", "lhy", "kinetic", "potential", "total"},
        {enes_now.e_ext, enes_now.e_int, enes_now.e_bmf, enes_now.e_kin, enes_now.e_pot, enes_now.e_total});
    OutputFormatter::printBoxedMessage("Differences");
    OutputFormatter::printScientificRow<double>(
        {"dip-dip", "contact", "lhy", "kinetic", "potential", "total"},
        {ext_diff, int_diff, bmf_diff, kin_diff, pot_diff, tot_diff});
    OutputFormatter::printBorderLine();
}
