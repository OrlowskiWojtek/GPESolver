#include "manager/sim_manager.hpp"
#include "output.hpp"

SimulationManager::SimulationManager()
    : file_manager(std::make_unique<FileManager>(this))
    , gpe_solver(std::make_unique<GrossPitaevskiSolver>(this))
    , params(PhysicalParameters::getInstance()) {
}

void SimulationManager::initialize() {
    try {
        file_manager->load_params();
    } catch (const std::exception &e) {
        OutputFormatter::printWarning("Could not load parameters from file:");

        throw std::runtime_error(e.what());
        params->set_default_values();
    }

    if (params->load_initial_state) {
        try {
            // file_manager->load_from_text_file("text_file");
        } catch (const std::exception &e) {
            OutputFormatter::printWarning("Could not load initial state from file.");
            OutputFormatter::printWarning(e.what());
        }
    } else {
        // init_with_multiple_gauss();
    }
}

void SimulationManager::run_simulation() {
    gpe_solver->solve();   
}

void SimulationManager::save_data(const wavefunction_t &wvf) {
    std::string filename = "bec_wavefunction";
    file_manager->save_to_text_file(wvf, filename);
}

void SimulationManager::save_checkpoint(const wavefunction_t & wvf) {
    std::string filename = "wavefunction_" + std::to_string(checkpoint_counter);
    file_manager->save_to_text_file(wvf, filename);

    checkpoint_counter++;
}

void SimulationManager::save_initial_state(const wavefunction_t & wvf) {
    std::string filename = "initial_state";
    file_manager->save_to_text_file(wvf, filename);
}

void SimulationManager::save_energies(const energies_container_t & enes) {
    file_manager->save_energies(enes);
}

void SimulationManager::on_params_loaded() {
    params->init_parameters();
    params->print();

    gpe_solver->initialize();
}

void SimulationManager::on_data_loaded(const wavefunction_t &t) {
    gpe_solver->load_buffer(t);
}
