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
        OutputFormatter::printWarning(e.what());
        OutputFormatter::printInfo("Using default parameters.");
        // rather kill

        return;
        params->set_default_values();
    }

    params->init_parameters();
    params->print();

    if (params->load_initial_state) {
        try {
            // TODO: add interface
            // file_manager->load_from_text_file("text_file");
        } catch (const std::exception &e) {
            OutputFormatter::printWarning("Could not load initial state from file.");
            OutputFormatter::printWarning(e.what());

            return;
            // init_with_gauss();
        }
    } else {
        // init_with_multiple_gauss();
    }
}

//TODO: implement connections
void SimulationManager::run_simulation() {
}

void SimulationManager::save_data(const wavefunction_t &) {
}

void SimulationManager::save_energies(const energies_container_t &) {
}

void SimulationManager::on_params_loaded() {
}

void SimulationManager::on_data_loaded(const wavefunction_t &) {
}
