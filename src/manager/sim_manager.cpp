#include "manager/sim_manager.hpp"
#include "output.hpp"

SimulationManager::SimulationManager()
    : m_file_manager(std::make_unique<FileManager>(this))
    , m_gpe_solver(std::make_unique<GrossPitaevskiSolver>(this))
    , m_initializer(std::make_unique<DataInitializer>(this))
    , params(PhysicalParameters::getInstance()) {
}

void SimulationManager::initialize() {
    try {
        m_file_manager->load_params();
    } catch (const std::exception &e) {
        OutputFormatter::printWarning("Could not load parameters from file.");
        throw std::runtime_error(e.what());
    }

    m_initializer->initialize_wavefunction();
}

void SimulationManager::run_simulation() {
    m_gpe_solver->solve();   
}

void SimulationManager::save_data(const wavefunction_t &wvf) {
    std::string filename = "bec_wavefunction";
    m_file_manager->save_to_text_file(wvf, filename);
}

void SimulationManager::save_checkpoint(const wavefunction_t & wvf) {
    std::string filename = "wavefunction_" + std::to_string(checkpoint_counter);
    m_file_manager->save_to_text_file(wvf, filename);

    checkpoint_counter++;
}

void SimulationManager::save_initial_state(const wavefunction_t & wvf) {
    std::string filename = "initial_state";
    m_file_manager->save_to_text_file(wvf, filename);
}

void SimulationManager::save_energies(const energies_container_t & enes) {
    m_file_manager->save_energies(enes);
}

void SimulationManager::on_params_loaded() {
    params->init_parameters();
    params->print();

    // TODO: initialize context
    m_gpe_solver->initialize();
}

void SimulationManager::on_data_loaded(const wavefunction_t &t) {
    m_gpe_solver->load_buffer(t);
}

void SimulationManager::on_data_initialized(const wavefunction_t &t){

}
