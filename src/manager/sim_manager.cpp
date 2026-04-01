#include "manager/sim_manager.hpp"
#include "context/context.hpp"
#include "output.hpp"

SimulationManager::SimulationManager()
    : params(PhysicalParameters::getInstance())
    , p_sctx(SimulationContext::getInstance())
    , m_file_manager(std::make_unique<FileManager>(this))
    , m_gpe_solver(std::make_unique<GrossPitaevskiSolver>(this))
    , m_initializer(std::make_unique<DataInitializer>(this)) {
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

void SimulationManager::save_checkpoint(const wavefunction_t &wvf) {
    std::string filename = "wavefunction_" + std::to_string(checkpoint_counter);
    m_file_manager->save_to_text_file(wvf, filename);

    checkpoint_counter++;
}

void SimulationManager::save_initial_state(const wavefunction_t &wvf) {
    std::string filename = "initial_state";
    m_file_manager->save_to_text_file(wvf, filename);
}

void SimulationManager::save_energies(const energies_container_t &enes) {
    m_file_manager->save_energies(enes);
}

void SimulationManager::on_params_loaded() {
    params->init_parameters();
    params->print();

    p_sctx->initialize();
    m_gpe_solver->initialize();
}

void SimulationManager::on_data_loaded(const wavefunction_t &wvf) {
    m_gpe_solver->load_buffer(wvf);
}

void SimulationManager::on_data_initialized(const wavefunction_t &wvf) {
    m_gpe_solver->load_buffer(wvf);
}

void SimulationManager::request_load_from_text(wavefunction_t &wvf) {
    m_file_manager->load_from_text_file(params->load_filename,
                                        [&wvf](wavefunction_t &loaded_buffer) {
                                            if (loaded_buffer.size() != 0) {
                                                wvf = loaded_buffer;
                                            }
                                        });
}

void SimulationManager::request_load_from_binary(wavefunction_t &wvf) {
    m_file_manager->load_from_binary_file(params->load_filename,
                                          [&wvf](wavefunction_t &loaded_buffer) {
                                              if (loaded_buffer.size() != 0) {
                                                  wvf = loaded_buffer;
                                              }
                                          });
}
