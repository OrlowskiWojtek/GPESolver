#include "manager/sim_manager.hpp"
#include "context/context.hpp"
#include "output.hpp"
#include "solver/solver_factory.hpp"
#include <regex>

SimulationManager::SimulationManager()
    : params(PhysicalParameters::getInstance())
    , p_sctx(SimulationContext::getInstance())
    , m_file_manager(std::make_unique<FileManager>(this))
    , m_initializer(std::make_unique<DataInitializer>(this)) {
}

void SimulationManager::initialize() {
    try {
        m_file_manager->load_params();
    } catch (const std::exception &e) {
        OutputFormatter::printWarning("Could not load parameters from file.");
        throw std::runtime_error(e.what());
    }
}

void SimulationManager::run_simulation() {
    m_gpe_solver->solve();
}

void SimulationManager::save_data(const wavefunction_t &wvf) {
    std::string filename = "bec_wavefunction";
    m_file_manager->save_to_binary_file(wvf, filename);
}

void SimulationManager::save_checkpoint(const wavefunction_t &wvf) {
    std::string filename = "wavefunction_" + std::to_string(checkpoint_counter);
    m_file_manager->save_to_binary_file(wvf, filename);

    checkpoint_counter++;
}

void SimulationManager::save_binary_checkpoint(
    const BinaryMixtureState &state) {
    save_binary_state(state, "wavefunction_");
    checkpoint_counter++;
}

void SimulationManager::save_initial_state(const wavefunction_t &wvf) {
    std::string filename = "initial_state";
    m_file_manager->save_to_text_file(wvf, filename);
}

void SimulationManager::save_binary_initial_state(
    const BinaryMixtureState &state) {
    m_file_manager->save_to_text_file(state.spec_a, "initial_state_a");
    m_file_manager->save_to_text_file(state.spec_b, "initial_state_b");
}

void SimulationManager::save_binary_state(const BinaryMixtureState &state,
                                          const std::string &filename) {
    m_file_manager->save_to_binary_file(
        state.spec_a, filename + "_a" + std::to_string(checkpoint_counter));
    m_file_manager->save_to_binary_file(
        state.spec_b, filename + "_b" + std::to_string(checkpoint_counter));
}

void SimulationManager::save_energies(const energies_container_t &enes) {
    m_file_manager->save_energies(enes);
    m_file_manager->save_energies_bin(enes);
}

void SimulationManager::on_params_loaded() {
    params->init_parameters();
    params->print();

    m_gpe_solver = SolverFactory::create(this);

    p_sctx->initialize();
    m_gpe_solver->initialize();

    PotentialInitConfig pote_conf;
    pote_conf.key = params->pote_key;
    pote_conf.add_random_noise = params->add_random_pote_noise;

    m_initializer->initialize_potential(pote_conf);

    if (params->calc_strategy.type ==
        CalcStrategy::Type::BINARY_IMAGINARY_TIME) {
        const auto *mix_params = BinaryMixtureParameters::getInstance();
        m_initializer->initialize_binary_wavefunctions(mix_params->init_a,
                                                       mix_params->init_b);
    } else {
        WavefunctionInitConfig wvf_conf;
        wvf_conf.key = params->wvf_key;
        wvf_conf.add_random_noise = params->add_random_noise;
        wvf_conf.load_filename = params->load_filename;
        wvf_conf.dd = params->dd;
        wvf_conf.n_gauss_max = params->n_gauss_max;
        wvf_conf.bec_droplets_x = params->bec_droplets_x;
        wvf_conf.bec_droplets_y = params->bec_droplets_y;
        wvf_conf.bec_droplets_z = params->bec_droplets_z;
        m_initializer->initialize_wavefunction(wvf_conf);
    }
}

void SimulationManager::on_data_loaded(const wavefunction_t &wvf) {
    m_gpe_solver->load_buffer(wvf);
}

void SimulationManager::on_data_initialized(const wavefunction_t &wvf) {
    m_gpe_solver->load_buffer(wvf);
}

void SimulationManager::on_binary_data_initialized(
    const BinaryMixtureState &state) {
    m_gpe_solver->load_binary_buffer(state);
}

void SimulationManager::on_pote_initialized(const potential_t &pote) {
    m_file_manager->save_pote_to_text_file(pote, "initial_potential");
    m_gpe_solver->load_pote(pote);
}

void SimulationManager::on_pote_changed(const potential_t &pote) {
    m_file_manager->save_pote_to_text_file(pote, "changed_potential");
    m_gpe_solver->load_pote(pote);
}

void SimulationManager::request_free_potential() {
    // Potential storage is currently owned by the initializer and solver.
}

void SimulationManager::request_load_from_text(const std::string &filename,
                                               wavefunction_t &wvf) {
    check_load_filename(filename);
    m_file_manager->load_from_text_file(filename,
                                        [&wvf](wavefunction_t &loaded_buffer) {
                                            if (loaded_buffer.size() != 0) {
                                                wvf = loaded_buffer;
                                            }
                                        });
}

void SimulationManager::request_load_from_binary(const std::string &filename,
                                                 wavefunction_t &wvf) {
    check_load_filename(filename);
    m_file_manager->load_from_binary_file(
        filename, [&wvf](wavefunction_t &loaded_buffer) {
            if (loaded_buffer.size() != 0) {
                wvf = loaded_buffer;
            }
        });
}

void SimulationManager::check_load_filename() {
    check_load_filename(params->load_filename);
}

void SimulationManager::check_load_filename(const std::string &filename) {
    if (filename.find("wavefunction_") == std::string::npos) {
        return;
    }

    std::regex pattern("wavefunction_(\\d+)");
    std::smatch match;

    if (std::regex_search(filename, match, pattern)) {
        std::string number_str = match[1].str();        // grabujemy cyfrę
        int number             = std::stoi(number_str); // konwertujemy na int

        checkpoint_counter = number + 1;
    }
}
