#ifndef SIMULATION_MANAGER_HPP
#define SIMULATION_MANAGER_HPP

#include "filemanager/file_manager.hpp"
#include "solver/solver.hpp"
#include "manager/sim_mediator.hpp"
#include "initializer/initializer.hpp"
#include <memory>

class SimulationManager : public AbstractSimulationMediator{
public:
    SimulationManager();
    void initialize();
    void run_simulation();

    //! From AbstractSimualtionMediator
    void on_params_loaded() override;
    void on_data_loaded(const wavefunction_t&) override;
    void on_data_initialized(const wavefunction_t&) override;
    void on_pote_initialized(const potential_t&) override;
    void on_pote_changed(const potential_t&) override;

    void request_load_from_binary(wavefunction_t&) override;
    void request_load_from_text(wavefunction_t&) override;
    void request_free_potential() override;
    void request_cradle_potential() override;

    void save_data(const wavefunction_t&) override;
    void save_checkpoint(const wavefunction_t&) override;
    void save_initial_state(const wavefunction_t&) override;
    void save_energies(const energies_container_t&) override;

private:
    PhysicalParameters* params;
    SimulationContext* p_sctx;
    std::unique_ptr<FileManager> m_file_manager;
    std::unique_ptr<GrossPitaevskiSolver> m_gpe_solver;
    std::unique_ptr<DataInitializer> m_initializer;

    int checkpoint_counter = 0;
};

#endif
