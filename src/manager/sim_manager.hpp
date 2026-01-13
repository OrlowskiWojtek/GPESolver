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

    void save_data(const wavefunction_t&) override;
    void save_checkpoint(const wavefunction_t&) override;
    void save_initial_state(const wavefunction_t&) override;
    void save_energies(const energies_container_t&) override;

private:
    PhysicalParameters* params;
    std::unique_ptr<FileManager> m_file_manager;
    std::unique_ptr<GrossPitaevskiSolver> m_gpe_solver;
    std::unique_ptr<DataInitializer> m_initializer;

    int checkpoint_counter = 0;
};

#endif
