#ifndef SIMULATION_MANAGER_HPP
#define SIMULATION_MANAGER_HPP

#include "filemanager/file_manager.hpp"
#include "solver/solver.hpp"
#include "manager/sim_mediator.hpp"
#include <memory>

class SimulationManager : public AbstractSimulationMediator{
public:
    SimulationManager();
    void initialize();
    void run_simulation();

    //! From AbstractSimualtionMediator
    void on_params_loaded() override;
    void on_data_loaded(const wavefunction_t&) override;

    void save_data(const wavefunction_t&) override;
    void save_energies(const energies_container_t&) override;

private:
    PhysicalParameters* params;
    std::unique_ptr<FileManager> file_manager;
    std::unique_ptr<GrossPitaevskiSolver> gpe_solver;

};

#endif
