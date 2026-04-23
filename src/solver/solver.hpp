#ifndef GPE_SOLVER_HPP
#define GPE_SOLVER_HPP

#include "parameters/parameters.hpp"

#include "solver/fft_solvers/fft_export.hpp"
#include "manager/sim_mediator.hpp"
#include "context/context.hpp"

#include <chrono>
#include <memory>

/*! Solver of time dependent Gross Pitaevski equation.
 *
 *
 */
class GrossPitaevskiSolver {
public:
    GrossPitaevskiSolver(AbstractSimulationMediator* mediator);
    void solve();

    void initialize();
    void load_buffer(const wavefunction_t&);
    void load_pote(const potential_t&);
private:
    // time of last checkpoint
    std::chrono::time_point<std::chrono::steady_clock> iter_time_ms;

    PhysicalParameters *params;

    //! Containers
    //! Wavefunction of bec + copy for calculations.
    wavefunction_t cpsii;
    wavefunction_t cpsi;

    //! Map of dipole-dipole potential
    potential_t fi3d;

    //! Map of external potential
    potential_t pote;

    energies_container_t enes;
    energies_t ene;

    // current norm of wavefunction
    double xnorma;

    static double imag_time_dt;
    static double real_time_dt;

    std::unique_ptr<AbstractPoissonSolver> poisson_solver;
    std::unique_ptr<AbstractRealTimeSplitSolver> rt_split_solver;
    AbstractSimulationMediator* p_mediator;
    SimulationContext* p_sctx;

    void calc_initial_state();
    void calc_evolution();
    void calc_cradle();
    void calc_energy();
    void run_speed_test();

    void init_containers();
    void free_potential_well();
    void imag_time_iter();
    void real_time_iter();
    void summarize_imag_iter();
    void summarize_real_iter();

    void calc_fi3d();
    void calc_norm();
    void normalize();
    void imag_iter_linear_step();
    void imag_iter_nonlinear_step();
    void real_fft_potential_half_step();
    void real_fft_kinetic_step();
};

#endif
