#ifndef GPE_SOLVER_HPP
#define GPE_SOLVER_HPP

#include "parameters/parameters.hpp"

#include "context/context.hpp"
#include "manager/sim_mediator.hpp"
#include "solver/fft_solvers/fft_export.hpp"
#include "solver/solver_data/cpu_solver_data.hpp"

#include <chrono>
#include <memory>

/*! Solver of time dependent Gross Pitaevski equation.
 *
 */
class AbstractGrossPitaevskiSolver {
public:
    AbstractGrossPitaevskiSolver(AbstractSimulationMediator *mediator);
    void solve();
    void initialize();

    void load_buffer(const wavefunction_t &);
    void load_pote(const potential_t &);

protected:
    // time of last checkpoint
    std::chrono::time_point<std::chrono::steady_clock> iter_time_ms;
    PhysicalParameters *params;

    energies_container_t enes;
    energies_t ene;

    //! current norm of wavefunction
    double xnorma;

    std::unique_ptr<AbstractPoissonSolver> poisson_solver;
    std::unique_ptr<AbstractRealTimeSplitSolver> rt_split_solver;
    AbstractSimulationMediator *p_mediator;
    SimulationContext *p_sctx;

    void calc_initial_state();
    void calc_evolution();

    void free_potential_well();
    void imag_time_iter();
    void real_time_iter();
    void summarize_imag_iter();
    void summarize_real_iter();

    //! CPU data buffer for file saving and program integration
    std::unique_ptr<CPUSolverData> buf_data;

    //! numerical methods
    virtual void calc_energy()                  = 0;
    virtual void init_containers()              = 0;
    virtual void calc_fi3d()                    = 0;
    virtual void calc_norm()                    = 0;
    virtual void normalize()                    = 0;
    virtual void imag_iter_linear_step()        = 0;
    virtual void imag_iter_nonlinear_step()     = 0;
    virtual void real_fft_potential_half_step() = 0;
    virtual void real_fft_kinetic_step()        = 0;

    //! prepare fft transformers
    virtual void prepare_fft() = 0;

    //! Data transfer between solver and base class data.
    virtual void export_data() = 0; // derived -> base
    virtual void import_data() = 0; // base -> derived
    virtual void import_pote() = 0; // base -> derived
};

#define MEASURE_TIME(func, ...)                                                                    \
    [&]() {                                                                                        \
        const auto MEASURE_TIME_START = std::chrono::high_resolution_clock::now();                 \
        func(__VA_ARGS__);                                                                         \
        const auto MEASURE_TIME_END = std::chrono::high_resolution_clock::now();                   \
        const std::chrono::duration<double, std::milli> MEASURE_TIME_DURATION =                    \
            MEASURE_TIME_END - MEASURE_TIME_START;                                                 \
        std::cout << "[MEASURE_TIME] " << #func << " took " << MEASURE_TIME_DURATION.count()       \
                  << " ms" << std::endl;                                                           \
    }()

#endif
