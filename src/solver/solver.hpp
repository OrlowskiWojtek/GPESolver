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
    //! time of last checkpoint
    std::chrono::time_point<std::chrono::steady_clock> iter_time_ms;
    //! time of the beginning of calculations
    std::chrono::time_point<std::chrono::steady_clock> start_time_ms;
    PhysicalParameters *params;

    energies_container_t enes;
    energies_t ene;

    //! current norm of wavefunction
    double xnorma = 0;

    AbstractSimulationMediator *p_mediator;
    SimulationContext *p_sctx;

    void summarize_energies();
    void summarize_iter(int current_iter);

    //! CPU data buffer for file saving and program integration
    std::unique_ptr<CPUSolverData> buf_data;

    virtual void iterate()                     = 0;
    virtual void adjust(int iter)              = 0;
    virtual void finish()                      = 0;
    virtual const int iter_per_summary() const = 0;

    //! numerical methods
    virtual void init_containers() = 0;
    virtual void calc_energy()     = 0;
    virtual void calc_fi3d()       = 0;
    virtual void calc_norm()       = 0;
    virtual void normalize()       = 0;

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
