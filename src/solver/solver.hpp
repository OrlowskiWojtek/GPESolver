#ifndef GPE_SOLVER_HPP
#define GPE_SOLVER_HPP

#include "parameters/parameters.hpp"
#include "solver/fft_rt_split_solver.hpp"
#include "solver/fft_poisson_solver.hpp"
#include "manager/sim_mediator.hpp"

#include "context/context.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <memory>

/*! Solver of time dependent Gross Pitaevski equation.
 *
 *
 */
class GrossPitaevskiSolver {
public:
    GrossPitaevskiSolver(AbstractSimulationMediator* mediator);
    void solve();

private:
    PhysicalParameters *params;

    //! Containers
    //! Wavefunction of bec - copy for calculations.
    StdMat3D<std::complex<double>> cpsii;

    //! Wavefunction of bec - copy for calculations.
    StdMat3D<std::complex<double>> cpsi;

    //! Map of dipole-dipole potential - copy
    StdMat3D<double> fi3d;

    //! Map of external potential
    StdMat3D<double> pote;

    BECEnergies _enes;

    // current norm of wavefunction
    double xnorma;

    std::unique_ptr<PoissonSolver> poisson_solver;
    std::unique_ptr<RealTimeSplitSolver> rt_split_solver;
    AbstractSimulationMediator* mediator;

    void calc_initial_state();
    void calc_evolution();
    void calc_energy();
    void run_speed_test();

    void init_containers();
    void init_with_cos();
    void init_with_gauss();
    void init_with_multiple_gauss();
    void init_potential();
    void free_potential_well();
    void imag_time_iter();
    void real_time_iter();

    void calc_fi3d();
    void calc_norm();
    void normalize();
    void imag_iter_linear_step();
    void imag_iter_nonlinear_step();
    void real_fft_potential_half_step();
    void real_fft_kinetic_step();

    double pote_value(int ix, int iy, int iz);
    double pote_released_value(int ix, int iy, int iz);
};

#endif
