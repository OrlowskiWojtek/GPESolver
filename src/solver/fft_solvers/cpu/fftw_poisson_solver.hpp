#ifndef FFTW_POISSON_SOLVER_HPP
#define FFTW_POISSON_SOLVER_HPP

#include "solver/fft_solvers/abstract_poisson_solver.hpp"
#include "solver/fft_solvers/cpu/fftw_abstract_cpu_solver.hpp"

/*! class PoissonSolver.
 *
 * Solver uses FFTW R2C and C2R transformations for pure real input data.
 */
class FFTWPoissonSolver : public AbstractPoissonSolver, public FFTWAbstractCPUSolver {
public:
    FFTWPoissonSolver(wavefunction_t* psi, potential_t* fi3d, potential_t* pote);
    ~FFTWPoissonSolver();

    void execute() override;

private:
    void prepare_transforms() override;
    void prepare_containers() override;
};

#endif
