#ifndef FFTW_POISSON_SOLVER_HPP
#define FFTW_POISSON_SOLVER_HPP

#include "solver/fft_solvers/abstract_poisson_solver.hpp"

/*! class PoissonSolver.
 *
 * Solver uses FFTW R2C and C2R transformations for pure real input data.
 */
class FFTWPoissonSolver : public AbstractPoissonSolver {
public:
    FFTWPoissonSolver();
    ~FFTWPoissonSolver();

    void execute() override;

private:
    void prepare_transforms() override;
    void prepare_containers() override;
};

#endif
