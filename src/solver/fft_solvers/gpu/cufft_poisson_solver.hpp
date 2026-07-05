#ifndef CUFFT_POISSON_SOLVER_HPP
#define CUFFT_POISSON_SOLVER_HPP

#include "solver/fft_solvers/abstract_poisson_solver.hpp"

/*! class PoissonSolver.
 *
 * \brief class implements fast poisson solver for mesh.
 *
 * Solver implements R2C and C2R transformations for pure real input data.
 * FFTW has been used in calculations.
 */
class CUFFTPoissonSolver : public AbstractPoissonSolver {
public:
    CUFFTPoissonSolver();
    ~CUFFTPoissonSolver();

    void execute() override;

private:
    void prepare_transforms() override;
    void prepare_containers() override;
};

#endif
