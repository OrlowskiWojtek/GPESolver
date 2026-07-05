#ifndef ABSTRACT_POISSON_SOLVER_HPP
#define ABSTRACT_POISSON_SOLVER_HPP

#include "solver/fft_solvers/fft_context.hpp"

/*! class AbstractPoissonSolver.
 *
 * \brief class presents interface for abstract poisson solver used for dipolar potential calculations.
 */
class AbstractPoissonSolver : public FFTContext {
public:
    //AbstractPoissonSolver();
    //virtual ~AbstractPoissonSolver();

protected:
    complex_type *h_Vdip_k;
    complex_type *d_Vdip_k;

    complex_type *h_rho_k;
    complex_type *d_rho_k;

    complex_type *d_psi;

    real_type *d_rho_r;
    real_type *h_rho_r;
};

#endif
