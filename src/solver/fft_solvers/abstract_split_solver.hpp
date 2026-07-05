#ifndef REAL_TIME_SPLIT_SOLVER_HPP
#define REAL_TIME_SPLIT_SOLVER_HPP

#include "solver/fft_solvers/fft_context.hpp"

/*! class RealTimeSplitSolver.
*
* \brief class implements split step method for real time evolution.
* FFTW / cuFFT has been used in calculations.
*/
class AbstractRealTimeSplitSolver : public FFTContext {
public:

protected:
    complex_type *h_rho_r;
    complex_type *d_rho_r;

    complex_type *h_rho_k;
    complex_type *d_rho_k;

    real_type *h_kinetic_factor;
    real_type *d_kinetic_factor;
};

#endif
