#ifndef REAL_TIME_SPLIT_SOLVER_HPP
#define REAL_TIME_SPLIT_SOLVER_HPP

#include "solver/fft_context.hpp"
#ifdef USE_CUDA
#include <cufft.h>
#include <cuda_runtime.h>
#else
#include <fftw3.h>
#endif

/*! class RealTimeSplitSolver.
*
* \brief class implements split step method for real time evolution.
*
* FFTW has been used in calculations.
*/
class RealTimeSplitSolver : public FFTContext {
public:
    RealTimeSplitSolver();
    ~RealTimeSplitSolver();

    void execute() override;
private:

    void prepare_transforms() override;
    void prepare_containers() override;

    complex_type *h_rho_r;
    complex_type *d_rho_r;
    double *kinetic_factor;
};

#endif
