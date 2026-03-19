#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "solver/fft_context.hpp"
#ifdef USE_CUDA
#include <cufft.h>
#include <cuda_runtime.h>
#else
#include <fftw3.h>
#endif

/*! class PoissonSolver.
*
* \brief class implements fast poisson solver for mesh.
*
* Solver implements R2C and C2R transformations for pure real input data.
* FFTW has been used in calculations.
*/
class PoissonSolver : public FFTContext {
public:
    PoissonSolver();
    ~PoissonSolver();

    void execute() override;
private:

    void prepare_transforms() override;
    void prepare_containers() override;

    complex_type* Vdip_k;

    complex_type *h_rho_k;
    complex_type *d_rho_k;

    double* d_rho_r;
    double* h_rho_r;
};

#endif
