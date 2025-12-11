#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "include/fft_context.hpp"
#include <fftw3.h>

/*! class PoissonSolver.
*
* \brief class implements fast poisson solver for mesh.
*
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

    fftw_complex *Vdip_k;
};

#endif
