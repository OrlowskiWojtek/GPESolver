#ifndef REAL_TIME_SPLIT_SOLVER_HPP
#define REAL_TIME_SPLIT_SOLVER_HPP

#include "solver/fft_context.hpp"
#include <fftw3.h>

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

    double *kinetic_factor;
};

#endif
