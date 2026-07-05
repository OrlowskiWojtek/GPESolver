#ifndef FFTW_REAL_TIME_SPLIT_SOLVER_HPP
#define FFTW_REAL_TIME_SPLIT_SOLVER_HPP

#include "solver/fft_solvers/abstract_split_solver.hpp"

/*! class RealTimeSplitSolver.
*
* \brief class implements split step method for real time evolution.
*
* FFTW / cuFFT has been used in calculations.
*/
class FFTWRealTimeSplitSolver : public AbstractRealTimeSplitSolver {
public:
    FFTWRealTimeSplitSolver();
    ~FFTWRealTimeSplitSolver();

    void execute() override;
private:

    void prepare_transforms() override;
    void prepare_containers() override;
};

#endif
