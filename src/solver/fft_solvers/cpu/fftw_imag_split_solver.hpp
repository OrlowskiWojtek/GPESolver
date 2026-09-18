#ifndef FFTW_IMAG_TIME_SPLIT_SOLVER_HPP
#define FFTW_IMAG_TIME_SPLIT_SOLVER_HPP

#include "solver/fft_solvers/abstract_split_solver.hpp"
#include "solver/fft_solvers/cpu/fftw_abstract_cpu_solver.hpp"

/*! class ImagTimeSplitSolver.
*
* \brief class implements split step method for imag time evolution.
*
* FFTW / cuFFT has been used in calculations.
*/
class FFTWImagTimeSplitSolver : public AbstractRealTimeSplitSolver, public FFTWAbstractCPUSolver {
public:
    FFTWImagTimeSplitSolver(wavefunction_t* psi, potential_t* fi3d);
    ~FFTWImagTimeSplitSolver();

    void execute() override;
private:

    void prepare_transforms() override;
    void prepare_containers() override;
};

#endif
