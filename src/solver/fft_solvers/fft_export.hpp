#include "solver/fft_solvers/fft_context.hpp"
#ifdef USE_CUDA
#include "solver/fft_solvers/gpu/cufft_poisson_solver.hpp"
#include "solver/fft_solvers/gpu/cufft_split_solver.hpp"
#else
#include "solver/fft_solvers/cpu/fftw_poisson_solver.hpp"
#include "solver/fft_solvers/cpu/fftw_split_solver.hpp"
#endif
