#include "solver/fft_context.hpp"
#include "parameters/parameters.hpp"
#ifdef USE_CUDA
#include <cufft.h>
#include <cuda_runtime.h>
#else
#include <fftw3.h>
#endif

int FFTContext::FFTW_N_THREADS = 4;

FFTContext::FFTContext()
    : p(PhysicalParameters::getInstance()) {
}

FFTContext::~FFTContext() {
#ifdef USE_CUDA
    if (plan_fwd) cufftDestroy(plan_fwd);
    if (plan_bwd) cufftDestroy(plan_bwd);
#else
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);
#endif
}

void FFTContext::prepare(wavefunction_t *cpsi, potential_t *_fi3d, potential_t *_pote) {
    psi  = cpsi;
    fi3d = _fi3d;
    pote = _pote;

#ifndef USE_CUDA
    FFTW_N_THREADS = p->fftw_n_threads;
    int res = fftw_init_threads();

    if (res == 0) {
        throw std::runtime_error("FFTW thread initialization failed!");
    }
#endif

    prepare_transforms();
    prepare_containers();
}
