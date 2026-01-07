#include "solver/fft_context.hpp"
#include "parameters/parameters.hpp"
#include <fftw3.h>

int FFTContext::FFTW_N_THREADS = 4;

FFTContext::FFTContext()
    : p(PhysicalParameters::getInstance()) {
}

FFTContext::~FFTContext() {
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);
}

void FFTContext::prepare(StdMat3D<std::complex<double>> *cpsi, StdMat3D<double> *_fi3d) {
    psi  = cpsi;
    fi3d = _fi3d;

    FFTW_N_THREADS = p->fftw_n_threads;
    int res = fftw_init_threads();
    if (res == 0) {
        throw std::runtime_error("FFTW thread initialization failed!");
    }

    prepare_transforms();
    prepare_containers();
}
