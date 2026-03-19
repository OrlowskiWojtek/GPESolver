#ifndef FFT_CONTEXT_HPP
#define FFT_CONTEXT_HPP

#include "parameters/parameters.hpp"
#include "context/context.hpp"
#ifdef USE_CUDA
#include <cufft.h>
#include <cuda_runtime.h>
#else
#include <fftw3.h>
#endif

#ifdef USE_CUDA
using plan_type = cufftHandle;
using complex_type = cufftDoubleComplex;
#else
using plan_type = fftw_plan;
using complex_type = fftw_complex;
#endif

inline double& real(complex_type& num){
#ifdef USE_CUDA
    return num.x;
#else
    return num[0];
#endif
}

inline double& imag(complex_type& num){
#ifdef USE_CUDA
    return num.y;
#else
    return num[1];
#endif
}

/*! class FFT_CONTEXT.
*
* \brief Abstract interface for FFTW calculations.
*/
class FFTContext{
public:
    FFTContext();
    virtual ~FFTContext();

    void prepare(wavefunction_t* psi, potential_t* fi3d);
    virtual void execute() = 0;

private:
    virtual void prepare_transforms() = 0;
    virtual void prepare_containers() = 0;

protected:
    wavefunction_t* psi;
    potential_t* fi3d;
    PhysicalParameters* p;

    static int FFTW_N_THREADS;

    plan_type plan_fwd;
    plan_type plan_bwd; 
};

#endif
