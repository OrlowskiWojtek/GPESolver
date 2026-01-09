#ifndef FFT_CONTEXT_HPP
#define FFT_CONTEXT_HPP

#include "parameters/parameters.hpp"
#include "context/context.hpp"
#include <fftw3.h>

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

    fftw_plan plan_fwd;
    fftw_plan plan_bwd;
    
    fftw_complex *rho_r;
    fftw_complex *rho_k;
};

#endif
