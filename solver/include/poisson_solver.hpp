#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "include/params.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <fftw3.h>

/*! class PoissonSolver.
*
* \brief class implements fast poisson solver for mesh.
*
* FFTW has been used in calculations.
*/
class PoissonSolver{
public:
    PoissonSolver();
    ~PoissonSolver();

    void prepare(StdMat3D<std::complex<double>>* psi, StdMat3D<double>* fi3d);
    void execute();
private:

    void prepare_transforms();
    PhysicalParameters* p;
    StdMat3D<std::complex<double>>* psi;
    StdMat3D<double>* fi3d;

    static int FFTW_N_THREADS;

    fftw_plan plan_fwd;
    fftw_plan plan_bwd;
    
    fftw_complex *rho_r;
    fftw_complex *rho_k;
    fftw_complex *Vdip_k;
};

#endif
