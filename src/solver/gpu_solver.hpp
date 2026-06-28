#ifndef GPU_GPE_SOLVER
#define GPU_GPE_SOLVER

#include "solver/solver.hpp"
#include "mat3d/gpu_3d_array.hpp"
#include "cuComplex.h"

class GpuGrossPitaevskiSolver: public AbstractGrossPitaevskiSolver {
public:
    GpuGrossPitaevskiSolver(AbstractSimulationMediator* mediator);

private:
    //! Containers
    //! Wavefunction of bec + copy for calculations.
    GpuArray<cuDoubleComplex> cpsii_gpu;
    GpuArray<cuDoubleComplex> cpsi_gpu;

    //! Map of dipole-dipole potential
    GpuArray<double> fi3d_gpu;

    //! Map of external potential
    GpuArray<double> pote_gpu;

    // TODO: for now need to override those
    void load_buffer(const wavefunction_t&) override;
    void load_pote(const potential_t&) override;
    void initialize() override;

    void load_psi() override;

    void calc_energy() override;
    void init_containers() override;
    void calc_fi3d() override;
    void calc_norm() override;
    void normalize() override;
    void imag_iter_linear_step() override;
    void imag_iter_nonlinear_step() override;
    void real_fft_potential_half_step() override;
    void real_fft_kinetic_step() override;
};

#endif
