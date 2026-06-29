#ifndef GPU_GPE_SOLVER
#define GPU_GPE_SOLVER

#include "solver/solver_data/gpu_solver_data.hpp"
#include "solver/solver.hpp"

class GpuGrossPitaevskiSolver: public AbstractGrossPitaevskiSolver {
public:
    GpuGrossPitaevskiSolver(AbstractSimulationMediator* mediator);

private:
    //! Data used in solver
    GPUSolverData m_data;

    void prepare_fft() override;
    void import_pote() override;
    void import_data() override;
    void export_data() override;

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
