#ifndef CPU_GPE_SOLVER
#define CPU_GPE_SOLVER

#include "solver/solver.hpp"

class CpuGrossPitaevskiSolver: public AbstractGrossPitaevskiSolver {
public:
    CpuGrossPitaevskiSolver(AbstractSimulationMediator* mediator);

private:
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
