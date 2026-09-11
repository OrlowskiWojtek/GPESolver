#ifndef CPU_IMAGINARY_GPE_SOLVER
#define CPU_IMAGINARY_GPE_SOLVER

#include "solver/solver.hpp"

class CpuITGrossPitaevskiSolver : public AbstractGrossPitaevskiSolver {
public:
    CpuITGrossPitaevskiSolver(AbstractSimulationMediator *mediator);

private:
    CPUSolverData m_data;

    std::unique_ptr<AbstractPoissonSolver> poisson_solver;

    void prepare_fft() override;
    void import_pote() override;
    void import_data() override;
    void export_data() override;

    void calc_energy() override;
    void init_containers() override;
    void calc_fi3d() override;
    void calc_norm() override;
    void normalize() override;
    void iterate() override;
    void adjust(int iter) override;
    void finish() override;
    const int iter_per_summary() const override;

    void imag_iter_linear_step();
    void imag_iter_nonlinear_step();
};

#endif
