#ifndef GPU_IMAGINARY_GPE_SOLVER
#define GPU_IMAGINARY_GPE_SOLVER

#include "solver/solver.hpp"
#include "solver/solver_data/gpu_solver_data.hpp"

class GpuITGrossPitaevskiSolver : public AbstractGrossPitaevskiSolver {
public:
    GpuITGrossPitaevskiSolver(AbstractSimulationMediator *mediator);
    ~GpuITGrossPitaevskiSolver();

private:
    //! Data used in solver
    GPUSolverData m_data;
    double *d_norm;
    double *d_kin_dev, *d_pot_dev, *d_int_dev, *d_ext_dev, *d_bmf_dev;

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

    void imag_iteration_full();

    void iterate() override;
    void adjust(int iter) override;
    void finish() override;
    const int iter_per_summary() const override;
};

#endif
