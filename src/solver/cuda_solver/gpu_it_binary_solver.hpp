#ifndef GPU_IMAGINARY_BINARY_GPE_SOLVER
#define GPU_IMAGINARY_BINARY_GPE_SOLVER

#include "solver/solver.hpp"
#include "solver/solver_data/gpu_solver_data.hpp"
#include "parameters/binary_parameters.hpp"

class GpuITBinaryGrossPitaevskiSolver : public AbstractGrossPitaevskiSolver {
public:
    GpuITBinaryGrossPitaevskiSolver(AbstractSimulationMediator *mediator);
    ~GpuITBinaryGrossPitaevskiSolver();

private:
    BinaryMixtureParameters* p_mix;

    //! Data for first component 'a' of a condensate
    GPUSolverData m_data_a_0; //!< first gpu
    GPUSolverData m_data_a_1; //!< second gpu
    //! Data for second component 'b' of a condensate
    GPUSolverData m_data_b_0;
    GPUSolverData m_data_b_1;

    double *d_norm_a;
    double *d_norm_b;

    std::unique_ptr<AbstractPoissonSolver> poisson_solver_a;
    std::unique_ptr<AbstractPoissonSolver> poisson_solver_b;

    void lhy_point(double n1, double n2, double &lhy1, double&lhy2);
    void calc_lhy();
    void copy_fi3d_from_gpu1();
    void copy_data_from_gpu1();
    void copy_data_from_gpu0();
    //! LHY potential of component a
    GpuArray<double> flhy_a;
    //! LHY potential of component b
    GpuArray<double> flhy_b;

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

    void imag_iter_full_step();
};

#endif
