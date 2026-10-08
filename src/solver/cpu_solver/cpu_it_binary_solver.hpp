#ifndef CPU_IMAGINARY_BINARY_GPE_SOLVER
#define CPU_IMAGINARY_BINARY_GPE_SOLVER

#include "solver/solver.hpp"
#include "parameters/binary_parameters.hpp"
#include <array>

class CpuITBinaryGrossPitaevskiSolver : public AbstractGrossPitaevskiSolver {
public:
    CpuITBinaryGrossPitaevskiSolver(AbstractSimulationMediator *mediator);

private:
    BinaryMixtureParameters* p_mix;

    //! Data for first component 'a' of a condensate
    CPUSolverData m_data_a;
    //! Data for second component 'b' of a condensate
    CPUSolverData m_data_b;

    //! norm for first component
    double norm_a = 0;
    //! norm for second component
    double norm_b = 0;

    std::unique_ptr<AbstractPoissonSolver> poisson_solver_a;
    std::unique_ptr<AbstractPoissonSolver> poisson_solver_b;

    //! arrays for gauss integration u -> ui, w -> wi
    static std::array<double, 8> ug;
    static std::array<double, 8> wg;

    void lhy_point(double n1, double n2, double &lhy1, double&lhy2);
    void calc_lhy();
    //! LHY potential of component a
    potential_t flhy_a;
    //! LHY potential of component b
    potential_t flhy_b;

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
