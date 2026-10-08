#ifndef CPU_IMAGINARY_BINARY_GPE_SOLVER
#define CPU_IMAGINARY_BINARY_GPE_SOLVER

#include "solver/solver.hpp"
#include <array>

struct BinaryMixtureParameters {
    BinaryMixtureParameters(const BinaryMixtureParameters &)            = delete;
    BinaryMixtureParameters &operator=(const BinaryMixtureParameters &) = delete;

    static BinaryMixtureParameters *getInstance() {
        if (!instance) {
            instance = new BinaryMixtureParameters();
        }

        return instance;
    }

    //! Mass of atom;
    double m_a;
    double m_b;

    // number of atoms
    double n_atoms_a;
    double n_atoms_b;

    double g11, g22, g12;       // contact interaction
    // cdd12 = \mu_0 \mu_a \mu_b
    double cdd11, cdd22, cdd12;   // dipole strengths

    //! Adds random noise at the beginning of calculations
    bool add_random_noise = false;

    void set_to_default();
private:

    BinaryMixtureParameters() {};

    static BinaryMixtureParameters *instance;
};

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
