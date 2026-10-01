#ifndef CPU_IMAGINARY_GPE_SOLVER
#define CPU_IMAGINARY_GPE_SOLVER

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

    std::unique_ptr<AbstractPoissonSolver> poisson_solver;

    //! arrays for gauss integration u -> ui, w -> wi
    static constexpr std::array<double, 8> ug = {
        0.019855071751231884,
        0.101666761293186630,
        0.237233795041835507,
        0.408282678752175098,
        0.591717321247824902,
        0.762766204958164493,
        0.898333238706813370,
        0.980144928248768116
    };

    static constexpr std::array<double, 8> wg = {
        0.050614268145188129,
        0.111190517226687235,
        0.156853322938943644,
        0.181341891689180991,
        0.181341891689180991,
        0.156853322938943644,
        0.111190517226687235,
        0.050614268145188129
    };

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
