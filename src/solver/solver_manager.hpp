#ifndef SOLVER_MANAGER_HPP
#define SOLVER_MANAGER_HPP

#include "context/context.hpp"

class AbstractGPESolver {
    // those thing need to be in solver, not solver managar
    virtual void calc_energy()                  = 0;
    virtual void calc_fi3d()                    = 0;
    virtual void calc_norm()                    = 0;
    virtual void normalize()                    = 0;
    virtual void imag_iter_linear_step()        = 0;
    virtual void imag_iter_nonlinear_step()     = 0;
    virtual void real_fft_potential_half_step() = 0;
    virtual void real_fft_kinetic_step()        = 0;
};

class SolverManager {
public:
    void solve();

    void initialize();
    void load_buffer(const wavefunction_t &);
    void load_pote(const potential_t &);

private:
    AbstractGPESolver* solver;

    void calc_initial_state();
    void calc_evolution();
    void calc_cradle();

    void init_containers();
    void free_potential_well();

    void imag_time_iter();
    void real_time_iter();
    void summarize_imag_iter();
    void summarize_real_iter();
};

#endif
