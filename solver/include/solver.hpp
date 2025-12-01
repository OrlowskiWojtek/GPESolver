#ifndef GPE_SOLVER_HPP
#define GPE_SOLVER_HPP

#include "include/params.hpp"
#include "include/poisson_solver.hpp"
#include "include/file_manager.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <memory>

/*! Solver of time dependent Gross Pitaevski equation.
 *
 *
 */
class GrossPitaevskiSolver {
public:
    GrossPitaevskiSolver();
    void solve();

private:
    PhysicalParameters *params;

    //! Containers
    //! Wavefunction of bec - copy for calculations.
    StdMat3D<std::complex<double>> cpsii;

    //! Wavefunction of bec - copy for calculations.
    StdMat3D<std::complex<double>> cpsi;

    //! Wavefunction of bec - copy for calculations.
    StdMat3D<std::complex<double>> cpsin;

    //! Map of dipole-dipole potential (from dipole - dipole interaction)
    StdMat3D<double> fi3do;

    //! Map of dipole-dipole potential - copy
    StdMat3D<double> fi3d;

    //! Map of external potential
    StdMat3D<double> pote;

    // current norm of wavefunction
    double xnorma;

    // Energy terms
    double e_kin;   //!< Kinetic energy
    double e_pot;   //!< Potential energy
    double e_int;   //!< Interaction energy
    double e_ext;   //!< Dipole-dipole interaction energy
    double e_bmf;   //!< Beyond mean-field energy
    double e_total; //!< Total energy

    std::unique_ptr<PoissonSolver> poisson_solver;
    std::unique_ptr<FileManager> file_manager;

    void calc_initial_state();
    void calc_evolution();
    void calc_energy();

    void init_containers();
    void init_with_cos();
    void init_with_gauss();
    void init_potential();
    void free_potential_well();
    void imag_time_iter();
    void real_time_iter();

    void save_xy_cut_to_file(int iter);
    void save_ene_to_file(int iter);
    void save_x_cut_to_file();
    void load_initial_state_from_file();
    void save_initial_state_to_file();

    void calc_fi3d();
    void calc_norm();
    void normalize();
    void imag_iter_linear_step();
    void imag_iter_nonlinear_step();
    void real_iter_linear_step();
    void real_iter_nonlinear_step();

    double pote_value(int ix, int iy, int iz);
    double pote_released_value(int ix, int iy, int iz);
};

#endif
