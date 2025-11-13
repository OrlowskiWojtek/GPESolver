#ifndef GPE_SOLVER_HPP
#define GPE_SOLVER_HPP

#include "include/params.hpp"
#include "include/poisson_solver.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <memory>

/*! Solver of time dependent Gross Pitaevski equation.
*
*
*/
class GrossPitaevskiSolver{
public:
    GrossPitaevskiSolver();
    void solve();

private:
    PhysicalParameters* params;

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

    double xnorma;

    std::unique_ptr<PoissonSolver> poisson_solver;

    void calc_initial_state();
    void calc_evolution();

    void init_containers();
    void init_with_cos();
    void init_potential();
    void free_potential_well();
    void imag_time_iter();
    void real_time_iter();

    // \todo: temporary - remove iter
    void save_xy_cut_to_file(int iter);
    void save_x_cut_to_file();

    void calc_fi3d();
    void calc_norm();

    double pote_value(int ix, int iy, int iz);
    double pote_released_value(int ix, int iy, int iz);
};

#endif
