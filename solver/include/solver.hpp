#ifndef GPE_SOLVER_HPP
#define GPE_SOLVER_HPP

#include "include/params.hpp"

//! \todo Just meanwhile sollution, Rather will use blaze (CUDA parallelization, 3D arrays)
#include <armadillo>
using CMat3D    = arma::Cube<std::complex<double>>;
using Mat3D     = arma::Cube<double>;

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
    CMat3D cpsii;

    //! Wavefunction of bec - copy for calculations.
    CMat3D cpsi;

    //! Wavefunction of bec - copy for calculations.
    CMat3D cpsin;

    //! Map of dipole-dipole potential (from dipole - dipole interaction)
    Mat3D fi3do;

    //! Map of dipole-dipole potential - copy
    Mat3D fi3d;

    //! Map of external potential
    Mat3D pote;

    void calc_initial_state();

    void init_containers();
    void init_with_cos();
    void init_potential();
    void free_potential_well();
    void imag_time_iter();

    void calc_fi3d();

    double pote_value(int ix, int iy, int iz);
    double pote_released_value(int ix, int iy, int iz);
};

#endif
