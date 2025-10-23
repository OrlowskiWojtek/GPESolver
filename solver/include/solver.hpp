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
    PhysicalParameters params;

    //! Containers
    //! Wavefunction of single particle - copy for calculations.
    CMat3D cpsii;

    //! Wavefunction of single particle - copy for calculations.
    CMat3D cpsi;

    //! Wavefunction of single particle - copy for calculations.
    CMat3D cpsin;

    //! Map of potential (from dipole - dipole interaction)
    Mat3D fi3do;

    //! Map of potential - copy
    Mat3D fi3d;

};

#endif
