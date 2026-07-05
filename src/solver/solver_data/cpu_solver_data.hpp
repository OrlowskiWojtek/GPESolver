#pragma once

#include "solver/solver_data/isolver_data.hpp"
#include "context/context.hpp"

/*! CPU-side solver data container.
 * 
 * stores data used in calculations in CPU format (wavefunction_t, potential_t).
 * Inherits ISolverData interface.
 */
class CPUSolverData : public ISolverData {
public:
    CPUSolverData() = default;

    void allocate(size_t nx, size_t ny, size_t nz) override {
        cpsi.resize(nx, ny, nz);
        cpsii.resize(nx, ny, nz);
        fi3d.resize(nx, ny, nz);
        pote.resize(nx, ny, nz);
    }

    //! Wavefunction of bec.
    wavefunction_t cpsii;
    //! Wavefunction of bec - copy for nonlinear calculations
    wavefunction_t cpsi;

    //! Map of dipole-dipole potential
    potential_t fi3d;

    //! Map of external potential
    potential_t pote;
};
