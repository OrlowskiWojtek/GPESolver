#pragma once

#include "context/context.hpp"

class FFTWAbstractCPUSolver {
public:
    virtual ~FFTWAbstractCPUSolver() {
        psi = nullptr;
        fi3d = nullptr;
    };
    FFTWAbstractCPUSolver(wavefunction_t *psi, potential_t *fi3d)
        : psi(psi)
        , fi3d(fi3d)
        {
    }

protected:
    wavefunction_t *psi;
    potential_t *fi3d;
};
