#pragma once

#include "context/context.hpp"

class FFTWAbstractCPUSolver {
public:
    virtual ~FFTWAbstractCPUSolver() {
        psi = nullptr;
        fi3d = nullptr;
        pote = nullptr;
    };
    FFTWAbstractCPUSolver(wavefunction_t *psi, potential_t *fi3d, potential_t *pote)
        : psi(psi)
        , fi3d(fi3d)
        , pote(pote) {
    }

protected:
    wavefunction_t *psi;
    potential_t *fi3d;
    potential_t *pote;
};
