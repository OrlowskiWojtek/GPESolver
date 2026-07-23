#pragma once

#include "context/context.hpp"

class CUFFTAbstractGPUSolver {
public:
    CUFFTAbstractGPUSolver(wavefun_gpu_t *psi, pote_gpu_t *fi3d)
        : psi(psi)
        , fi3d(fi3d) {

        };

    virtual ~CUFFTAbstractGPUSolver() {
        psi  = nullptr;
        fi3d = nullptr;
    }

protected:
    wavefun_gpu_t *psi;
    pote_gpu_t *fi3d;
};
