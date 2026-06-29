#pragma once

#include "solver/solver_data/isolver_data.hpp"
#include "context/context.hpp"

/*! GPU-side solver data.
 * 
 * Public GPU arrays.
 */
class GPUSolverData : public ISolverData {
public:
    void allocate(size_t nx, size_t ny, size_t nz) override {
        cpsi_gpu  = GpuArray<cuDoubleComplex>(nx, ny, nz);
        cpsii_gpu = GpuArray<cuDoubleComplex>(nx, ny, nz);
        fi3d_gpu  = GpuArray<double>(nx, ny, nz);
        pote_gpu  = GpuArray<double>(nx, ny, nz);
    }

    GpuArray<cuDoubleComplex> cpsi_gpu;
    GpuArray<cuDoubleComplex> cpsii_gpu;
    GpuArray<double> fi3d_gpu;
    GpuArray<double> pote_gpu;
};
