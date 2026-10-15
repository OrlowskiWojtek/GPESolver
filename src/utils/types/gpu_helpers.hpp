#ifndef GPU_HELPERS_HPP

#include <cstdio>
#include <cuda_runtime.h>
#include <stdexcept>

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t _err = call;                                               \
        if (_err != cudaSuccess) {                                             \
            fprintf(stderr,                                                    \
                    "CUDA error at %s:%d: %s\n",                               \
                    __FILE__,                                                  \
                    __LINE__,                                                  \
                    cudaGetErrorString(_err));                                 \
            throw std::runtime_error(cudaGetErrorString(_err));                \
        }                                                                      \
    } while (0)

// Energies of GPE stored on GPU
class GpuEnergy {
public:
    GpuEnergy() {
        CUDA_CHECK(cudaMalloc(&d_enes_dev, n_enes * sizeof(double)));
    }

    ~GpuEnergy() noexcept {
        cudaFree(d_enes_dev);
    }

    double *d_kin() noexcept {
        return d_enes_dev + 0;
    }
    double *d_pot() noexcept {
        return d_enes_dev + 1;
    }
    double *d_contact() noexcept {
        return d_enes_dev + 2;
    }
    double *d_dipol_dipol() noexcept {
        return d_enes_dev + 3;
    }
    double *d_lhy() noexcept {
        return d_enes_dev + 4;
    }

    double *get_data() noexcept{
        return d_enes_dev;
    }

private:
    //! Values to be stored on gpu
    //! 0 -> kinetic
    //! 1 -> potential
    //! 2 -> contact
    //! 3 -> dipol-dipol
    //! 4 -> LHY
    double *d_enes_dev  = nullptr;
    const size_t n_enes = 5;
};

#endif
