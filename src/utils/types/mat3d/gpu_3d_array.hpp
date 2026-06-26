#ifndef GPU_3D_ARRAY_HPP
#define GPU_3D_ARRAY_HPP

#include <cstddef>
#include <stdexcept>
#include <utility> // std::exchange
#include <vector>
#include <iostream>

#include "cuda_runtime.h"

// Helper macro lub funkcja do sprawdzania błędów
#define CUDA_CHECK(call)                                                                           \
    do {                                                                                           \
        cudaError_t _err = call;                                                                   \
        if (_err != cudaSuccess) {                                                                 \
            fprintf(stderr,                                                                        \
                    "CUDA error at %s:%d: %s\n",                                                   \
                    __FILE__,                                                                      \
                    __LINE__,                                                                      \
                    cudaGetErrorString(_err));                                                     \
            throw std::runtime_error(cudaGetErrorString(_err));                                    \
        }                                                                                          \
    } while (0)

template <typename T>
class GpuArray {
public:
    GpuArray() = default;

    explicit GpuArray(std::size_t size)
        : nx_(size)
        , ny_(1)
        , nz_(1) {
        allocate(size);
    }

    GpuArray(std::size_t nx, std::size_t ny, std::size_t nz)
        : nx_(nx)
        , ny_(ny)
        , nz_(nz) {
        allocate(nx * ny * nz);
    }

    ~GpuArray() {
        if (data_) {
            cudaFree(data_);
        }
    }

    // Move constructor
    GpuArray(GpuArray &&other) noexcept
        : data_(std::exchange(other.data_, nullptr))
        , nx_(std::exchange(other.nx_, 0))
        , ny_(std::exchange(other.ny_, 0))
        , nz_(std::exchange(other.nz_, 0)) {
    }

    // Move assignment
    GpuArray &operator=(GpuArray &&other) {
        if (this != &other) {
            if (data_) {
                CUDA_CHECK(cudaFree(data_));
            }
            data_ = std::exchange(other.data_, nullptr);
            nx_   = std::exchange(other.nx_, 0);
            ny_   = std::exchange(other.ny_, 0);
            nz_   = std::exchange(other.nz_, 0);
        }
        return *this;
    }

    // No copy (expensive!)
    GpuArray(const GpuArray &)            = delete;
    GpuArray &operator=(const GpuArray &) = delete;

    // Dostęp
    T *data() {
        return data_;
    }
    const T *data() const {
        return data_;
    }

    std::size_t size() const {
        return nx_ * ny_ * nz_;
    }
    std::size_t nx() const {
        return nx_;
    }
    std::size_t ny() const {
        return ny_;
    }
    std::size_t nz() const {
        return nz_;
    }

    bool empty() const {
        return data_ == nullptr;
    }

    // Swap - przydatne dla niektórych kontenerów
    void swap(GpuArray &other) noexcept {
        std::swap(data_, other.data_);
        std::swap(nx_, other.nx_);
        std::swap(ny_, other.ny_);
        std::swap(nz_, other.nz_);
    }

    // Reset - zwolnij pamięć
    void reset() {
        if (data_) {
            CUDA_CHECK(cudaFree(data_));
        }
        data_ = nullptr;
        nx_ = ny_ = nz_ = 0;
    }

    // Indeksacja 3D -> 1D (row-major)
    std::size_t idx(std::size_t i, std::size_t j, std::size_t k) const {
        return (i * ny_ + j) * nz_ + k;
    }

    // Synchroniczne kopiowanie
    void copy_from_host(const T *host_data) {
        CUDA_CHECK(cudaMemcpy(data_, host_data, size() * sizeof(T), cudaMemcpyHostToDevice));
    }

    void copy_from_host(const std::vector<T> &host_data) {
        copy_from_host(host_data.data());
    }

    void copy_to_host(T *host_data) const {
        CUDA_CHECK(cudaMemcpy(host_data, data_, size() * sizeof(T), cudaMemcpyDeviceToHost));
    }

    std::vector<T> copy_to_host() const {
        std::vector<T> result(size());
        copy_to_host(result.data());
        return result;
    }

private:
    void allocate(std::size_t size) {
        std::cout << "Alocating memory" << std::endl;
        if (size > 0) {
            CUDA_CHECK(cudaMalloc(&data_, size * sizeof(T)));
        } else {
            throw std::runtime_error("Wrong size of alloc");
        }
    }

    T *data_        = nullptr;
    std::size_t nx_ = 0, ny_ = 0, nz_ = 0;
};

#endif
