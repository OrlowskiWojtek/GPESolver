#pragma once

#include <cstddef>
#include <utility>

template <class T>
class IMat3D{
public:
    IMat3D() = default;

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
    IMat3D(const IMat3D &)            = delete;
    IMat3D &operator=(const IMat3D &) = delete;

    // Data 
    virtual T *data() = 0;
    virtual const T *data() const;

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
    void swap(IMat3D &other) noexcept {
        std::swap(nx_, other.nx_);
        std::swap(ny_, other.ny_);
        std::swap(nz_, other.nz_);

        swap_data();
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
        if (size > 0) {
            CUDA_CHECK(cudaMalloc(&data_, size * sizeof(T)));
            CUDA_CHECK(cudaMemset(data_, 0, size * sizeof(T)));
        } else {
            throw std::runtime_error("Wrong size of alloc");
        }
    }

    T *data_        = nullptr;
    std::size_t nx_ = 0, ny_ = 0, nz_ = 0;
};
