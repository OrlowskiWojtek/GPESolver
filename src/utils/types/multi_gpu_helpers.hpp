#ifndef MULTI_GPU_HELPERS_HPP
#define MULTI_GPU_HELPERS_HPP

#include <cuda_runtime.h>
#include "mat3d/gpu_3d_array.hpp"
#include <iostream>

class GpuDeviceGuard {
public:
    explicit GpuDeviceGuard(int dev_id) {
        CUDA_CHECK(cudaGetDevice(&_previous));
        CUDA_CHECK(cudaSetDevice(dev_id));
    }

    ~GpuDeviceGuard() noexcept {
        cudaSetDevice(_previous);
    }

    GpuDeviceGuard(const GpuDeviceGuard &)           = delete;
    GpuDeviceGuard operator=(const GpuDeviceGuard &) = delete;

private:
    int _previous = {};
};

class GpuStream {
public:
    explicit GpuStream(int device) : device_(device) {
        GpuDeviceGuard guard(device_);
        CUDA_CHECK(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
    }

    ~GpuStream() noexcept {
        if (stream_) {
            GpuDeviceGuard guard(device_);
            cudaStreamDestroy(stream_);
        }
    }

    GpuStream(const GpuStream&) = delete;
    GpuStream& operator=(const GpuStream&) = delete;

    cudaStream_t get() const noexcept { return stream_; }
    int device() const noexcept { return device_; }

private:
    int device_ = {};
    cudaStream_t stream_ = {};
};

template <typename T>
class DeviceArray {
public:
    DeviceArray(int device, std::size_t nx, std::size_t ny, std::size_t nz)
        : device_(device), array_(make_array(device, nx, ny, nz)) {}

    ~DeviceArray() noexcept {
        try {
            GpuDeviceGuard guard(device_);
            array_.reset(); // Free while the owning device is selected.
        } catch (...) {
            // Log in production; destructors must not throw.
        }
    }

    T* data() noexcept { return array_.data(); }
    const T* data() const noexcept { return array_.data(); }
    int device() const noexcept { return device_; }

    void copy_from_host(const T* src) {
        GpuDeviceGuard guard(device_);
        array_.copy_from_host(src);
    }

private:
    static GpuArray<T> make_array(
        int device, std::size_t nx, std::size_t ny, std::size_t nz) {
        GpuDeviceGuard guard(device);
        return GpuArray<T>(nx, ny, nz);
    }

    int device_;
    GpuArray<T> array_;
};

#endif
