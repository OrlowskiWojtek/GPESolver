#ifndef STD_MAT_3D_HPP
#define STD_MAT_3D_HPP

#include <cstddef>
#include <vector>

template <class T>
class StdMat3D {
public:
    // Constructors
    StdMat3D();
    StdMat3D(size_t nx, size_t ny, size_t nz);

    // Resize
    void resize(size_t nx, size_t ny, size_t nz);

    // at operator
    T &operator()(size_t nx, size_t ny, size_t nz);

private:
    size_t _nx, _ny, _nz;
    std::vector<std::vector<std::vector<T>>> data;
};

template <class T>
StdMat3D<T>::StdMat3D() {
    _nx = 0;
    _ny = 0;
    _nz = 0;
}

template <class T>
StdMat3D<T>::StdMat3D(size_t nx, size_t ny, size_t nz) {
    _nx = nx;
    _ny = ny;
    _nz = nz;

    data.resize(nx);
    for (size_t i = 0; i < nx; i++) {
        data[i].resize(ny);
        for (size_t j = 0; j < ny; j++) {
            data[i][j].resize(nz);
        }
    }
}

template <class T>
void StdMat3D<T>::resize(size_t nx, size_t ny, size_t nz) {
    _nx = nx;
    _ny = ny;
    _nz = nz;

    data.resize(nx);
    for (size_t i = 0; i < nx; i++) {
        data[i].resize(ny);
        for (size_t j = 0; j < ny; j++) {
            data[i][j].resize(nz);
        }
    }
}

template <class T>
T &StdMat3D<T>::operator()(size_t i, size_t j, size_t k) {
    return data[i][j][k];
}

#endif
