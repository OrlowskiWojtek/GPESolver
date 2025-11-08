#ifndef MAT_3D_HPP
#define MAT_3D_HPP

#include <cstddef>
#include <utility>

template <class T>
class Mat3D {
public:
    // Constructors
    Mat3D();
    Mat3D(size_t nx, size_t ny, size_t nz);

    // Destructor
    ~Mat3D();

    // Copy constructor
    Mat3D(const Mat3D &);

    // Copy assignment / move operator
    Mat3D &operator=(Mat3D);

    // Move constructor
    Mat3D(Mat3D &&);

    // Resize
    void resize(size_t nx, size_t ny, size_t nz);

    // at operator
    T &operator()(size_t nx, size_t ny, size_t nz);

private:
    size_t _nx, _ny, _nz;
    T ***data;
};

template <class T>
Mat3D<T>::Mat3D() {
    _nx = 0;
    _ny = 0;
    _nz = 0;
    data = nullptr;
}

template <class T>
Mat3D<T>::Mat3D(size_t nx, size_t ny, size_t nz) {
    _nx = nx;
    _ny = ny;
    _nz = nz;

    data = new T **[nx];
    for (size_t i = 0; i < nx; i++) {
        data[i] = new T *[ny];
        for (size_t j = 0; j < ny; j++) {
            data[i][j] = new T[nz];
        }
    }
}

template <class T>
Mat3D<T>::~Mat3D() {
    if(!data) return;

    for (size_t i = 0; i < _nx; i++) {
        for (size_t j = 0; j < _ny; j++) {
            delete[] data[i][j];
        }
        delete[] data[i];
    }
    delete[] data;

    _nx = 0;
    _ny = 0;
    _nz = 0;
}

template <class T>
Mat3D<T>::Mat3D(const Mat3D &_other) {
    _nx = _other._nx;
    _ny = _other._ny;
    _nz = _other._nz;

    data = new T **[_nx];
    for (size_t i = 0; i < _nx; i++) {
        data[i] = new T *[_ny];
        for (size_t j = 0; j < _ny; j++) {
            data[i][j] = new T[_nz];
        }
    }

    for (size_t i = 0; i < _nx; i++) {
        for (size_t j = 0; j < _ny; j++) {
            for (size_t k = 0; k < _nz; k++) {
                data[i][j][k] = _other.data[i][j][k];
            }
        }
    }
}

template <class T>
Mat3D<T>::Mat3D(Mat3D &&_other) {
    _nx  = _other._nx;
    _ny  = _other._ny;
    _nz  = _other._nz;
    data = _other.data;

    _other._nx = 0;
    _other._ny = 0;
    _other._nz = 0;
    _other.data = nullptr;
}

template <class T>
Mat3D<T> &Mat3D<T>::operator=(Mat3D _other) {
    std::swap(_nx, _other._nx);
    std::swap(_ny, _other._ny);
    std::swap(_nz, _other._nz);
    std::swap(data, _other.data);

    return *this;
}

template <class T>
void Mat3D<T>::resize(size_t nx, size_t ny, size_t nz){
    *this = Mat3D(nx, ny, nz);
}

template <class T>
T &Mat3D<T>::operator()(size_t i, size_t j, size_t k) {
    return data[i][j][k];
}

#endif
