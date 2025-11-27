#ifndef STD_MAT_3D_HPP
#define STD_MAT_3D_HPP

#include <vector>

template <class T>
class StdMat3D {
public:
    // Constructors
    StdMat3D();
    StdMat3D(int nx, int ny, int nz);

    // Resize
    void resize(int nx, int ny, int nz);

    // at operator
    inline T &operator()(int nx, int ny, int nz) noexcept;

    inline int get_index(int i, int j, int k) noexcept {
        return (i * _ny + j) * _nz + k;
    }

    inline T& operator()(int idx) noexcept;

private:
    int _nx, _ny, _nz;
    int N;
    std::vector<T> data;
};

template <class T>
T& StdMat3D<T>::operator()(int idx) noexcept {
    return data[idx];
}

template <class T>
StdMat3D<T>::StdMat3D() {
    _nx = 0;
    _ny = 0;
    _nz = 0;

    N = 0;
}

template <class T>
StdMat3D<T>::StdMat3D(int nx, int ny, int nz) {
    _nx = nx;
    _ny = ny;
    _nz = nz;
    N = nx * ny * nz;

    data.resize(nx * ny * nz, 0);
}

template <class T>
void StdMat3D<T>::resize(int nx, int ny, int nz) {
    _nx = nx;
    _ny = ny;
    _nz = nz;
    N = nx * ny * nz;

    data.resize(nx * ny * nz, 0);
}

template <class T>
inline T &StdMat3D<T>::operator()(int i, int j, int k) noexcept {
    return data[get_index(i, j, k)];
}

#endif
