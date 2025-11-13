#ifndef SOLVER_MESH_HPP
#define SOLVER_MESH_HPP

#include <cstddef>
#include "mat3d/stdmat3d.hpp"

template <class T>
class AbstractMesh{
protected:
    StdMat3D<T> data;

    static size_t nx;
    static size_t ny;
    static size_t nz;

    virtual void fill() = 0;
};

#endif
