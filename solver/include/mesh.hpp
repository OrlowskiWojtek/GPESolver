#ifndef SOLVER_MESH_HPP
#define SOLVER_MESH_HPP

#include <cstddef>
#include <functional>
#include "mat3d/stdmat3d.hpp"

// Ideas:
// Sizes are static, to make sure that all meshes has same size
// Meshes can be divided into separates types because few copies are needed usually
// Can also have params like dx dy and dz
// Easier to parallyze (can divide iteration inside this class with some copied operation)
template <class T>
class AbstractMesh{
public:
    inline T operator()(int i, int j, int k);

    void iterate(std::function<void(double)>&&);
    void iterate_no_boundary(std::function<void(double)>);

protected:
    StdMat3D<T> data;

    static size_t nx;
    static size_t ny;
    static size_t nz;

    virtual void fill() = 0;
};

#endif
