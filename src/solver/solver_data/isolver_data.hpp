#pragma once

#include <cstddef>

/*! Abstract interface for solver data.
 *
 * Separates the data layer from the computation layer.
 */
class ISolverData {
public:
    virtual ~ISolverData() = default;

    //! Allocate memory for all containers
    virtual void allocate(size_t nx, size_t ny, size_t nz) = 0;
};
