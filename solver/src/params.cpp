#include "include/params.hpp"
#include "units.hpp"
#include <cmath>

void PhysicalParameters::init_values() {
    wzl    = 120 * 4.1356e-12 / 27211.6;
    wrl    = 60. * 4.1356e-12 / 27211.6;

    y0     = 0;
    omega0 = 0;
    m      = 0;
}

void PhysicalParameters::init_containers() {
}

void PhysicalParameters::init_rdy() {
    rdy.resize(4 * nx, 4 * ny, 4 * nz);
}
