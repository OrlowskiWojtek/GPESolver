#include "include/params.hpp"
#include "units.hpp"
#include <cmath>

void PhysicalParameters::init_values() {
    wzl    = 120 * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - z direction
    wrl    = 60. * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - y direction

    double dd = 1500 / 0.05292; // TODO: do not know what it is - used in aa
    
    aa = m * std::pow(wrl, 2) / 4. / std::pow(dd,2);
    b =  0.5 * m * std::pow(wrl, 2);

    y0     = 0;
    omega0 = 0;
    m      = 0;
}

void PhysicalParameters::init_containers() {
       
}

void PhysicalParameters::init_rdy() {
    rdy.resize(4 * nx, 4 * ny, 4 * nz);
}
