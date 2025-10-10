#include "include/params.hpp"
#include "units.hpp"

void PhysicalParameters::init_values() {
    wzl    = 0;
    wrl    = 0;
    y0     = 0;
    omega0 = 0;
    m      = 0;
}

void PhysicalParameters::init_containers() {
    rdy.resize(0, 0, 0);
}
