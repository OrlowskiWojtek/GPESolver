#include "context/context.hpp"

SimulationContext* SimulationContext::instance = nullptr;

SimulationContext* SimulationContext::getInstance(){
    if(instance == nullptr){
        instance = new SimulationContext;
    }

    return instance;
}

SimulationContext::SimulationContext()
    : params(PhysicalParameters::getInstance()) {
}

void SimulationContext::initialize() {
    x_vec.resize(params->nx);
    y_vec.resize(params->ny);
    z_vec.resize(params->nz);

    for (int i = 0; i < params->nx; i++) {
        double x = (i - (static_cast<int>(params->nx / 2.) + 1)) * params->dx;
        x_vec[i] = x;
    }
    for (int j = 0; j < params->ny; j++) {
        double y = (j - (static_cast<int>(params->ny / 2.) + 1)) * params->dy;
        y_vec[j] = y;
    }
    for (int k = 0; k < params->nz; k++) {
        double z = (k - (static_cast<int>(params->nz / 2.) + 1)) * params->dz;
        z_vec[k] = z;
    }
}
