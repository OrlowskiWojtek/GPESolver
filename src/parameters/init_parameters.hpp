#ifndef INITIALIZING_PARAMETERS_HPP
#define INITIALIZING_PARAMETERS_HPP

#include <string>

struct WavefunctionInitConfig {
    std::string key;
    bool add_random_noise = false;
    std::string load_filename;
    double dd = 0.0;
    unsigned int random_seed = 5489u;

    int n_gauss_max = 0;
    int bec_droplets_x = 0;
    int bec_droplets_y = 0;
    int bec_droplets_z = 0;
};

struct PotentialInitConfig {
    std::string key;
    bool add_random_noise = false;
};

struct InitializingParameters {
    WavefunctionInitConfig wavefunction;
    PotentialInitConfig potential;
};

#endif
