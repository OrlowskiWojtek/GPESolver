#include "initializer/initializer.hpp"
#include "output.hpp"
#include "parameters/potentials.hpp"
#include "parameters/wavefunctions.hpp"
#include <random>

DataInitializer::DataInitializer(AbstractSimulationMediator *_mediator)
    : params(PhysicalParameters::getInstance())
    , p_mediator(_mediator)
    , p_sctx(SimulationContext::getInstance()) {
}

void DataInitializer::initialize_wavefunction(
    const WavefunctionInitConfig &conf) {
    OutputFormatter::printInfo("Initializing wavefunction");

    if (conf.key == "BINARY_FILE") {
        init_from_binary_file();
    } else if (conf.key == "TEXT_FILE") {
        init_from_text_file();
    } else {
        init_wavefunction(conf);
    }

    if (conf.add_random_noise) {
        add_initial_noise();
    }

    p_mediator->on_data_initialized(_data);
}

void DataInitializer::initialize_potential(const PotentialInitConfig &conf) {
    OutputFormatter::printInfo("Initializing potential");

    init_pote(conf);

    if (conf.add_random_noise) {
        OutputFormatter::printInfo("Adding random noise to potential");
        add_initial_pote_noise();
    }

    p_mediator->on_pote_initialized(_pote);
}

void DataInitializer::init_from_text_file() {
    p_mediator->request_load_from_text(_data);
}

void DataInitializer::init_from_binary_file() {
    p_mediator->request_load_from_binary(_data);
}

void DataInitializer::init_pote(const PotentialInitConfig &conf) {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    _pote_func = PotentialRegistry::instance().get_function(conf.key);
    _pote.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double x = p_sctx->get_x(i);
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                _pote(i, j, k) = _pote_func(x, y, z);
            }
        }
    }
}

void DataInitializer::init_wavefunction(const WavefunctionInitConfig &conf) {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    _data_func = InitializerRegistry::instance().get_function(conf.key);
    _data.resize(nx, ny, nz);

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double x = p_sctx->get_x(i);
                double y = p_sctx->get_y(j);
                double z = p_sctx->get_z(k);

                _data(i, j, k) = _data_func(x, y, z);
            }
        }
    }
}

void DataInitializer::add_initial_noise() {
    const int nx = params->nx;
    const int ny = params->ny;
    const int nz = params->nz;

    // default seed explicitly added
    const unsigned int seed  = 5489;
    const double noise_level = 0.05;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-noise_level, noise_level);

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double rand_real = dist(rng);

                _data(i, j, k) *= (1 + rand_real);
            }
        }
    }
}

void DataInitializer::add_initial_pote_noise() {
    const int nx = params->nx;
    const int ny = params->ny;
    const int nz = params->nz;

    // default seed explicitly added
    const unsigned int seed  = 5489u;
    const double noise_level = 0.05;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-noise_level, noise_level);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double rand_real = dist(rng);

                _pote(i, j, k) *= (1. + rand_real);
            }
        }
    }
}
