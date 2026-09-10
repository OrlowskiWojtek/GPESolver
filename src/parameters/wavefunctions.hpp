#ifndef WAVEFUNCTIONS_HPP
#define WAVEFUNCTIONS_HPP
/*! File containing initial wavefunctions for solver */

#include "context/context.hpp"
#include "parameters/parameters.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class InitializerRegistry {
public:
    using WavRegisterKey = std::string;
    using WavID          = unsigned int;
    using init_func_t    = std::function<std::complex<double>(double, double, double)>;

    struct InitializerInfo {
        WavRegisterKey key;
        WavID id;
        init_func_t func;
    };

    static InitializerRegistry &instance() {
        static InitializerRegistry registry;
        return registry;
    }

    // Register a potential with key, id, and function
    void register_potential(WavRegisterKey key, WavID id, init_func_t func) {
        initializers_[key] = {key, id, func};
        id_map_[id]        = key;
    }

    void print_options() {
        std::cout << "INITIALIZATION OPTIONS ARE: \n";
        std::for_each(initializers_.begin(), initializers_.end(), [](const auto &info) {
            std::cout << " " << info.second.key << " ";
        });
        std::cout << std::endl;
    }

    // Getters
    init_func_t get_function(const WavRegisterKey &key) const {
        auto it = initializers_.find(key);

        if (it == initializers_.end()) {
            std::for_each(initializers_.begin(), initializers_.end(), [](const auto &info) {
                std::cout << " " << info.second.key << " ";
            });
            std::cout << std::endl;

            std::cout << "Wrong initial wavefunction, options are listed above, you used: " << key; 

            throw std::runtime_error("Wrong init wavefunction");
        }

        return (it != initializers_.end()) ? it->second.func
                                           : initializers_.find("GAUSS")->second.func;
    }

    init_func_t get_function(WavID id) const {
        auto it = id_map_.find(id);
        return (it != id_map_.end()) ? get_function(it->second)
                                     : initializers_.find("GAUSS")->second.func;
    }

    bool contains(const WavRegisterKey &key) {
        auto it = initializers_.find(key);

        return it != initializers_.end();
    }

private:
    std::unordered_map<WavRegisterKey, InitializerInfo> initializers_;
    std::unordered_map<WavID, WavRegisterKey> id_map_;
};

static unsigned int WavID = 0;

#define REGISTER_INITIALIZER(KEY, FUNC)                                                            \
    static auto _registrar_##KEY = []() {                                                          \
        InitializerRegistry::instance().register_potential(#KEY, WavID++, FUNC);                   \
        return 0.;                                                                                 \
    }();

REGISTER_INITIALIZER(COS, [](double x, double y, double z) -> std::complex<double> {
    auto params = PhysicalParameters::getInstance();

    double rrr   = (static_cast<int>(params->nx / 2) * params->dx);
    double cos_x = std::cos(2 * M_PI * (x - params->dd) / rrr);
    double cos_y = std::cos(2 * M_PI * y / rrr);
    double cos_z = std::cos(1 * M_PI * z / rrr);

    return std::complex<double>(cos_x * cos_y * cos_z, 0.);
});

REGISTER_INITIALIZER(GAUSS, [](double x, double y, double z) -> std::complex<double> {
    auto params    = PhysicalParameters::getInstance();
    double sigma_x = (params->nx * params->dx) / 10.;
    double sigma_y = (params->ny * params->dy) / 10.;
    double sigma_z = (params->nz * params->dz) / 10.;

    double val = std::exp(-0.5 * (x * x) / (sigma_x * sigma_x)) *
                 std::exp(-0.5 * (y * y) / (sigma_y * sigma_y)) *
                 std::exp(-0.5 * (z * z) / (sigma_z * sigma_z));
    return std::complex<double>(val, 0.);
});

REGISTER_INITIALIZER(MULTIPLE_GAUSS, [](double x, double y, double z) {
    auto params = PhysicalParameters::getInstance();
    auto sctx   = SimulationContext::getInstance();

    int n_maximas = params->n_gauss_max;
    std::vector<double> centers_x(n_maximas);
    std::vector<double> centers_y(n_maximas);

    for (int idx = 0; idx < n_maximas; idx++) {
        centers_x[idx] = idx % 2 ? params->dd : -params->dd;

        if (n_maximas % 2 == 1) {
            int center_y_idx = n_maximas / 2;
            double y_offset =
                (idx - center_y_idx) * (params->ny * params->dy) / 2. / (n_maximas + 1.);
            centers_y[idx] = y_offset;
        }

        if (n_maximas % 2 == 0) {
            int y_idx       = (idx / 2 + 1);
            int y_maximas   = (n_maximas / 2 + 1);
            double y_offset = y_idx * (params->ny * params->dy) / y_maximas;
            centers_y[idx]  = sctx->get_y(0) + y_offset;
        }
    }

    double sigma_x = (params->nx * params->dx) / 20.;
    double sigma_y = (params->ny * params->dy) / 20.;
    double sigma_z = (params->nz * params->dz) / 20.;

    double val = 0.0;
    for (int m = 0; m < n_maximas; m++) {
        val += std::exp(-0.5 * std::pow(x - centers_x[m], 2) / (sigma_x * sigma_x)) *
               std::exp(-0.5 * std::pow(y - centers_y[m], 2) / (sigma_y * sigma_y)) *
               std::exp(-0.5 * z * z / (sigma_z * sigma_z));
    }
    return std::complex<double>(val, 0.);
});

REGISTER_INITIALIZER(
    SETUP_GAUSS, ([](double x, double y, double z) {
        auto params = PhysicalParameters::getInstance();

        int n_maximas = params->bec_droplets_x * params->bec_droplets_y * params->bec_droplets_z;
        std::vector<std::tuple<double, double, double>> centers(n_maximas);

        double x_tot = params->nx * params->dx;
        double y_tot = params->ny * params->dy;
        double z_tot = params->nz * params->dz;

        for (int i = 0; i < params->bec_droplets_x; i++) {
            for (int j = 0; j < params->bec_droplets_y; j++) {
                for (int k = 0; k < params->bec_droplets_z; k++) {
                    int idx = i * params->bec_droplets_y * params->bec_droplets_z +
                              j * params->bec_droplets_z + k;

                    double spacing_x = x_tot / static_cast<double>(params->bec_droplets_x + 1);
                    double spacing_y = y_tot / static_cast<double>(params->bec_droplets_y + 1);
                    double spacing_z = z_tot / static_cast<double>(params->bec_droplets_z + 1);

                    double cx = -x_tot / 2 + spacing_x * (i + 1);
                    double cy = -y_tot / 2 + spacing_y * (j + 1);
                    double cz = -z_tot / 2 + spacing_z * (k + 1);

                    centers[idx] = {cx, cy, cz};
                }
            }
        }

        double sigma_x = (params->nx * params->dx) / 20.;
        double sigma_y = (params->ny * params->dy) / 20.;
        double sigma_z = (params->nz * params->dz) / 20.;

        double val = 0.0;
        for (int m = 0; m < n_maximas; m++) {
            val += std::exp(-0.5 * std::pow(x - std::get<0>(centers[m]), 2) / (sigma_x * sigma_x)) *
                   std::exp(-0.5 * std::pow(y - std::get<1>(centers[m]), 2) / (sigma_y * sigma_y)) *
                   std::exp(-0.5 * std::pow(z - std::get<2>(centers[m]), 2) / (sigma_z * sigma_z));
        }
        return std::complex<double>(val, 0.);
    }));

//! Initialize with Gaussians on a circle with R = (lx + ly) / 20
REGISTER_INITIALIZER(CYLINDRICAL_GAUSS, [](double x, double y, double z) {
    auto params = PhysicalParameters::getInstance();

    int n_maximas = params->n_gauss_max;
    std::vector<double> centers_x(n_maximas);
    std::vector<double> centers_y(n_maximas);

    double dtheta = 2. * M_PI / n_maximas;
    double lx     = params->dx * params->nx;
    double ly     = params->dy * params->ny;
    double r      = (lx + ly / 2.) / 10.;

    for (int idx = 0; idx < n_maximas; idx++) {
        double theta = dtheta * idx;
        double x     = r * cos(theta);
        double y     = r * sin(theta);

        centers_x[idx] = x;
        centers_y[idx] = y;
    }

    double sigma_x = (params->nx * params->dx) / 20.;
    double sigma_y = (params->ny * params->dy) / 20.;
    double sigma_z = (params->nz * params->dz) / 20.;

    double val = 0.0;
    for (int m = 0; m < n_maximas; m++) {
        val += std::exp(-0.5 * std::pow(x - centers_x[m], 2) / (sigma_x * sigma_x)) *
               std::exp(-0.5 * std::pow(y - centers_y[m], 2) / (sigma_y * sigma_y)) *
               std::exp(-0.5 * z * z / (sigma_z * sigma_z));
    }
    return std::complex<double>(val, 0.);
});

//! Initialize with the centers inside the circle with R = (lx + ly) / 20.
REGISTER_INITIALIZER(RANDOM_GAUSS, [](double x, double y, double z) {
    auto params = PhysicalParameters::getInstance();

    int n_maximas = params->n_gauss_max;
    static bool init = true;

    static std::vector<double> centers_x(n_maximas);
    static std::vector<double> centers_y(n_maximas);
    if(init){
        double lx     = params->dx * params->nx;
        double ly     = params->dy * params->ny;
        double R      = (lx + ly / 2.) / 10.;

        for (int idx = 0; idx < n_maximas; idx++) {
            double r = rand() / (RAND_MAX + 1.) * R;
            double theta = rand() / (RAND_MAX + 1.) * 2. * M_PI;

            double x = r * std::cos(theta);
            double y = r * std::sin(theta);

            centers_x[idx] = x;
            centers_y[idx] = y;
        }

        init = false;
    }

    double sigma_x = (params->nx * params->dx) / 20.;
    double sigma_y = (params->ny * params->dy) / 20.;
    double sigma_z = (params->nz * params->dz) / 20.;

    double val = 0.0;
    for (int m = 0; m < n_maximas; m++) {
        val += std::exp(-0.5 * std::pow(x - centers_x[m], 2) / (sigma_x * sigma_x)) *
               std::exp(-0.5 * std::pow(y - centers_y[m], 2) / (sigma_y * sigma_y)) *
               std::exp(-0.5 * z * z / (sigma_z * sigma_z));
    }
    return std::complex<double>(val, 0.);
});

//! TODO: workaround this structure, these are valid, however need different place to be
REGISTER_INITIALIZER(TEXT_FILE, [](double x, double y, double z) -> std::complex<double> {

    return std::complex<double>(0., 0.);
});

REGISTER_INITIALIZER(BINARY_FILE, [](double x, double y, double z) -> std::complex<double> {
    return std::complex<double>(0., 0.);
});

#endif
