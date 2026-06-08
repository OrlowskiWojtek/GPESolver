#ifndef POTENTIALS_HPP
#define POTENTIALS_HPP

#include <string>
#include <functional>
#include <unordered_map>
#include <cmath>
#include "parameters/parameters.hpp"

class PotentialRegistry {
public:
    using RegisterKey = std::string;
    using ID = unsigned int;
    using pote_func_t = std::function<double(double, double, double)>;

    struct PotentialInfo {
        RegisterKey key;
        ID id;
        pote_func_t func;
    };

    static PotentialRegistry& instance() {
        static PotentialRegistry registry;
        return registry;
    }

    // Register a potential with key, id, and function
    void register_potential(RegisterKey key, ID id, pote_func_t func) {
        potentials_[key] = {key, id, func};
        id_map_[id] = key;
    }

    // Getters
    pote_func_t get_function(const RegisterKey& key) const {
        auto it = potentials_.find(key);

        return (it != potentials_.end()) ? it->second.func : potentials_.find("FREE")->second.func;
    }

    pote_func_t get_function(ID id) const {
        auto it = id_map_.find(id);
        return (it != id_map_.end()) ? get_function(it->second) : potentials_.find("FREE")->second.func;
    }

private:
    std::unordered_map<RegisterKey, PotentialInfo> potentials_;
    std::unordered_map<ID, RegisterKey> id_map_;
};

static unsigned int ID = 0;

#define REGISTER_POTENTIAL(KEY, FUNC) \
    static auto _registrar_##KEY = []() { \
        PotentialRegistry::instance().register_potential(#KEY, ID++, FUNC); \
        return 0.; \
    }();

//! Classic harmonic potential with V(x,y,z) = m / 2 * (\omega_x^2 * x + \omega_y^2 * y + \omega_z^2 * z)
REGISTER_POTENTIAL(REGULAR, [](double x, double y, double z){
        auto params = PhysicalParameters::getInstance();

        double vx = 0.5 * params->m * std::pow(x, 2) * std::pow(params->omega_x, 2);
        double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
        double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

        return vx + vy + vz;
});

//! Mexican hat potential with minimas separate by d
REGISTER_POTENTIAL(MEXICAN, [](double x, double y, double z){
        auto params = PhysicalParameters::getInstance();

        double vx = -params->b * std::pow(x, 2) + params->aa * std::pow(x, 4);
        double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->omega_y, 2);
        double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

        return vx + vy + vz;
});

//! Cylindrical potential with form V(r, z), where r = sqrt(x,y)
REGISTER_POTENTIAL(CYLINDRICAL, [](double x, double y, double z){
        auto params = PhysicalParameters::getInstance();

        double r = std::sqrt(x*x + y*y);

        double vr = 0.5 * params->m * std::pow(r, 2) * std::pow((params->omega_y + params->omega_x) / 2., 2);
        double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

        return vr + vz;
});

//! Bound droplets only in z plane;
REGISTER_POTENTIAL(FREE, [](double x, double y, double z){
        auto params = PhysicalParameters::getInstance();
        double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->omega_z, 2);

        return vz;
});


#endif
