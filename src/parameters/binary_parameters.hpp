#ifndef BINARY_MIXTURE_PARAMETERS_HPP
#define BINARY_MIXTURE_PARAMETERS_HPP

#include <cmath>

//! \brief SpeciesPamateres for single species of atoms.
//! All in atomic units.
struct SpeciesParameters {
    double mass            = 0.;
    int n_atoms            = 0.;
    double magnetic_moment = 0.0;
};

struct BinaryMixtureParameters {
    BinaryMixtureParameters(const BinaryMixtureParameters &) = delete;
    BinaryMixtureParameters &
    operator=(const BinaryMixtureParameters &) = delete;

    static BinaryMixtureParameters *getInstance() {
        if (!instance) {
            instance = new BinaryMixtureParameters();
        }

        return instance;
    }

    SpeciesParameters spec_a;
    SpeciesParameters spec_b;

    double s11 = 0.0; // Scattering lengths in a0 (atomic units of length)
    double s22 = 0.0;
    double s12 = 0.0;

    double g11 = 0.0;
    double g22 = 0.0;
    double g12 = 0.0;

    double cdd11 = 0.0;
    double cdd22 = 0.0;
    double cdd12 = 0.0;

    // Precompute factors
    double four_third_pi2 = 4.0 / (3.0 * M_PI * M_PI);
    double factor_m1;
    double factor_m2;

    void set_to_default();
    void precompute();

private:
    BinaryMixtureParameters() {};

    static BinaryMixtureParameters *instance;
};

#endif
