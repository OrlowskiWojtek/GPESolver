#include "parameters/binary_parameters.hpp"
#include "utils/units/units.hpp"
#include <cmath>

BinaryMixtureParameters *BinaryMixtureParameters::instance = nullptr;

void BinaryMixtureParameters::set_to_default() {
    spec_a.mass = UnitConverter::mass_Da_to_au(164.0);
    spec_b.mass = UnitConverter::mass_Da_to_au(162.0);

    spec_a.n_atoms = 4.0e4;
    spec_b.n_atoms = 4.0e4;

    // chyba chodzi o pole magnetyczne, reszte przepisuje z parametryzacji z
    // fortrana https://arxiv.org/pdf/1803.10676
    const double b    = 26.0;
    const double a164 = 91.0 * (1.0 - 30.9 / (b - 76.9) - 23.6 / (b - 178.8));
    const double a162 = 220.0 * (1.0 - 1.9 / (b - 21.91) - 0.14 / (b - 26.902));
    const double a6264 = 105.0 * (1.0 - 1.0 / (b - 10.8));

    const double reduced_mass =
        spec_a.mass * spec_b.mass / (spec_a.mass + spec_b.mass);

    g11 = 4.0 * M_PI * a164 / spec_a.mass;
    g22 = 4.0 * M_PI * a162 / spec_b.mass;
    g12 = 4.0 * M_PI * a6264 / reduced_mass;

    //! to think - what is a material parameter and what is depending on mass
    cdd11 = 12.0 * M_PI * 131.0 / spec_a.mass;
    cdd22 = 12.0 * M_PI * 131.0 / spec_b.mass;
    cdd12 = 12.0 * M_PI * 131.0 / reduced_mass;

    precompute();
}

void BinaryMixtureParameters::precompute() {
    four_third_pi2 = 4.0 / (3.0 * M_PI * M_PI);
    // fine structure constant
    constexpr double alpha = 7.2973525693e-3;
    constexpr double mu_0 = 4 * M_PI * alpha * alpha;

    const double reduced_mass =
        spec_a.mass * spec_b.mass / (spec_a.mass + spec_b.mass);

    four_third_pi2 = 4.0 / (3.0 * M_PI * M_PI);
    factor_m1      = four_third_pi2 * std::pow(spec_a.mass, 1.5);
    factor_m2      = four_third_pi2 * std::pow(spec_b.mass, 1.5);

    g11 = 4.0 * M_PI * s11 / spec_a.mass;
    g22 = 4.0 * M_PI * s22 / spec_b.mass;
    g12 = 2.0 * M_PI * s12 / reduced_mass;

    cdd11 =
        mu_0 * spec_a.magnetic_moment * spec_a.magnetic_moment;
    cdd22 =
        mu_0 * spec_b.magnetic_moment * spec_b.magnetic_moment;
    cdd12 =
        mu_0 * spec_a.magnetic_moment * spec_b.magnetic_moment;
}
