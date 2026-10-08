#include "parameters/binary_parameters.hpp"
#include "utils/units/units.hpp"
#include <cmath>

BinaryMixtureParameters *BinaryMixtureParameters::instance = nullptr;

void BinaryMixtureParameters::set_to_default() {
    m_a = UnitConverter::mass_Da_to_au(164.0);
    m_b = UnitConverter::mass_Da_to_au(162.0);

    n_atoms_a = 4.0e4;
    n_atoms_b = 4.0e4;

    // chyba chodzi o pole magnetyczne, reszte przepisuje z parametryzacji z fortrana
    // https://arxiv.org/pdf/1803.10676
    const double b     = 26.0;
    const double a164  = 91.0 * (1.0 - 30.9 / (b - 76.9) - 23.6 / (b - 178.8));
    const double a162  = 220.0 * (1.0 - 1.9 / (b - 21.91) - 0.14 / (b - 26.902));
    const double a6264 = 105.0 * (1.0 - 1.0 / (b - 10.8));

    const double reduced_mass = 2.0 * m_a * m_b / (m_a + m_b);

    g11 = 4.0 * M_PI * a164 / m_a;
    g22 = 4.0 * M_PI * a162 / m_b;
    g12 = 4.0 * M_PI * a6264 / reduced_mass;

    //! to think - what is a material parameter and what is depending on mass
    cdd11 = 12.0 * M_PI * 131.0 / m_a;
    cdd22 = 12.0 * M_PI * 131.0 / m_b;
    cdd12 = 12.0 * M_PI * 131.0 / reduced_mass;
}
