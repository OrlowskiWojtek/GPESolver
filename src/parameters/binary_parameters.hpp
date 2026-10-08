#ifndef BINARY_MIXTURE_PARAMETERS_HPP
#define BINARY_MIXTURE_PARAMETERS_HPP

struct BinaryMixtureParameters {
    BinaryMixtureParameters(const BinaryMixtureParameters &)            = delete;
    BinaryMixtureParameters &operator=(const BinaryMixtureParameters &) = delete;

    static BinaryMixtureParameters *getInstance() {
        if (!instance) {
            instance = new BinaryMixtureParameters();
        }

        return instance;
    }

    //! Mass of atom;
    double m_a;
    double m_b;

    // number of atoms
    double n_atoms_a;
    double n_atoms_b;

    double g11, g22, g12;       // contact interaction
    // cdd12 = \mu_0 \mu_a \mu_b
    double cdd11, cdd22, cdd12;   // dipole strengths

    //! Adds random noise at the beginning of calculations
    bool add_random_noise = false;

    void set_to_default();
private:

    BinaryMixtureParameters() {};

    static BinaryMixtureParameters *instance;
};

#endif
