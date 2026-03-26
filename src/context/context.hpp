#ifndef SIMULATION_CONTEXT_HPP
#define SIMULATION_CONTEXT_HPP

#include "mat3d/stdmat3d.hpp"
#include <complex>
#include <parameters/parameters.hpp>
#include <vector>

//! struct BECEnergies.
//!
//! Stores energies of Bose-Einstein Condensate
struct BECEnergies {
    double e_kin;   //!< Kinetic energy
    double e_pot;   //!< Potential energy
    double e_int;   //!< Interaction energy
    double e_ext;   //!< Dipole-dipole interaction energy
    double e_bmf;   //!< Beyond mean-field energy
    double e_total; //!< Total energy

    void sum() {
        e_total = e_kin + e_pot + e_int + e_ext + e_bmf;
    }
};

using wavefunction_t       = StdMat3D<std::complex<double>>;
using potential_t          = StdMat3D<double>;
using energies_container_t = std::vector<BECEnergies>;
using energies_t           = BECEnergies;

//! Class SimulationContext.
//!
//! \brief  used to store basic simulation containters like arrays of grid spaces;
class SimulationContext {
public:
    static SimulationContext* getInstance();

    // void calc_norm();
    // void normalize();
    // void init_potential();
    // void free_potential_well();
    // double pote_value(int ix, int iy, int iz);
    // double pote_released_value(int ix, int iy, int iz);
    // PhysicalParameters* params;
    // StdMat3D<std::complex<double>> cpsi;
    // StdMat3D<std::complex<double>> cpsii;
    // StdMat3D<double> potential;
    // StdMat3D<double> fi3d;

    //! Initialize vectors based on loaded parameters.
    void initialize();

    //! Return value of x at given index.
    inline double get_x(int ix) const {
        return x_vec[ix];
    }
    //! Return value of y at given index.
    inline double get_y(int iy) const {
        return y_vec[iy];
    }
    //! Return value of z at given index.
    inline double get_z(int iz) const {
        return z_vec[iz];
    }

private:
    PhysicalParameters *params;

    std::vector<double> x_vec;
    std::vector<double> y_vec;
    std::vector<double> z_vec;

    SimulationContext();
    static SimulationContext* instance;
};

#endif
