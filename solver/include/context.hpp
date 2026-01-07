#ifndef SOLVER_CONTEXT_HPP
#define SOLVER_CONTEXT_HPP

#include "include/params.hpp"
#include "mat3d/stdmat3d.hpp"
#include <complex>

struct BECEnergies{
    double e_kin;   //!< Kinetic energy
    double e_pot;   //!< Potential energy
    double e_int;   //!< Interaction energy
    double e_ext;   //!< Dipole-dipole interaction energy
    double e_bmf;   //!< Beyond mean-field energy
    double e_total; //!< Total energy
    
    void sum(){
        e_total = e_kin + e_pot + e_int + e_ext + e_bmf;
    }
};


class SimulationContext {
public:
    SimulationContext();
    void init_containers();

    void init_with_cos();
    void init_with_gauss();
    void init_from_file();

    void calc_norm();
    void normalize();

    void init_potential();
    void free_potential_well();
    double pote_value(int ix, int iy, int iz);
    double pote_released_value(int ix, int iy, int iz);

    PhysicalParameters* params;
    StdMat3D<std::complex<double>> cpsi;
    StdMat3D<std::complex<double>> cpsii;
    StdMat3D<double> potential;
    StdMat3D<double> fi3d;

    double norm;
};

#endif
