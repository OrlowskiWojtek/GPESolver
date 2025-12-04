#ifndef SOLVER_CONTEXT_HPP
#define SOLVER_CONTEXT_HPP

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

#endif
