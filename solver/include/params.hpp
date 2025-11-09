#ifndef SOLVER_PARAMETERS_HPP
#define SOLVER_PARAMETERS_HPP

//! \todo Just meanwhile sollution, Rather will use blaze (CUDA parallelization, 3D arrays)
#include <armadillo>
#include "mat3d/stdmat3d.hpp"

/*! Struct PhysicalParameters.
 *  \brief contains physical parameters of simulation.
 *
 *  For now original names of parameters are kept.
 *  References for used values:
 *  [1] 
 */
struct PhysicalParameters {
    PhysicalParameters(const PhysicalParameters &)            = delete;
    PhysicalParameters &operator=(const PhysicalParameters &) = delete;

    static PhysicalParameters *getInstance() {
        if (!instance) {
            instance = new PhysicalParameters();
        }

        return instance;
    }

    //! Oscilator omega in z direction.
    double wzl;
    //! Oscilator omega in y direction.
    double wrl;
    //! Smallest trap frequency in system.
    double omega0;
    //! Mass of atom;
    double m;
    //! Number of used atoms;
    double n_atoms;
    //! Used in mexican hat potential in x direction.
    //! x^4 term
    double aa;
    //! Used in mexican hat potential in x direction.
    //! x^2 term
    double b;

    double ggp11;
    double gamma;

    //! dipol dipol interaction constant
    double cdd;
    //! still do not know \todo
    double add;

    //! Number of nodes in calculations - x direction.
    int nx;
    //! Number of nodes in calculations - y direction.
    int ny;
    //! Number of nodes in calculations - z direction.
    int nz;

    //! Distance per node - x direction
    double dx;
    //! Distance per node - y direction.
    double dy;
    //! Distance per node - z direction.
    double dz;

    double get_x(int ix);
    double get_y(int iy);
    double get_z(int iz);
    double get_r(int ix, int iy, int iz);

    double get_dxdydz();
private:
    void init_values();
    void init_r();

    StdMat3D<double> r_matrix;
    std::vector<double> x_vec;
    std::vector<double> y_vec;
    std::vector<double> z_vec;

    double dxdydz;

    PhysicalParameters() {
        init_values();
    };

    static PhysicalParameters *instance;
};

#endif
