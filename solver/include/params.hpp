#ifndef SOLVER_PARAMETERS_HPP
#define SOLVER_PARAMETERS_HPP

//! \todo Just meanwhile sollution, Rather will use blaze (CUDA parallelization, 3D arrays)
#include <armadillo>
using CMat3D = arma::Cube<std::complex<double>>;
using Mat3D  = arma::Cube<double>;

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

    //! Center of simulation in y direction. Barely used.
    double y0;

    //! Smallest trap frequency in system.
    double omega0;

    //! Mass of atom - using Dysprosium mass;
    double m;

    //! Used in mexican hat potential in x direction.
    //! x^4 term
    double aa;

    //! Used in mexican hat potential in x direction.
    //! x^2 term
    double b;

    //! Number of nodes in calculations - x direction.
    int nx;

    //! Number of nodes in calculations - y direction.
    int ny;

    //! Number of nodes in calculations - z direction.
    int nz;

    //! Used for boundary condition.
    //!
    //! Original comment: rdy sluzy do liczenia warunku brzegowego, calkowania gestosci, ktora na
    //! brzegu jest 0 tak czy inaczej
    Mat3D rdy;

private:
    void init_values();
    void init_containers();

    void init_rdy();

    PhysicalParameters() {
        init_values();
        init_containers();
    };
    static PhysicalParameters *instance;
};

#endif
