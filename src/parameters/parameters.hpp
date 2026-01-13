#ifndef SOLVER_PARAMETERS_HPP
#define SOLVER_PARAMETERS_HPP

#include <vector>
#include "parameters/complex_param.hpp"

/*! Struct PhysicalParameters.
 *  \brief contains parameters of simulation.
 *
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
    // Distance between potential minima
    double dd;
    //! Used in mexican hat potential in x direction.
    //! x^2 term
    double b;

    double ggp11;
    double gamma;

    double cdd;
    double add;
    double edd;

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

    //! Number of gaussian maximas to initialize
    int n_gauss_max;

    // Number of threads used in FFTW calculations.
    int fftw_n_threads = 1;

    //! Calculation strategy (options)
    CalcStrategy calc_strategy;

    //! Initialization option
    InitializationOption init_strategy;

    double get_dxdydz();

    void print();
    void set_default_values();
    void init_parameters();

private:
    std::vector<double> x_vec;
    std::vector<double> y_vec;
    std::vector<double> z_vec;

    double dxdydz;

    PhysicalParameters() {};

    static PhysicalParameters *instance;
};

#endif
