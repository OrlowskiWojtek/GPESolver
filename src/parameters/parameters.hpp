#ifndef SOLVER_PARAMETERS_HPP
#define SOLVER_PARAMETERS_HPP

#include <array>
#include <string>
#include <vector>

struct CalcStrategy {
    enum class Type {
        IMAGINARY_TIME,
        REAL_TIME,
        FULL,
        SPEED_TEST
    };

    static const std::array<const char *, 4> TypeNames;

    std::string to_string() {
        return TypeNames[static_cast<size_t>(type)];
    }

    void from_string(const std::string &str) {
        for (size_t i = 0; i < TypeNames.size(); ++i) {
            if (str == TypeNames[i]) {
                type = static_cast<Type>(i);
                return;
            }
        }

        type = CalcStrategy::Type::FULL;
    }

    Type type = CalcStrategy::Type::FULL;
};

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

    //! Calculation strategy (options)
    CalcStrategy calc_strategy;

    //! Number of gaussian maximas to initialize
    int n_gauss_max;

    // Number of threads used in FFTW calculations.
    int fftw_n_threads = 1;

    bool load_initial_state;

    double get_x(int ix);
    double get_y(int iy);
    double get_z(int iz);
    double get_r(int ix, int iy, int iz);

    double get_dxdydz();

    void print();
    void set_default_values();
    void init_parameters();

private:
    void init_r();

    std::vector<double> x_vec;
    std::vector<double> y_vec;
    std::vector<double> z_vec;

    double dxdydz;

    PhysicalParameters() {};

    static PhysicalParameters *instance;
};

#endif
