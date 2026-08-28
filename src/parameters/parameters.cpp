#include "parameters/parameters.hpp"
#include "output.hpp"
#include "units.hpp"
#include <cmath>

PhysicalParameters *PhysicalParameters::instance = nullptr;

void PhysicalParameters::set_default_values() {
    n_atoms = 4e4;
    m       = UnitConverter::mass_Da_to_au(163.929); // mass of Erb 164

    dd = UnitConverter::len_nm_to_au(1500.0);

    omega_x = 60. * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - x direction
    omega_y = 60. * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - y direction
    omega_z =
        120 * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - z direction fixed

    nx = 40 * 2 + 1;
    ny = 40 * 2 + 1;
    nz = 20 * 2 + 1;

    iter_imag = 10000;
    iter_real = 300000;

    edd = 1.45;

    dx = UnitConverter::len_nm_to_au(150);
    dy = UnitConverter::len_nm_to_au(150);
    dz = UnitConverter::len_nm_to_au(200);

    fftw_n_threads = 4;

    init_parameters();
}

void PhysicalParameters::init_parameters() {
    aa = m * std::pow(omega_x, 2) / 4. / std::pow(dd, 2);
    b  = 0.5 * m * std::pow(omega_x, 2);

    dxdydz = dx * dy * dz;

    add = 131.;
    cdd = 12. * M_PI * add / m;

    double a = add / edd; // ![1]

    ggp11 = 4. * M_PI * a / m;
    gamma = 128. * std::sqrt(M_PI) * std::pow(a, 2.5) / 3. / m * (1. + 1.5 * std::pow(edd, 2));
    w_15  = std::pow(n_atoms, 1.5);


    double lx = nx * dx;
    double ly = ny * dy;
    double lz = nz * dz;
    Rc = std::max(lx, std::max(ly, lz));

    assert(Rc < (2 * lx) && Rc < (2 * ly) && Rc < (2 * lz));
}

void PhysicalParameters::update_edd(uint iter) {
    edd = edd_start + (edd_stop - edd_start) * iter / iter_real;

    double a = add / edd;
    ggp11    = 4 * M_PI * a / m;
    gamma    = 128. * std::sqrt(M_PI) * std::pow(a, 2.5) / 3. / m * (1. + 1.5 * std::pow(edd, 2));
}

double PhysicalParameters::get_dxdydz() {
    return dxdydz;
}

void PhysicalParameters::print_initialization_parameters() {
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Initialization method:");
    OutputFormatter::printBoxedMessage(init_strategy.to_string());
    if (init_strategy.type == InitializationOption::Type::MULTIPLE_GAUSS) {
        OutputFormatter::printBoxedMessage("Number of initial maximas:");
        OutputFormatter::printBoxedMessage(n_gauss_max);
    }
    if (init_strategy.type == InitializationOption::Type::FROM_BINARY_FILE) {
        OutputFormatter::printBoxedMessage("Binary file name:");
        OutputFormatter::printBoxedMessage(load_filename);
    }
    if (init_strategy.type == InitializationOption::Type::FROM_TEXT_FILE) {
        OutputFormatter::printBoxedMessage("Text file name:");
        OutputFormatter::printBoxedMessage(load_filename);
    }
    if (init_strategy.type == InitializationOption::Type::SETUP_GAUSS) {
        OutputFormatter::printBoxedMessage("Gaussians configuration:");
        OutputFormatter::printBoxedMessage(std::to_string(bec_droplets_x) + " x " +
                                           std::to_string(bec_droplets_y) + " x " +
                                           std::to_string(bec_droplets_z));
    }

    OutputFormatter::printBorderLine();
}

void PhysicalParameters::print_box_parameters() {
    OutputFormatter::printBorderLine();

    OutputFormatter::printBoxedMessage("COMPUTATIONAL BOX PARAMETERS");

    OutputFormatter::printBoxedMessage("Grid size");
    OutputFormatter::printBoxedMessage("nx: ", nx);
    OutputFormatter::printBoxedMessage("ny: ", ny);
    OutputFormatter::printBoxedMessage("nz: ", nz);

    OutputFormatter::printBoxedMessage("Grid spacing");
    OutputFormatter::printBoxedMessage("dx (nm): ", UnitConverter::len_au_to_nm(dx));
    OutputFormatter::printBoxedMessage("dy (nm): ", UnitConverter::len_au_to_nm(dy));
    OutputFormatter::printBoxedMessage("dz (nm): ", UnitConverter::len_au_to_nm(dz));

    OutputFormatter::printBoxedMessage("total grid size:");
    OutputFormatter::printBoxedMessage("x (nm): ", UnitConverter::len_au_to_nm(dx) * nx);
    OutputFormatter::printBoxedMessage("y (nm): ", UnitConverter::len_au_to_nm(dy) * ny);
    OutputFormatter::printBoxedMessage("z (nm): ", UnitConverter::len_au_to_nm(dz) * nz);

    OutputFormatter::printBorderLine();
}

//! \todo: make dependency like parameters -> not need to print all parameters all the time
void PhysicalParameters::print_simulation_parameters() {
    OutputFormatter::printBorderLine();

    OutputFormatter::printBoxedMessage("SIMULATION PARAMETERS");

    OutputFormatter::printBoxedMessage("Mass (Da): ", UnitConverter::mass_au_to_Da(m));
    OutputFormatter::printBoxedMessage("Number of atoms: ", n_atoms);
    OutputFormatter::printBoxedMessage("Scattering length (nm): ",
                                       UnitConverter::len_au_to_nm(add));

    OutputFormatter::printBoxedMessage("Dipole-dipole interaction");
    OutputFormatter::printBoxedMessage("Cdd (au): ", cdd);
    OutputFormatter::printBoxedMessage("epsilon_dd: ", edd);
    OutputFormatter::printBoxedMessage("epsilon_dd_start: ", edd_start);
    OutputFormatter::printBoxedMessage("epsilon_dd_stop: ", edd_stop);

    OutputFormatter::printBoxedMessage("ITERATIONS");
    OutputFormatter::printBoxedMessage("IMAGINARY -  " + std::to_string(iter_imag));
    OutputFormatter::printBoxedMessage("REAL -  " + std::to_string(iter_real));

    std::stringstream ss;
    ss << std::scientific << imag_time_dt;
    OutputFormatter::printBoxedMessage("imag_dt -  " + ss.str());
    ss.str(std::string());
    ss << std::scientific << real_time_dt;
    OutputFormatter::printBoxedMessage("real_dt -  " + ss.str());
    
    OutputFormatter::printBoxedMessage("Calculation strategy:", calc_strategy.to_string());

    OutputFormatter::printBorderLine();
}

void PhysicalParameters::print_potential_parameters(){
    OutputFormatter::printBorderLine();

    OutputFormatter::printBoxedMessage("POTENTIAL_PARAMETERS");

    OutputFormatter::printBoxedMessage("Potential:", pote_key);
    OutputFormatter::printBoxedMessage("Trap frequencies:");

    OutputFormatter::printBoxedMessage("omega_x (Hz): ", UnitConverter::freq_au_to_Hz(omega_x));
    OutputFormatter::printBoxedMessage("omega_y (Hz): ", UnitConverter::freq_au_to_Hz(omega_y));
    OutputFormatter::printBoxedMessage("omega_z (Hz): ", UnitConverter::freq_au_to_Hz(omega_z));

    OutputFormatter::printBoxedMessage("Potential minimas dd: ", UnitConverter::len_au_to_nm(dd));

    OutputFormatter::printBorderLine();
}

void PhysicalParameters::print() {
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("GPE SOLVER PARAMETERS");
    OutputFormatter::printBorderLine();

    print_initialization_parameters();
    print_simulation_parameters();
    print_potential_parameters();
    print_box_parameters();
}
