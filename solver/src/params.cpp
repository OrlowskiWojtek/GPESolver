#include "include/params.hpp"
#include "include/output.hpp"
#include "units.hpp"
#include <cassert>
#include <cmath>

PhysicalParameters *PhysicalParameters::instance = nullptr;

void PhysicalParameters::set_default_values() {
    n_atoms = 4e4;
    m       = UnitConverter::mass_Da_to_au(163.929); // mass of Erb 164

    dd = UnitConverter::len_nm_to_au(1500.0);

    nx = 40 * 2 + 1;
    ny = 40 * 2 + 1;
    nz = 20 * 2 + 1;

    edd = 1.45;
    load_initial_state = false;

    dx = UnitConverter::len_nm_to_au(150);
    dy = UnitConverter::len_nm_to_au(150);
    dz = UnitConverter::len_nm_to_au(200);

    init_parameters();
}

void PhysicalParameters::init_parameters(){
    assert(nx % 2 == 1);
    assert(ny % 2 == 1);
    assert(nz % 2 == 1);

    wzl = 120 * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - z direction
    wrl = 60. * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - y direction
    
    aa = m * std::pow(wrl, 2) / 4. / std::pow(dd, 2);
    b  = 0.5 * m * std::pow(wrl, 2);

    dxdydz = dx * dy * dz;

    add = 131.;
    cdd = 12. * M_PI * add / m;

    double a   = add / edd; // ![1]

    ggp11 = 4. * M_PI * a / m;
    gamma = 128. * std::sqrt(M_PI) * std::pow(a, 2.5) / 3. / m * (1. + 1.5 * std::pow(edd, 2));

    init_r();
}

void PhysicalParameters::init_r() {
    r_matrix.resize(nx, ny, nz);

    x_vec.resize(nx);
    y_vec.resize(ny);
    z_vec.resize(nz);

    for (int i = 0; i < nx; i++) {
        double x = (i - (static_cast<int>(nx / 2.) + 1)) * dx;
        x_vec[i] = x;
        for (int j = 0; j < ny; j++) {
            double y = (j - (static_cast<int>(ny / 2.) + 1)) * dy;
            y_vec[j] = y;
            for (int k = 0; k < nz; k++) {
                double z          = (k - (static_cast<int>(nz / 2.) + 1)) * dz;
                r_matrix(i, j, k) = std::sqrt(x * x + y * y + z * z);
                z_vec[k]          = z;
            }
        }
    }
}

double PhysicalParameters::get_x(int ix) {
    return x_vec[ix];
}

double PhysicalParameters::get_y(int iy) {
    return y_vec[iy];
}

double PhysicalParameters::get_z(int iz) {
    return z_vec[iz];
}

double PhysicalParameters::get_dxdydz() {
    return dxdydz;
}

double PhysicalParameters::get_r(int ix, int iy, int iz) {
    return r_matrix(ix, iy, iz);
}

void PhysicalParameters::print() {
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Parameters");
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Load initial state:");
    OutputFormatter::printBoxedMessage(load_initial_state ? "Yes" : "No");
    OutputFormatter::printBoxedMessage("Calculation strategy:");
    OutputFormatter::printBoxedMessage(calc_strategy.to_string());
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Mass (Da): ", UnitConverter::mass_au_to_Da(m));
    OutputFormatter::printBoxedMessage("Number of atoms: ", n_atoms);
    OutputFormatter::printBoxedMessage("Scattering length (nm): ", UnitConverter::len_au_to_nm(add));
    OutputFormatter::printBoxedMessage("Potential minimas dd: ", UnitConverter::len_au_to_nm(dd));
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Dipole-dipole interaction");
    OutputFormatter::printBoxedMessage("Cdd (au): ", cdd);
    OutputFormatter::printBoxedMessage("epsilon_dd: ", edd);
    OutputFormatter::printBorderLine();
    OutputFormatter::printBoxedMessage("Trap frequencies:");
    OutputFormatter::printBoxedMessage("wrl (y direction) (au): ", wrl);
    OutputFormatter::printBoxedMessage("wzl (z direction) (au): ", wzl);

    OutputFormatter::printBorderLine();
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
