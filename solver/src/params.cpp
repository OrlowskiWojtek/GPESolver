#include "include/params.hpp"
#include "units.hpp"
#include <cassert>
#include <cmath>

PhysicalParameters *PhysicalParameters::instance = nullptr;

void PhysicalParameters::init_values() {
    wzl = 120 * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - z direction
    wrl = 60. * 4.1356e-12 / 27211.6; // angular frequency of harmonic potential - y direction

    n_atoms   = 5e4;
    double dd = 1500 / 0.05292;                        // TODO: do not know what it is - used in aa
    m         = UnitConverter::mass_Da_to_au(163.929); // mass of Erb 164

    aa = m * std::pow(wrl, 2) / 4. / std::pow(dd, 2);
    b  = 0.5 * m * std::pow(wrl, 2);

    omega0 = 0;
    add    = 131;
    cdd    = 12 * M_PI * add / m;

    double edd = 1.5;       // ![1] zgodnie z podpisem powinno być dla jednostek zredukowanych
    double a   = add / edd; // ![1]

    gpp11 = 4. * M_PI * a / m;
    gamma = 128 * std::sqrt(M_PI) * std::pow(a, 2.5) / 3. / m * (1. + 1.5 * std::pow(edd,2));

    // 81 nodes in one direction + one central node
    nx = 81 * 2 + 1;
    ny = 81 * 2 + 1;
    nz = 12 * 2 + 1;

    dx = UnitConverter::len_nm_to_au(60);
    dy = UnitConverter::len_nm_to_au(60);
    dz = UnitConverter::len_nm_to_au(350);

    dxdydz = dx * dy * dz;

    assert(nx % 2 == 1);
    assert(ny % 2 == 1);
    assert(nz % 2 == 1);

    init_r();
}

void PhysicalParameters::init_r() {
    r_matrix.resize(nx, ny, nz);

    x_vec.resize(nx);
    y_vec.resize(ny);
    z_vec.resize(nz);

    for (int i = 0; i < nx; i++) {
        double x = (i - static_cast<int>(nx / 2.) + 1) * dx;
        x_vec[i] = x;
        for (int j = 0; j < ny; j++) {
            double y = (j - static_cast<int>(ny / 2.) + 1) * dy;
            y_vec[j] = y;
            for (int k = 0; k < nz; k++) {
                double z          = (k - static_cast<int>(nz / 2.) + 1) * dz;
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
