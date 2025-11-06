#include "include/solver.hpp"
#include "include/numerical_params.hpp"
#include "include/params.hpp"

GrossPitaevskiSolver::GrossPitaevskiSolver()
    : params(PhysicalParameters::getInstance()) {

    init_containers();
}

void GrossPitaevskiSolver::solve() {
    calc_initial_state();
}

void GrossPitaevskiSolver::calc_initial_state() {
    for (size_t iter = 0; iter < NumericalParameters::iter_imag_evo; iter++) {
        std::cout << iter << std::endl;
        imag_time_iter();
    }
}

void GrossPitaevskiSolver::imag_time_iter() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    // todo: push into double psi_norm();
    double xnorma = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                xnorma += std::norm(cpsi(i, j, k)) * params->get_dxdydz();
            }
        }
    }

    // todo
    calc_fi3d();

    //! WARNING: assumption dx = dy
    // \todo hide into imaginary step
    for (int i = 1; i < nx -1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double v = pote(i, j, k);
                std::complex<double> c1 =
                    -0.5 / params->m / std::pow(params->dx, 2) *
                        (cpsi(i - 1, j, k) + cpsi(i + 1, j, k) + cpsi(i, j - 1, k) +
                         cpsi(i, j + 1, k) - 4. * cpsi(i, j, k)) -
                    0.5 / params->m / std::pow(params->dz, 2) *
                        (cpsi(i, j, k - 1) + cpsi(i, j, k + 1) - 2. * cpsi(i, j, k)) +
                    cpsi(i, j, k) * (v + params->cdd * fi3d(i, j, k));
                cpsii(i, j, k) = cpsi(i, j, k) - NumericalParameters::imag_time_dt * c1;
            }
        }
    }

    double w = params->n_atoms / xnorma;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                cpsii(i, j, k) =
                    cpsii(i, j, k) - NumericalParameters::imag_time_dt *
                                         ((params->gpp11 - params->cdd / 3.) * std::norm(cpsi(i, j, k)) *
                                              cpsi(i, j, k) * (w) +
                                          params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
                                              cpsi(i, j, k) * std::pow(w, 1.5));
            }
        }
    }

    cpsi = cpsii;
}

void GrossPitaevskiSolver::init_containers() {
    init_potential();

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    cpsii.resize(nx, ny, nz);
    cpsi.resize(nx, ny, nz);
    cpsin.resize(nx, ny, nz);

    fi3do.resize(nx, ny, nz);
    fi3d.resize(nx, ny, nz);

    init_with_cos();
}

void GrossPitaevskiSolver::init_potential() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    pote.resize(nx, ny, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                pote(i, j, k) = pote_value(i, j, k);
            }
        }
    }
}

double GrossPitaevskiSolver::pote_value(int ix, int iy, int iz) {
    double x = params->get_x(ix);
    double y = params->get_y(iy);
    double z = params->get_z(iz);

    //! \todo: change wrl and wzl
    double vx = -params->b * std::pow(x, 2) + params->aa * std::pow(x, 4);
    double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->wrl, 2);
    double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->wzl, 2);

    return vx + vy + vz;
}

double GrossPitaevskiSolver::pote_released_value(int ix, int iy, int iz) {
    double x = params->get_x(ix);
    double y = params->get_y(iy);
    double z = params->get_z(iz);

    //! \todo: change wrl and wzl
    double vx = params->aa * std::pow(x, 4);
    double vy = 0.5 * params->m * std::pow(y, 2) * std::pow(params->wrl, 2);
    double vz = 0.5 * params->m * std::pow(z, 2) * std::pow(params->wzl, 2);

    return vx + vy + vz;
}

void GrossPitaevskiSolver::init_with_cos() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                double x = params->get_x(i);
                double y = params->get_y(j);
                double z = params->get_z(k);

                double rrr = params->get_x(nx - 1);

                double cos_x = std::cos(2 * M_PI * x / rrr);
                double cos_y = std::cos(2 * M_PI * y / rrr);
                double cos_z = std::cos(1 * M_PI * z / rrr);

                cpsi(i, j, k) = cos_x * cos_y * cos_z;
            }
        }
    }
}

void GrossPitaevskiSolver::calc_fi3d() {
}
