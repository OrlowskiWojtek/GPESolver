#include "solver/cpu_solver/cpu_it_solver.hpp"

CpuITGrossPitaevskiSolver::CpuITGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator) {
}

void CpuITGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data.allocate(nx, ny, nz);
}

void CpuITGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver->execute();
}

void CpuITGrossPitaevskiSolver::calc_norm() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    xnorma = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                xnorma += std::norm(m_data.cpsi(i, j, k));
            }
        }
    }

    xnorma *= params->get_dxdydz();
}

void CpuITGrossPitaevskiSolver::normalize() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                m_data.cpsi(i, j, k) /= std::sqrt(xnorma);
            }
        }
    }
}

void CpuITGrossPitaevskiSolver::imag_iter_linear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    wavefunction_t &cpsi  = m_data.cpsi;
    wavefunction_t &cpsii = m_data.cpsii;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double v = m_data.pote(i, j, k);
                std::complex<double> c1 =
                    -0.5 / (params->m * std::pow(params->dx, 2)) *
                        (cpsi(i - 1, j, k) + cpsi(i + 1, j, k) - 2. * cpsi(i, j, k)) -
                    0.5 / (params->m * std::pow(params->dy, 2)) *
                        (cpsi(i, j - 1, k) + cpsi(i, j + 1, k) - 2. * cpsi(i, j, k)) -
                    0.5 / (params->m * std::pow(params->dz, 2)) *
                        (cpsi(i, j, k - 1) + cpsi(i, j, k + 1) - 2. * cpsi(i, j, k)) +
                    cpsi(i, j, k) * (v + params->cdd * m_data.fi3d(i, j, k));
                cpsii(i, j, k) = cpsi(i, j, k) - params->time_dt * c1;
            }
        }
    }
}

void CpuITGrossPitaevskiSolver::imag_iter_nonlinear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    wavefunction_t &cpsi  = m_data.cpsi;
    wavefunction_t &cpsii = m_data.cpsii;
    double w              = params->n_atoms;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                cpsii(i, j, k) = cpsii(i, j, k) -
                                 params->time_dt *
                                     (params->ggp11 * std::norm(cpsi(i, j, k)) * cpsi(i, j, k) * w +
                                      params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
                                          cpsi(i, j, k) * params->w_15);
            }
        }
    }

    cpsi = cpsii;
}

void CpuITGrossPitaevskiSolver::calc_energy() {
    ene.e_kin = 0.;
    ene.e_pot = 0.;
    ene.e_int = 0.;
    ene.e_ext = 0.;
    ene.e_bmf = 0.;

    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    std::complex<double> grad_psi_x;
    std::complex<double> grad_psi_y;
    std::complex<double> grad_psi_z;

    wavefunction_t &cpsi = m_data.cpsi;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                // Kinetic energy
                grad_psi_x = -(cpsi(i + 1, j, k) + cpsi(i - 1, j, k) - 2. * cpsi(i, j, k)) /
                             (std::pow(params->dx, 2));
                grad_psi_y = -(cpsi(i, j + 1, k) + cpsi(i, j - 1, k) - 2. * cpsi(i, j, k)) /
                             (std::pow(params->dy, 2));
                grad_psi_z = -(cpsi(i, j, k + 1) + cpsi(i, j, k - 1) - 2. * cpsi(i, j, k)) /
                             (std::pow(params->dz, 2));
                ene.e_kin +=
                    ((grad_psi_x + grad_psi_y + grad_psi_z) * std::conj(cpsi(i, j, k))).real();

                // Potential energy
                ene.e_pot += m_data.pote(i, j, k) * std::norm(cpsi(i, j, k));

                // Interaction energy
                ene.e_int +=
                    0.5 * params->ggp11 * std::norm(cpsi(i, j, k)) * std::norm(cpsi(i, j, k));

                // Dipole-dipole interaction energy
                ene.e_ext += 0.5 * params->cdd * m_data.fi3d(i, j, k) * std::norm(cpsi(i, j, k)) *
                             params->n_atoms;

                // beyond mean-field energy
                ene.e_bmf += 2. / 5. * params->gamma * std::pow(std::abs(cpsi(i, j, k)), 5);
            }
        }
    }

    ene.e_kin *= params->get_dxdydz() / (2 * params->m) * params->n_atoms;
    ene.e_pot *= params->get_dxdydz() * params->n_atoms;
    ene.e_int *= params->get_dxdydz() * std::pow(params->n_atoms, 2);
    ene.e_ext *= params->get_dxdydz();
    ene.e_bmf *= params->get_dxdydz() * std::pow(params->n_atoms, 2.5);

    ene.sum();

    enes.emplace_back(ene);
}

void CpuITGrossPitaevskiSolver::prepare_fft() {
    poisson_solver = std::make_unique<FFTWPoissonSolver>(&m_data.cpsi, &m_data.fi3d);
};

void CpuITGrossPitaevskiSolver::import_pote() {
    m_data.pote = buf_data->pote;
};

void CpuITGrossPitaevskiSolver::import_data() {
    m_data.cpsi  = buf_data->cpsi;
    m_data.cpsii = buf_data->cpsii;
};

void CpuITGrossPitaevskiSolver::export_data() {
    buf_data->cpsi = m_data.cpsi;
};

const int CpuITGrossPitaevskiSolver::iter_per_summary() const {
    return 100;
}

void CpuITGrossPitaevskiSolver::iterate() {
    calc_fi3d();
    imag_iter_linear_step();
    imag_iter_nonlinear_step();
    calc_norm();
    normalize();
}

// Do not need to adjust these calculations
void CpuITGrossPitaevskiSolver::adjust(int /*iter*/) {
}

void CpuITGrossPitaevskiSolver::finish() {
    export_data();
    p_mediator->save_initial_state(buf_data->cpsi);
}
