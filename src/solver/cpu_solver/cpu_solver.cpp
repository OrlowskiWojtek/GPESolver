#include "solver/cpu_solver/cpu_solver.hpp"

CpuGrossPitaevskiSolver::CpuGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator) {
}

void CpuGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data.allocate(nx, ny, nz);
}

void CpuGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver->execute();
}

void CpuGrossPitaevskiSolver::calc_norm() {
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

void CpuGrossPitaevskiSolver::normalize() {
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

void CpuGrossPitaevskiSolver::imag_iter_linear_step() {
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
                cpsii(i, j, k) = cpsi(i, j, k) - params->imag_time_dt * c1;
            }
        }
    }
}

void CpuGrossPitaevskiSolver::imag_iter_nonlinear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    wavefunction_t &cpsi  = m_data.cpsi;
    wavefunction_t &cpsii = m_data.cpsii;
    double w              = params->n_atoms;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                cpsii(i, j, k) =
                    cpsii(i, j, k) -
                    params->imag_time_dt * ((params->ggp11 - params->cdd / 3) *
                                                std::norm(cpsi(i, j, k)) * cpsi(i, j, k) * w +
                                            params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
                                                cpsi(i, j, k) * params->w_15);
            }
        }
    }

    cpsi = cpsii;
}

void CpuGrossPitaevskiSolver::real_fft_potential_half_step() {
    int nx           = params->nx;
    int ny           = params->ny;
    int nz           = params->nz;
    double w         = params->n_atoms;
    double dt_factor = params->real_time_dt / 2.;

    wavefunction_t &cpsi = m_data.cpsi;
    double psi_re;
    double psi_im;
    double s;
    double c;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double v_ext   = m_data.pote(i, j, k);
                double density = std::norm(cpsi(i, j, k));

                double v_int = (params->ggp11 - params->cdd / 3) * density * w +
                               params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) * params->w_15;

                double total_potential = v_ext + params->cdd * m_data.fi3d(i, j, k) + v_int;

                psi_re = cpsi(i, j, k).real();
                psi_im = cpsi(i, j, k).imag();
                sincos(-dt_factor * total_potential, &s, &c);

                cpsi(i, j, k).imag(psi_re * s + psi_im * c);
                cpsi(i, j, k).real(psi_re * c - psi_im * s);
            }
        }
    }
}

void CpuGrossPitaevskiSolver::real_fft_kinetic_step() {
    rt_split_solver->execute();
}

void CpuGrossPitaevskiSolver::calc_energy() {
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
                ene.e_ext -= params->cdd / 3. * std::pow(std::norm(cpsi(i, j, k)), 2) / 2. *
                             std::pow(params->n_atoms, 2);

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

void CpuGrossPitaevskiSolver::prepare_fft() {
    poisson_solver  = std::make_unique<FFTWPoissonSolver>(&m_data.cpsi, &m_data.fi3d, &m_data.pote);
    rt_split_solver = std::make_unique<FFTWRealTimeSplitSolver>(&m_data.cpsi, &m_data.fi3d, &m_data.pote);
};

void CpuGrossPitaevskiSolver::import_pote() {
    m_data.pote = buf_data->pote;
};

void CpuGrossPitaevskiSolver::import_data() {
    m_data.cpsi = buf_data->cpsi;
    m_data.cpsii = buf_data->cpsii;
};

void CpuGrossPitaevskiSolver::export_data() {
    buf_data->cpsi = m_data.cpsi;
};
