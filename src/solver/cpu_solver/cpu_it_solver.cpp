#include "solver/cpu_solver/cpu_it_solver.hpp"
#include <iostream>

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
    const int nx = params->nx;
    const int ny = params->ny;
    const int nz = params->nz;

    const auto *__restrict__ p_cpsi = m_data.cpsi.get_data_restrict();

    xnorma = 0.0;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                const int idx = i * ny * nz + j * nz + k;
                xnorma += std::norm(p_cpsi[idx]);
            }
        }
    }

    xnorma *= params->get_dxdydz();
}

void CpuITGrossPitaevskiSolver::normalize() {
    const int nx = params->nx;
    const int ny = params->ny;
    const int nz = params->nz;
    const double normsq = std::sqrt(xnorma);

    auto *__restrict__ p_cpsi = m_data.cpsi.get_data_restrict();

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                const int idx = i * ny * nz + j * nz + k;
                p_cpsi[idx] /= normsq;
            }
        }
    }
}

void CpuITGrossPitaevskiSolver::imag_iter_linear_step() {
    const int nx     = params->nx;
    const int ny     = params->ny;
    const int nz     = params->nz;
    const double dx  = params->dx;
    const double dy  = params->dy;
    const double dz  = params->dz;
    const double m   = params->m;
    const double dt  = params->time_dt;
    const double cdd = params->cdd;

    const auto *__restrict__ p_cpsi = m_data.cpsi.get_data_restrict();
    const auto *__restrict__ p_fi3d = m_data.fi3d.get_data_restrict();
    const auto *__restrict__ p_pote = m_data.pote.get_data_restrict();
    auto *__restrict__ p_cpsii      = m_data.cpsii.get_data_restrict();

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                const int idx    = i * ny * nz + j * nz + k;
                const int idx_ip = (i + 1) * ny * nz + j * nz + k;
                const int idx_im = (i - 1) * ny * nz + j * nz + k;
                const int idx_jp = i * ny * nz + (j + 1) * nz + k;
                const int idx_jm = i * ny * nz + (j - 1) * nz + k;
                const int idx_kp = i * ny * nz + j * nz + k + 1;
                const int idx_km = i * ny * nz + j * nz + k - 1;
                const auto wav   = p_cpsi[idx];

                double v = p_pote[idx];
                std::complex<double> c1 =
                    -0.5 / (m * std::pow(dx, 2)) * (p_cpsi[idx_im] + p_cpsi[idx_ip] - 2. * wav) -
                    0.5 / (m * std::pow(dy, 2)) * (p_cpsi[idx_jm] + p_cpsi[idx_jp] - 2. * wav) -
                    0.5 / (m * std::pow(dz, 2)) * (p_cpsi[idx_km] + p_cpsi[idx_kp] - 2. * wav) +
                    wav * (v + cdd * p_fi3d[idx]);
                p_cpsii[idx] = wav - dt * c1;
            }
        }
    }
}

void CpuITGrossPitaevskiSolver::imag_iter_nonlinear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    wavefunction_t &cpsi            = m_data.cpsi;
    wavefunction_t &cpsii           = m_data.cpsii;
    const double w                  = params->n_atoms;
    const double ggp11              = params->ggp11;
    const double gamma              = params->ggp11;
    const double w_15               = params->ggp11;
    const double dt                 = params->time_dt;
    const auto *__restrict__ p_cpsi = m_data.cpsi.get_data_restrict();
    auto *__restrict__ p_cpsii      = m_data.cpsii.get_data_restrict();

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                const int idx      = i * ny * nz + j * nz + k;
                const auto wav     = p_cpsi[idx];
                const auto density = std::norm(wav);
                p_cpsii[idx]       = p_cpsii[idx] -
                               dt * density * wav * (ggp11 * w + gamma * std::sqrt(density) * w_15);
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
