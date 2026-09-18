#include "solver/cpu_solver/cpu_it_split_solver.hpp"
#ifdef USE_OPENMP
#include <omp.h>
#endif

CpuITSplitGrossPitaevskiSolver::CpuITSplitGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator) {
}

void CpuITSplitGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data.allocate(nx, ny, nz);
}

void CpuITSplitGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver->execute();
}

void CpuITSplitGrossPitaevskiSolver::calc_norm() {
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

void CpuITSplitGrossPitaevskiSolver::normalize() {
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

void CpuITSplitGrossPitaevskiSolver::calc_energy() {
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

void CpuITSplitGrossPitaevskiSolver::prepare_fft() {
    poisson_solver = std::make_unique<FFTWPoissonSolver>(&m_data.cpsi, &m_data.fi3d);
    rt_imag_split_solver = std::make_unique<FFTWImagTimeSplitSolver>(&m_data.cpsi, &m_data.fi3d);
};

void CpuITSplitGrossPitaevskiSolver::import_pote() {
    m_data.pote = buf_data->pote;
};

void CpuITSplitGrossPitaevskiSolver::import_data() {
    m_data.cpsi  = buf_data->cpsi;
    m_data.cpsii = buf_data->cpsii;
};

void CpuITSplitGrossPitaevskiSolver::export_data() {
    buf_data->cpsi = m_data.cpsi;
};

const int CpuITSplitGrossPitaevskiSolver::iter_per_summary() const {
    return 100;
}

// Do not need to adjust these calculations
void CpuITSplitGrossPitaevskiSolver::adjust(int /*iter*/) {
}

void CpuITSplitGrossPitaevskiSolver::finish() {
    export_data();
    p_mediator->save_initial_state(buf_data->cpsi);
}

void CpuITSplitGrossPitaevskiSolver::iterate() {
    imag_fft_potential_half_step();
    calc_fi3d();
    imag_fft_kinetic_step();
    imag_fft_potential_half_step();

    calc_norm();
    normalize();
}

void CpuITSplitGrossPitaevskiSolver::imag_fft_potential_half_step() {
    const int nx           = params->nx;
    const int ny           = params->ny;
    const int nz           = params->nz;
    const double w         = params->n_atoms;
    const double dt_factor = params->time_dt / 2.;
    const double ggp11     = params->ggp11;
    const double gamma     = params->gamma;
    const double w15       = params->w_15;
    const double cdd       = params->cdd;

    wavefunction_t &cpsi = m_data.cpsi;

#ifdef USE_OPENMP
    #pragma omp parallel for simd collapse(3)
#endif
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                double v_ext   = m_data.pote(i, j, k);
                double density = std::norm(cpsi(i, j, k));

                double v_int = ggp11 * density * w +
                               gamma * density * std::sqrt(density) * w15;

                double total_potential = v_ext + cdd * m_data.fi3d(i, j, k) + v_int;
                cpsi(i,j,k) *= std::exp( - dt_factor * total_potential);
            }
        }
    }
}

void CpuITSplitGrossPitaevskiSolver::imag_fft_kinetic_step(){
    rt_imag_split_solver->execute();
}
