#include "solver/cpu_solver/cpu_it_binary_solver.hpp"
#include "units.hpp"
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

BinaryMixtureParameters *BinaryMixtureParameters::instance = nullptr;

void BinaryMixtureParameters::set_to_default() {
    m_a = UnitConverter::mass_Da_to_au(164.0);
    m_b = UnitConverter::mass_Da_to_au(162.0);

    n_atoms_a = 4.0e4;
    n_atoms_b = 4.0e4;

    // chyba chodzi o pole magnetyczne, reszte przepisuje z parametryzacji z fortrana
    // https://arxiv.org/pdf/1803.10676
    const double b     = 26.0;
    const double a164  = 91.0 * (1.0 - 30.9 / (b - 76.9) - 23.6 / (b - 178.8));
    const double a162  = 220.0 * (1.0 - 1.9 / (b - 21.91) - 0.14 / (b - 26.902));
    const double a6264 = 105.0 * (1.0 - 1.0 / (b - 10.8));

    const double reduced_mass = 2.0 * m_a * m_b / (m_a + m_b);

    g11 = 4.0 * M_PI * a164 / m_a;
    g22 = 4.0 * M_PI * a162 / m_b;
    g12 = 4.0 * M_PI * a6264 / reduced_mass;

    cdd11 = 12.0 * M_PI * 131.0 / m_a;
    cdd22 = 12.0 * M_PI * 131.0 / m_b;
    cdd12 = 12.0 * M_PI * 131.0 / reduced_mass;
}

std::array<double, 8> CpuITBinaryGrossPitaevskiSolver::ug = {0.019855071751231884,
                                                             0.101666761293186630,
                                                             0.237233795041835507,
                                                             0.408282678752175098,
                                                             0.591717321247824902,
                                                             0.762766204958164493,
                                                             0.898333238706813370,
                                                             0.980144928248768116};
std::array<double, 8> CpuITBinaryGrossPitaevskiSolver::wg = {0.050614268145188129,
                                                             0.111190517226687235,
                                                             0.156853322938943644,
                                                             0.181341891689180991,
                                                             0.181341891689180991,
                                                             0.156853322938943644,
                                                             0.111190517226687235,
                                                             0.050614268145188129};

CpuITBinaryGrossPitaevskiSolver::CpuITBinaryGrossPitaevskiSolver(
    AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator)
    , p_mix(BinaryMixtureParameters::getInstance()) {

    p_mix->set_to_default();
}

void CpuITBinaryGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data_a.allocate(nx, ny, nz);
    m_data_b.allocate(nx, ny, nz);

    flhy_a.resize(nx, ny, nz);
    flhy_b.resize(nx, ny, nz);
}

void CpuITBinaryGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver_a->execute();
    poisson_solver_b->execute();
}

void CpuITBinaryGrossPitaevskiSolver::calc_norm() {
    const int nx = params->nx;
    const int ny = params->ny;
    const int nz = params->nz;

    const auto *__restrict__ p_cpsi_a = m_data_a.cpsi.get_data_restrict();
    const auto *__restrict__ p_cpsi_b = m_data_b.cpsi.get_data_restrict();

    norm_a = 0.0;
    norm_b = 0.0;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                const int idx = i * ny * nz + j * nz + k;
                norm_a += std::norm(p_cpsi_a[idx]);
                norm_b += std::norm(p_cpsi_b[idx]);
            }
        }
    }

    norm_a *= params->get_dxdydz();
    norm_b *= params->get_dxdydz();
}

void CpuITBinaryGrossPitaevskiSolver::normalize() {
    const int nx          = params->nx;
    const int ny          = params->ny;
    const int nz          = params->nz;
    const double normsq_a = std::sqrt(norm_a);
    const double normsq_b = std::sqrt(norm_b);

    auto *__restrict__ p_cpsi_a = m_data_a.cpsi.get_data_restrict();
    auto *__restrict__ p_cpsi_b = m_data_b.cpsi.get_data_restrict();

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                const int idx = i * ny * nz + j * nz + k;
                p_cpsi_a[idx] /= normsq_a;
                p_cpsi_b[idx] /= normsq_b;
            }
        }
    }
}

void CpuITBinaryGrossPitaevskiSolver::imag_iter_linear_step() {
    const int nx    = params->nx;
    const int ny    = params->ny;
    const int nz    = params->nz;
    const double dx = params->dx;
    const double dy = params->dy;
    const double dz = params->dz;
    const double dt = params->time_dt;

    const double m_a   = p_mix->m_a;
    const double m_b   = p_mix->m_b;
    const double cdd11 = p_mix->cdd11;
    const double cdd12 = p_mix->cdd12;
    const double cdd22 = p_mix->cdd22;

    const auto *__restrict__ p_cpsi_a = m_data_a.cpsi.get_data_restrict();
    const auto *__restrict__ p_fi3d_a = m_data_a.fi3d.get_data_restrict();
    const auto *__restrict__ p_pote_a = m_data_a.pote.get_data_restrict();
    auto *__restrict__ p_cpsii_a      = m_data_a.cpsii.get_data_restrict();

    const auto *__restrict__ p_cpsi_b = m_data_b.cpsi.get_data_restrict();
    const auto *__restrict__ p_fi3d_b = m_data_b.fi3d.get_data_restrict();
    const auto *__restrict__ p_pote_b = m_data_b.pote.get_data_restrict();
    auto *__restrict__ p_cpsii_b      = m_data_b.cpsii.get_data_restrict();

    // TODO: now there is strong assumption fi3d = both from a and b (as mu_a == mu_b)
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

                const auto wav_a = p_cpsi_a[idx];
                const auto wav_b = p_cpsi_b[idx];
                const double v_a = p_pote_a[idx] + cdd11 * p_fi3d_a[idx] + cdd12 * p_fi3d_b[idx];
                const double v_b = p_pote_b[idx] + cdd12 * p_fi3d_a[idx] + cdd22 * p_fi3d_b[idx];

                std::complex<double> c1a =
                    -0.5 / (m_a * dx * dx) * (p_cpsi_a[idx_im] + p_cpsi_a[idx_ip] - 2. * wav_a) -
                    0.5 / (m_a * dy * dy) * (p_cpsi_a[idx_jm] + p_cpsi_a[idx_jp] - 2. * wav_a) -
                    0.5 / (m_a * dz * dz) * (p_cpsi_a[idx_km] + p_cpsi_a[idx_kp] - 2. * wav_a) +
                    wav_a * v_a;
                p_cpsii_a[idx] = wav_a - dt * c1a;

                std::complex<double> c1b =
                    -0.5 / (m_b * dx * dx) * (p_cpsi_b[idx_im] + p_cpsi_b[idx_ip] - 2. * wav_b) -
                    0.5 / (m_b * dy * dy) * (p_cpsi_b[idx_jm] + p_cpsi_b[idx_jp] - 2. * wav_b) -
                    0.5 / (m_b * dz * dz) * (p_cpsi_b[idx_km] + p_cpsi_b[idx_kp] - 2. * wav_b) +
                    wav_a * v_b;
                p_cpsii_b[idx] = wav_b - dt * c1b;
            }
        }
    }
}

void CpuITBinaryGrossPitaevskiSolver::imag_iter_nonlinear_step() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    const double dt                   = params->time_dt;
    const auto *__restrict__ p_cpsi_a = m_data_a.cpsi.get_data_restrict();
    const auto *__restrict__ p_cpsi_b = m_data_b.cpsi.get_data_restrict();
    auto *__restrict__ p_cp_a         = m_data_a.cpsii.get_data_restrict();
    auto *__restrict__ p_cp_b         = m_data_b.cpsii.get_data_restrict();

    const double ggp11 = p_mix->g11; // contact 4πħ²a_{11}/m_a
    const double ggp22 = p_mix->g22; // contact 4πħ²a_{22}/m_b
    const double ggp12 = p_mix->g12; // contact 4πħ²a_{12}/m_{ab}
    const double na    = p_mix->n_atoms_a;
    const double nb    = p_mix->n_atoms_b;

    const auto *__restrict__ p_flhy_a = flhy_a.get_data_restrict();
    const auto *__restrict__ p_flhy_b = flhy_b.get_data_restrict();

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++) {
                const int idx = i * ny * nz + j * nz + k;

                const auto wa = p_cpsi_a[idx];
                const auto wb = p_cpsi_b[idx];

                const double density_a = std::norm(wa);
                const double density_b = std::norm(wb);

                // same-species + cross-species contact (czynniki N wchodzą w g_ij)
                const double V_a_contact = ggp11 * na * density_a + ggp12 * nb * density_b;
                const double V_b_contact = ggp22 * nb * density_b + ggp12 * na * density_a;

                const double V_a_lhy = p_flhy_a[idx];
                const double V_b_lhy = p_flhy_b[idx];

                p_cp_a[idx] -= dt * wa * (V_a_contact + V_a_lhy);
                p_cp_b[idx] -= dt * wb * (V_b_contact + V_b_lhy);
            }
        }
    }

    // \todo data swap (only pointers)
    m_data_a.cpsi = m_data_a.cpsii;
    m_data_b.cpsi = m_data_b.cpsii;
}

void CpuITBinaryGrossPitaevskiSolver::calc_energy() {
}

void CpuITBinaryGrossPitaevskiSolver::prepare_fft() {
    // ok, for now, normalization is broken, however assuming, that all number of atoms = 40000 (
    poisson_solver_a = std::make_unique<FFTWPoissonSolver>(&m_data_a.cpsi, &m_data_a.fi3d);
    poisson_solver_b = std::make_unique<FFTWPoissonSolver>(&m_data_b.cpsi, &m_data_b.fi3d);
};

void CpuITBinaryGrossPitaevskiSolver::import_pote() {
    m_data_a.pote = buf_data->pote;
    m_data_b.pote = buf_data->pote;
};

void CpuITBinaryGrossPitaevskiSolver::import_data() {
    static int wavefunction_load_count = 0;

    const bool load_to_a = (wavefunction_load_count % 2 == 0);
    auto &target         = load_to_a ? m_data_a : m_data_b;

    target.cpsi  = buf_data->cpsi;
    target.cpsii = buf_data->cpsi;

    wavefunction_load_count++;

    std::cout << std::boolalpha << load_to_a << std::endl;
};

//! \todo Task -> how to export cpsi_b?
void CpuITBinaryGrossPitaevskiSolver::export_data() {
    static int export_a_counter = 0;

    const bool export_to_a = (export_a_counter % 2 == 0);

    auto &target   = export_to_a ? m_data_a : m_data_b;
    buf_data->cpsi = target.cpsi;

    export_a_counter++;
};

const int CpuITBinaryGrossPitaevskiSolver::iter_per_summary() const {
    return 100;
}

void CpuITBinaryGrossPitaevskiSolver::iterate() {
    calc_fi3d();
    calc_lhy();
    imag_iter_linear_step();
    imag_iter_nonlinear_step();
    calc_norm();
    normalize();
}

// Do not need to adjust these calculations
void CpuITBinaryGrossPitaevskiSolver::adjust(int /*iter*/) {
}

void CpuITBinaryGrossPitaevskiSolver::finish() {
    export_data();
    p_mediator->save_initial_state(buf_data->cpsi);
}

void CpuITBinaryGrossPitaevskiSolver::calc_lhy() {
    const int nx = params->nx;
    const int ny = params->ny;
    const int nz = params->nz;

    const auto *pa = m_data_a.cpsi.get_data_restrict();
    const auto *pb = m_data_b.cpsi.get_data_restrict();

    auto *fa = flhy_a.get_data_restrict();
    auto *fb = flhy_b.get_data_restrict();

    const auto na = p_mix->n_atoms_a;
    const auto nb = p_mix->n_atoms_b;

    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int k = 1; k < nz - 1; ++k) {
                const int idx = i * ny * nz + j * nz + k;
                double mu1, mu2;
                lhy_point(std::norm(pa[idx]) * na, std::norm(pb[idx]) * nb, mu1, mu2);
                fa[idx] = mu1;
                fb[idx] = mu2;
            }
        }
    }
}

void CpuITBinaryGrossPitaevskiSolver::lhy_point(double n1, double n2, double &mu1, double &mu2) {
    // Easily precomputable
    constexpr double four_third_pi2 = 4.0 / (3.0 * M_PI * M_PI);
    const double factor_m1          = four_third_pi2 * std::pow(p_mix->m_a, 1.5);
    const double factor_m2          = four_third_pi2 * std::pow(p_mix->m_b, 1.5);

    mu1 = mu2 = 0.0;
    for (int k = 0; k < 8; ++k) {
        const double u       = ug[k];
        const double angular = 3.0 * u * u - 1.0;

        const double g11 = p_mix->g11 + p_mix->cdd11 * angular / 3.0;
        const double g22 = p_mix->g22 + p_mix->cdd22 * angular / 3.0;
        const double g12 = p_mix->g12 + p_mix->cdd12 * angular / 3.0;

        const double diff     = g11 * n1 - g22 * n2;
        const double D        = std::sqrt(diff * diff + 4.0 * g12 * g12 * n1 * n2);
        const double lambda_p = 0.5 * (g11 * n1 + g22 * n2 + D);
        const double lambda_m = 0.5 * (g11 * n1 + g22 * n2 - D);
        const double s_p      = std::pow(lambda_p, 1.5);
        const double s_m      = std::pow(std::max(lambda_m, 0.0), 1.5); // ignore |.|^{3/2} branch

        const double dep_p_n1 =
            0.5 * (g11 + (g11 * diff + 2.0 * g12 * g12 * n2) / std::max(D, 1e-30));
        const double dep_m_n1 =
            0.5 * (g11 - (g11 * diff + 2.0 * g12 * g12 * n2) / std::max(D, 1e-30));
        const double dep_p_n2 =
            0.5 * (g22 + (-g22 * diff + 2.0 * g12 * g12 * n1) / std::max(D, 1e-30));
        const double dep_m_n2 =
            0.5 * (g22 - (-g22 * diff + 2.0 * g12 * g12 * n1) / std::max(D, 1e-30));

        // mu1 += wg[k] * factor_m1 * (s_p * dep_p_n1 + s_m * dep_m_n1);
        // mu2 += wg[k] * factor_m2 * (s_p * dep_p_n2 + s_m * dep_m_n2);
        mu1 += wg[k] * (factor_m1 * s_m * dep_m_n1 + factor_m2 * s_p * dep_p_n1);
        mu2 += wg[k] * (factor_m2 * s_m * dep_m_n2 + factor_m1 * s_p * dep_p_n2);
    }
}
