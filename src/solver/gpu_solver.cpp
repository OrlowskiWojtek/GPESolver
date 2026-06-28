#include "solver/gpu_solver.hpp"
#include "solver/cuda_kernels.hpp"

#include <iostream>

GpuGrossPitaevskiSolver::GpuGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator) {
    poisson_solver  = std::make_unique<CUFFTPoissonSolver>();
    rt_split_solver = std::make_unique<CUFFTRealTimeSplitSolver>();
}

void GpuGrossPitaevskiSolver::init_containers() {
    int nx = params->nx;
    int ny = params->ny;
    int nz = params->nz;

    fi3d_gpu = GpuArray<double>(nx, ny, nz);
    cudaMemset(fi3d_gpu.data(), 0, sizeof(double) * fi3d_gpu.size());
}

void GpuGrossPitaevskiSolver::initialize() {
    init_containers();

    //TODO: fix dependency on first allocating memory with copy-operator as this changes adress, look on load_buffer function
    poisson_solver->prepare_gpu(&cpsi_gpu, &fi3d_gpu, &pote_gpu);
    rt_split_solver->prepare_gpu(&cpsi_gpu, &fi3d_gpu, &pote_gpu);

    calc_norm();
    normalize();
}

void GpuGrossPitaevskiSolver::load_buffer(const wavefunction_t &wav) {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    if(wav.size() != nx * ny * nz)
        throw std::runtime_error("bad initialization");
    
    std::cerr << "loading buffer with size: " << wav.size() << std::endl;
    cpsi_gpu  = GpuArray<cuDoubleComplex>(wav.size());
    cpsii_gpu = GpuArray<cuDoubleComplex>(wav.size());

    cpsi.resize(params->nx, params->ny, params->nz);

    auto wav_copy = wav;
    cpsi_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(wav_copy.get_data()));
    cpsii_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(wav_copy.get_data()));
}

void GpuGrossPitaevskiSolver::load_pote(const potential_t &pote) {
    std::cerr << "loading pote" << std::endl;
    this->pote_gpu = GpuArray<double>(pote.size());

    auto pote_copy = pote;
    this->pote_gpu.copy_from_host(pote_copy.get_data());
}

void GpuGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver->execute();
}

void GpuGrossPitaevskiSolver::calc_norm() {
    int N = params->nx * params->ny * params->nz;
    xnorma = launch_kernel_calc_norm(cpsi_gpu.data(), N) * params->get_dxdydz();

    // int nx = params->nx;
    // int ny = params->ny;
    // int nz = params->nz;

    // xnorma = 0.0;
    // for (int i = 0; i < nx; i++) {
    //     for (int j = 0; j < ny; j++) {
    //         for (int k = 0; k < nz; k++) {
    //             xnorma += std::norm(cpsi(i, j, k));
    //         }
    //     }
    // }

    // xnorma *= params->get_dxdydz();
}

void GpuGrossPitaevskiSolver::normalize() {
    int N = params->nx * params->ny * params->nz;
    launch_kernel_normalize(cpsi_gpu.data(), N, xnorma);

    // int nx = params->nx;
    // int ny = params->ny;
    // int nz = params->nz;

    // for (int i = 0; i < nx; i++) {
    //     for (int j = 0; j < ny; j++) {
    //         for (int k = 0; k < nz; k++) {
    //             cpsi(i, j, k) /= std::sqrt(xnorma);
    //         }
    //     }
    // }
}

void GpuGrossPitaevskiSolver::imag_iter_linear_step() {
    launch_kernel_imag_time_iteration(cpsi_gpu.data(),
                                      pote_gpu.data(),
                                      fi3d_gpu.data(),
                                      cpsii_gpu.data(),
                                      params->nx,
                                      params->ny,
                                      params->nz,
                                      params->dx,
                                      params->dy,
                                      params->dz,
                                      params->m,
                                      params->imag_time_dt,
                                      params->cdd,
                                      params->ggp11,
                                      params->gamma,
                                      params->n_atoms,
                                      params->w_15);

    // int nx = params->nx;
    // int ny = params->ny;
    // int nz = params->nz;

    // for (int i = 1; i < nx - 1; i++) {
    //     for (int j = 1; j < ny - 1; j++) {
    //         for (int k = 1; k < nz - 1; k++) {
    //             double v = pote(i, j, k);
    //             std::complex<double> c1 =
    //                 -0.5 / (params->m * std::pow(params->dx, 2)) *
    //                     (cpsi(i - 1, j, k) + cpsi(i + 1, j, k) - 2. * cpsi(i, j, k)) -
    //                 0.5 / (params->m * std::pow(params->dy, 2)) *
    //                     (cpsi(i, j - 1, k) + cpsi(i, j + 1, k) - 2. * cpsi(i, j, k)) -
    //                 0.5 / (params->m * std::pow(params->dz, 2)) *
    //                     (cpsi(i, j, k - 1) + cpsi(i, j, k + 1) - 2. * cpsi(i, j, k)) +
    //                 cpsi(i, j, k) * (v + params->cdd * fi3d(i, j, k));
    //             cpsii(i, j, k) = cpsi(i, j, k) - params->imag_time_dt * c1;
    //         }
    //     }
    // }
}

void GpuGrossPitaevskiSolver::imag_iter_nonlinear_step() {
    // int nx = params->nx;
    // int ny = params->ny;
    // int nz = params->nz;

    // double w = params->n_atoms;
    // for (int i = 1; i < nx - 1; i++) {
    //     for (int j = 1; j < ny - 1; j++) {
    //         for (int k = 1; k < nz - 1; k++) {
    //             cpsii(i, j, k) =
    //                 cpsii(i, j, k) -
    //                 params->imag_time_dt * ((params->ggp11 - params->cdd / 3) *
    //                                             std::norm(cpsi(i, j, k)) * cpsi(i, j, k) * w +
    //                                         params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3)
    //                                         *
    //                                             cpsi(i, j, k) * params->w_15);
    //         }
    //     }
    // }

    cpsi_gpu.swap(cpsii_gpu);
    // cpsi = cpsii;
}

void GpuGrossPitaevskiSolver::real_fft_potential_half_step() {
    launch_kernel_potential_half_step(
        cpsi_gpu.data(),
        pote_gpu.data(),
        fi3d_gpu.data(),
        params->real_time_dt,
        params->cdd,
        params->ggp11,
        params->gamma,
        params->n_atoms,
        params->w_15,
        params->nx,
        params->ny,
        params->nz
    );

    // int nx           = params->nx;
    // int ny           = params->ny;
    // int nz           = params->nz;
    // double w         = params->n_atoms;
    // double dt_factor = params->real_time_dt / 2.;

    // double psi_re;
    // double psi_im;
    // double s;
    // double c;

    // for (int i = 1; i < nx - 1; i++) {
    //     for (int j = 1; j < ny - 1; j++) {
    //         for (int k = 1; k < nz - 1; k++) {
    //             double v_ext   = pote(i, j, k);
    //             double density = std::norm(cpsi(i, j, k));

    //            double v_int = (params->ggp11 - params->cdd / 3) * density * w +
    //                           params->gamma * std::pow(std::abs(cpsi(i, j, k)), 3) *
    //                           params->w_15;

    //            double total_potential = v_ext + params->cdd * fi3d(i, j, k) + v_int;

    //            psi_re = cpsi(i, j, k).real();
    //            psi_im = cpsi(i, j, k).imag();
    //            sincos(-dt_factor * total_potential, &s, &c);

    //            cpsi(i, j, k).imag(psi_re * s + psi_im * c);
    //            cpsi(i, j, k).real(psi_re * c - psi_im * s);
    //        }
    //    }
    //}
}

void GpuGrossPitaevskiSolver::real_fft_kinetic_step() {
    rt_split_solver->execute();
}

void GpuGrossPitaevskiSolver::calc_energy() {
    
    ene.e_kin = 0.;
    ene.e_pot = 0.;
    ene.e_int = 0.;
    ene.e_ext = 0.;
    ene.e_bmf = 0.;

    //int nx = params->nx;
    //int ny = params->ny;
    //int nz = params->nz;

    //std::complex<double> grad_psi_x;
    //std::complex<double> grad_psi_y;
    //std::complex<double> grad_psi_z;

    //for (int i = 1; i < nx - 1; i++) {
    //    for (int j = 1; j < ny - 1; j++) {
    //        for (int k = 1; k < nz - 1; k++) {
    //            // Kinetic energy
    //            grad_psi_x = -(cpsi(i + 1, j, k) + cpsi(i - 1, j, k) - 2. * cpsi(i, j, k)) /
    //                         (std::pow(params->dx, 2));
    //            grad_psi_y = -(cpsi(i, j + 1, k) + cpsi(i, j - 1, k) - 2. * cpsi(i, j, k)) /
    //                         (std::pow(params->dy, 2));
    //            grad_psi_z = -(cpsi(i, j, k + 1) + cpsi(i, j, k - 1) - 2. * cpsi(i, j, k)) /
    //                         (std::pow(params->dz, 2));
    //            ene.e_kin +=
    //                ((grad_psi_x + grad_psi_y + grad_psi_z) * std::conj(cpsi(i, j, k))).real();

    //            // Potential energy
    //            ene.e_pot += pote(i, j, k) * std::norm(cpsi(i, j, k));

    //            // Interaction energy
    //            ene.e_int +=
    //                0.5 * params->ggp11 * std::norm(cpsi(i, j, k)) * std::norm(cpsi(i, j, k));

    //            // Dipole-dipole interaction energy
    //            ene.e_ext +=
    //                0.5 * params->cdd * fi3d(i, j, k) * std::norm(cpsi(i, j, k)) * params->n_atoms;
    //            ene.e_ext -= params->cdd / 3. * std::pow(std::norm(cpsi(i, j, k)), 2) / 2. *
    //                         std::pow(params->n_atoms, 2);

    //            // beyond mean-field energy
    //            ene.e_bmf += 2. / 5. * params->gamma * std::pow(std::abs(cpsi(i, j, k)), 5);
    //        }
    //    }
    //}

    // ene.e_kin *= params->get_dxdydz() / (2 * params->m) * params->n_atoms;
    // ene.e_pot *= params->get_dxdydz() * params->n_atoms;
    // ene.e_int *= params->get_dxdydz() * std::pow(params->n_atoms, 2);
    // ene.e_ext *= params->get_dxdydz();
    // ene.e_bmf *= params->get_dxdydz() * std::pow(params->n_atoms, 2.5);
    // ene.sum();

    launch_kernel_calc_energies(
        cpsi_gpu.data(),
        pote_gpu.data(),
        fi3d_gpu.data(),
        ene,
        params->nx,
        params->ny,
        params->nz,
        params->dx,
        params->dy,
        params->dz,
        params->m,
        params->ggp11,
        params->cdd,
        params->gamma,
        params->n_atoms,
        params->w_15
    );

    enes.emplace_back(ene);
}

void GpuGrossPitaevskiSolver::load_psi(){
    cpsi_gpu.copy_to_host(reinterpret_cast<cuDoubleComplex *>(cpsi.get_data()));
}
