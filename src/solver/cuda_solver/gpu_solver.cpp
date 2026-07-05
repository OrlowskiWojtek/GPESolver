#include "solver/cuda_solver/gpu_solver.hpp"
#include "solver/cuda_solver/cuda_kernels.hpp"

GpuGrossPitaevskiSolver::GpuGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator) {
    poisson_solver  = std::make_unique<CUFFTPoissonSolver>();
    rt_split_solver = std::make_unique<CUFFTRealTimeSplitSolver>();
}

void GpuGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data.allocate(nx, ny, nz);

    cudaMalloc(&d_norm, sizeof(double));
}

void GpuGrossPitaevskiSolver::prepare_fft() {
    poisson_solver->prepare_gpu(&m_data.cpsi_gpu, &m_data.fi3d_gpu, &m_data.pote_gpu);
    rt_split_solver->prepare_gpu(&m_data.cpsi_gpu, &m_data.fi3d_gpu, &m_data.pote_gpu);
}

// void GpuGrossPitaevskiSolver::load_buffer(const wavefunction_t &wav) {
//     size_t nx = params->nx;
//     size_t ny = params->ny;
//     size_t nz = params->nz;
//
//     if(wav.size() != nx * ny * nz)
//         throw std::runtime_error("bad initialization");
//
//     std::cerr << "loading buffer with size: " << wav.size() << std::endl;
//     cpsi_gpu  = GpuArray<cuDoubleComplex>(wav.size());
//     cpsii_gpu = GpuArray<cuDoubleComplex>(wav.size());
//
//     cpsi.resize(params->nx, params->ny, params->nz);
//
//     auto wav_copy = wav;
//     cpsi_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(wav_copy.get_data()));
//     cpsii_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(wav_copy.get_data()));
// }

// void GpuGrossPitaevskiSolver::load_pote(const potential_t &pote) {
//     std::cerr << "loading pote" << std::endl;
//     this->pote_gpu = GpuArray<double>(pote.size());
//
//     auto pote_copy = pote;
//     this->pote_gpu.copy_from_host(pote_copy.get_data());
// }

void GpuGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver->execute();
}

void GpuGrossPitaevskiSolver::calc_norm() {
    int N  = params->nx * params->ny * params->nz;
    xnorma = launch_kernel_calc_norm(m_data.cpsi_gpu.data(), d_norm, N) * params->get_dxdydz();
}

void GpuGrossPitaevskiSolver::normalize() {
    int N = params->nx * params->ny * params->nz;
    launch_kernel_normalize(m_data.cpsi_gpu.data(), N, xnorma);
}

void GpuGrossPitaevskiSolver::imag_iter_linear_step() {
    launch_kernel_imag_time_iteration(m_data.cpsi_gpu.data(),
                                      m_data.pote_gpu.data(),
                                      m_data.fi3d_gpu.data(),
                                      m_data.cpsii_gpu.data(),
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
}

void GpuGrossPitaevskiSolver::imag_iter_nonlinear_step() {
    m_data.cpsi_gpu.swap(m_data.cpsii_gpu);
}

void GpuGrossPitaevskiSolver::real_fft_potential_half_step() {
    launch_kernel_potential_half_step(m_data.cpsi_gpu.data(),
                                      m_data.pote_gpu.data(),
                                      m_data.fi3d_gpu.data(),
                                      params->real_time_dt,
                                      params->cdd,
                                      params->ggp11,
                                      params->gamma,
                                      params->n_atoms,
                                      params->w_15,
                                      params->nx,
                                      params->ny,
                                      params->nz);
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

    launch_kernel_calc_energies(m_data.cpsi_gpu.data(),
                                m_data.pote_gpu.data(),
                                m_data.fi3d_gpu.data(),
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
                                params->w_15);

    ene.e_kin *= params->get_dxdydz() * params->n_atoms;
    ene.e_pot *= params->get_dxdydz() * params->n_atoms;
    ene.e_int *= params->get_dxdydz() * std::pow(params->n_atoms, 2);
    ene.e_ext *= params->get_dxdydz();
    ene.e_bmf *= params->get_dxdydz() * std::pow(params->n_atoms, 2.5);
    ene.sum();

    enes.emplace_back(ene);
}

void GpuGrossPitaevskiSolver::export_data() {
    cudaDeviceSynchronize();
    m_data.cpsi_gpu.copy_to_host(reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
}

void GpuGrossPitaevskiSolver::import_data() {
    if (m_data.cpsi_gpu.size() != buf_data->cpsi.size() ||
        m_data.cpsii_gpu.size() != buf_data->cpsii.size())
        throw std::runtime_error("bad wavefunction import");

    m_data.cpsi_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
    m_data.cpsii_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
}

void GpuGrossPitaevskiSolver::import_pote() {
    if (m_data.pote_gpu.size() != buf_data->pote.size())
        throw std::runtime_error("bad potential import");

    m_data.pote_gpu.copy_from_host(buf_data->pote.get_data());
}
