#include "solver/cuda_solver/gpu_it_binary_solver.hpp"
#include "solver/cuda_solver/cuda_binary_kernels.cuh"
#include "solver/cuda_solver/cuda_kernels.cuh"
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

GpuITBinaryGrossPitaevskiSolver::GpuITBinaryGrossPitaevskiSolver(
    AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator)
    , p_mix(BinaryMixtureParameters::getInstance()) {

    cudaMalloc(&d_norm_a, sizeof(double));
    cudaMalloc(&d_norm_b, sizeof(double));

    initialize_gauss();
    p_mix->set_to_default();
}

GpuITBinaryGrossPitaevskiSolver::~GpuITBinaryGrossPitaevskiSolver(){
    cudaFree(&d_norm_a);
    cudaFree(&d_norm_b);
}

void GpuITBinaryGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data_a.allocate(nx, ny, nz);
    m_data_b.allocate(nx, ny, nz);

    flhy_a = GpuArray<double>(nx, ny, nz);
    flhy_b = GpuArray<double>(nx, ny, nz);
}

void GpuITBinaryGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver_a->execute();
    poisson_solver_b->execute();
}

void GpuITBinaryGrossPitaevskiSolver::calc_norm() {
    const int N = params->nx * params->ny * params->nz;
    launch_kernel_calc_norm(
        m_data_a.cpsi_gpu.data(), d_norm_a, params->get_dxdydz(), N);
    launch_kernel_calc_norm(
        m_data_b.cpsi_gpu.data(), d_norm_b, params->get_dxdydz(), N);
}

void GpuITBinaryGrossPitaevskiSolver::normalize() {
    const int N = params->nx * params->ny * params->nz;

    launch_kernel_normalize(m_data_a.cpsi_gpu.data(), d_norm_a, N);
    launch_kernel_normalize(m_data_b.cpsi_gpu.data(), d_norm_b, N);
}

void GpuITBinaryGrossPitaevskiSolver::calc_energy() {
}

void GpuITBinaryGrossPitaevskiSolver::prepare_fft() {
    // ok, for now, normalization is broken, however assuming, that all number
    // of atoms = 40000 (
    //
    poisson_solver_a = std::make_unique<CUFFTPoissonSolver>(&m_data_a.cpsi_gpu,
                                                            &m_data_a.fi3d_gpu);
    poisson_solver_b = std::make_unique<CUFFTPoissonSolver>(&m_data_b.cpsi_gpu,
                                                            &m_data_b.fi3d_gpu);
};

void GpuITBinaryGrossPitaevskiSolver::import_pote() {
    if (m_data_a.pote_gpu.size() != buf_data->pote.size())
        throw std::runtime_error("bad potential import");

    m_data_a.pote_gpu.copy_from_host(buf_data->pote.get_data());
};

void GpuITBinaryGrossPitaevskiSolver::import_data() {
    static int wavefunction_load_count = 0;

    const bool load_to_a = (wavefunction_load_count % 2 == 0);
    auto &target         = load_to_a ? m_data_a : m_data_b;

    if (target.cpsi_gpu.size() != buf_data->cpsi.size() ||
        target.cpsii_gpu.size() != buf_data->cpsii.size())
        throw std::runtime_error("bad wavefunction import");

    std::cerr << "IMPORTING FIRST, size: " << target.cpsi_gpu.size() << std::endl;

    target.cpsi_gpu.copy_from_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
    std::cerr << "IMPORTING FIRST, size: " << target.cpsi_gpu.size() << std::endl;
    target.cpsii_gpu.copy_from_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));

    wavefunction_load_count++;

    std::cout << std::boolalpha << load_to_a << std::endl;
};

//! \todo Task -> how to export cpsi_b?
void GpuITBinaryGrossPitaevskiSolver::export_data() {
    static int export_a_counter = 0;

    const bool export_to_a = (export_a_counter % 2 == 0);

    auto &target = export_to_a ? m_data_a : m_data_b;

    cudaDeviceSynchronize();
    target.cpsi_gpu.copy_to_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));

    export_a_counter++;
};

const int GpuITBinaryGrossPitaevskiSolver::iter_per_summary() const {
    return 100;
}

void GpuITBinaryGrossPitaevskiSolver::iterate() {
    calc_fi3d();
    calc_lhy();
    imag_iter_full_step();
    calc_norm();
    normalize();
}

void GpuITBinaryGrossPitaevskiSolver::imag_iter_full_step() {
    launch_kernel_imag_time_iteration(m_data_a.cpsi_gpu.data(),
                                      m_data_b.cpsi_gpu.data(),
                                      m_data_a.pote_gpu.data(),
                                      m_data_a.fi3d_gpu.data(),
                                      m_data_b.fi3d_gpu.data(),
                                      flhy_a.data(),
                                      flhy_b.data(),
                                      p_mix->m_a,
                                      p_mix->m_b,
                                      p_mix->n_atoms_a,
                                      p_mix->n_atoms_b,
                                      p_mix->n_atoms_a * 1.5,
                                      p_mix->n_atoms_b * 1.5,
                                      p_mix->cdd11,
                                      p_mix->cdd12,
                                      p_mix->cdd22,
                                      p_mix->g11,
                                      p_mix->g12,
                                      p_mix->g22,
                                      params->nx,
                                      params->ny,
                                      params->nz,
                                      params->dx,
                                      params->dy,
                                      params->dz,
                                      params->time_dt);
}

// Do not need to adjust these calculations
void GpuITBinaryGrossPitaevskiSolver::adjust(int /*iter*/) {
}

void GpuITBinaryGrossPitaevskiSolver::finish() {
    export_data();
    p_mediator->save_initial_state(buf_data->cpsi);
}

void GpuITBinaryGrossPitaevskiSolver::calc_lhy() {
    launch_kernel_calc_lhy(m_data_a.cpsi_gpu.data(),
                           m_data_b.cpsi_gpu.data(),
                           flhy_a.data(),
                           flhy_b.data(),
                           p_mix->m_a,
                           p_mix->m_b,
                           p_mix->n_atoms_a,
                           p_mix->n_atoms_b,
                           p_mix->g11,
                           p_mix->g12,
                           p_mix->g22,
                           p_mix->cdd11,
                           p_mix->cdd12,
                           p_mix->cdd22,
                           params->nx,
                           params->ny,
                           params->nz);
}
