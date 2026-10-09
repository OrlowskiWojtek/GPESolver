#include "solver/cuda_solver/gpu_it_binary_solver.hpp"
#include "solver/cuda_solver/cuda_binary_kernels.cuh"
#include "solver/cuda_solver/cuda_kernels.cuh"
#include "utils/output/output.hpp"
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

GpuITBinaryGrossPitaevskiSolver::GpuITBinaryGrossPitaevskiSolver(
    AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator)
    , p_mix(BinaryMixtureParameters::getInstance()) {

    cudaDeviceEnablePeerAccess(0, 0);

    cudaSetDevice(0);
    cudaStreamCreateWithFlags(&copy_stream_0, cudaStreamNonBlocking);
    cudaMalloc(&d_norm_a, sizeof(double));

    cudaSetDevice(1);
    cudaStreamCreateWithFlags(&copy_stream_1, cudaStreamNonBlocking);
    cudaMalloc(&d_norm_b, sizeof(double));

    cudaSetDevice(0);

    if (!initialize_gauss()) {
        throw std::runtime_error("Can't copy static arrays onto GPU");
    }

    p_mix->set_to_default();
}

GpuITBinaryGrossPitaevskiSolver::~GpuITBinaryGrossPitaevskiSolver() {
    cudaSetDevice(0);
    cudaStreamDestroy(copy_stream_0);
    cudaFree(d_norm_a);

    cudaSetDevice(1);
    cudaStreamDestroy(copy_stream_1);
    cudaFree(d_norm_b);

    cudaSetDevice(0);
}

void GpuITBinaryGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    cudaSetDevice(0);
    m_data_a_0.allocate(nx, ny, nz);
    m_data_b_0.allocate(nx, ny, nz);

    cudaSetDevice(1);
    m_data_a_1.allocate(nx, ny, nz);
    m_data_b_1.allocate(nx, ny, nz);

    cudaSetDevice(0);
    flhy_a = GpuArray<double>(nx, ny, nz);
    flhy_b = GpuArray<double>(nx, ny, nz);
}

void GpuITBinaryGrossPitaevskiSolver::calc_fi3d() {
    cudaSetDevice(0); // dont know if this is needed
    poisson_solver_a
        ->execute(); // this uses dens_a, output pote_a, works ok gpu 0

    cudaSetDevice(1);
    poisson_solver_b
        ->execute(); // this uses dens_b, output pote_b, works on gpu 1

    cudaSetDevice(0);
}

void GpuITBinaryGrossPitaevskiSolver::calc_norm() {
    const int N = params->nx * params->ny * params->nz;

    cudaSetDevice(0);
    launch_kernel_calc_norm(
        m_data_a_0.cpsi_gpu.data(), d_norm_a, params->get_dxdydz(), N);

    cudaSetDevice(1);
    launch_kernel_calc_norm(
        m_data_b_1.cpsi_gpu.data(), d_norm_b, params->get_dxdydz(), N);

    cudaSetDevice(0);
}

void GpuITBinaryGrossPitaevskiSolver::normalize() {
    const int N = params->nx * params->ny * params->nz;

    cudaSetDevice(0);
    launch_kernel_normalize(m_data_a_0.cpsi_gpu.data(), d_norm_a, N);

    cudaSetDevice(1);
    launch_kernel_normalize(m_data_b_1.cpsi_gpu.data(), d_norm_b, N);

    cudaSetDevice(0);
}

void GpuITBinaryGrossPitaevskiSolver::calc_energy() {
}

void GpuITBinaryGrossPitaevskiSolver::prepare_fft() {
    // ok, for now, normalization is broken, however assuming, that all number
    // of atoms = 40000 (
    cudaSetDevice(0);
    poisson_solver_a = std::make_unique<CUFFTPoissonSolver>(
        &m_data_a_0.cpsi_gpu, &m_data_a_0.fi3d_gpu);
    cudaSetDevice(1);
    poisson_solver_b = std::make_unique<CUFFTPoissonSolver>(
        &m_data_b_1.cpsi_gpu, &m_data_b_1.fi3d_gpu);
    cudaSetDevice(0);
};

void GpuITBinaryGrossPitaevskiSolver::import_pote() {
    if (m_data_a_0.pote_gpu.size() != buf_data->pote.size())
        throw std::runtime_error("bad potential import");

    m_data_a_0.pote_gpu.copy_from_host(buf_data->pote.get_data());
    m_data_a_1.pote_gpu.copy_from_host(buf_data->pote.get_data());
};

void GpuITBinaryGrossPitaevskiSolver::import_data() {
    static int wavefunction_load_count = 0;

    const bool load_to_a = (wavefunction_load_count % 2 == 0);
    auto &target_0       = load_to_a ? m_data_a_0 : m_data_b_0;
    auto &target_1       = load_to_a ? m_data_a_1 : m_data_b_1;

    if (target_0.cpsi_gpu.size() != buf_data->cpsi.size() ||
        target_0.cpsii_gpu.size() != buf_data->cpsii.size())
        throw std::runtime_error("bad wavefunction import");

    target_0.cpsi_gpu.copy_from_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
    target_0.cpsii_gpu.copy_from_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));

    target_1.cpsi_gpu.copy_from_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
    target_1.cpsii_gpu.copy_from_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));

    wavefunction_load_count++;

    OutputFormatter::printInfo(std::string("Loading to wavefunction: ") +
                               std::string(load_to_a ? "a" : "b"));
    std::cout << std::boolalpha << load_to_a << std::endl;
};

//! \todo Task -> how to export cpsi_b?
void GpuITBinaryGrossPitaevskiSolver::export_data() {
    static int export_a_counter = 0;

    const bool export_to_a = (export_a_counter % 2 == 0);

    auto &target = export_to_a ? m_data_a_0 : m_data_b_1;

    cudaDeviceSynchronize();
    target.cpsi_gpu.copy_to_host(
        reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));

    export_a_counter++;
};

const int GpuITBinaryGrossPitaevskiSolver::iter_per_summary() const {
    return 100;
}

void GpuITBinaryGrossPitaevskiSolver::iterate() {
    copy_data_from_gpu1();
    MEASURE_NVTX(calc_fi3d);

    cudaSetDevice(0);
    cudaStreamSynchronize(copy_stream_0);

    copy_fi3d_from_gpu1();
    MEASURE_NVTX(calc_lhy); // this uses dens_a and dens_b, now from both gpus
    //
    
    cudaStreamSynchronize(copy_stream_0);

    MEASURE_NVTX(imag_iter_full_step);

    copy_data_from_gpu0();
    MEASURE_NVTX(calc_norm);
    MEASURE_NVTX(normalize);
}

void GpuITBinaryGrossPitaevskiSolver::imag_iter_full_step() {
    launch_kernel_imag_time_iteration(m_data_a_0.cpsi_gpu.data(),
                                      m_data_b_0.cpsi_gpu.data(),
                                      m_data_a_0.pote_gpu.data(),
                                      m_data_a_0.fi3d_gpu.data(),
                                      m_data_b_0.fi3d_gpu.data(),
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
    launch_kernel_calc_lhy(m_data_a_0.cpsi_gpu.data(),
                           m_data_b_0.cpsi_gpu.data(),
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
                           params->nz,
                           p_mix->four_third_pi2,
                           p_mix->factor_m1,
                           p_mix->factor_m2);
}

void GpuITBinaryGrossPitaevskiSolver::copy_fi3d_from_gpu1() {
    cudaSetDevice(1);

    const size_t fi3d_bytes = m_data_b_1.fi3d_gpu.size() * sizeof(double);

    cudaError_t error = cudaMemcpyPeerAsync(m_data_b_0.fi3d_gpu.data(),
                                            0,
                                            m_data_b_1.fi3d_gpu.data(),
                                            1,
                                            fi3d_bytes,
                                            copy_stream_0);

    if (error != cudaSuccess) {
        throw std::runtime_error(
            std::string("Failed to copy b.fi3d from GPU 1 to GPU 0: ") +
            cudaGetErrorString(error));
    }

    cudaSetDevice(0);
}

void GpuITBinaryGrossPitaevskiSolver::copy_data_from_gpu1() {
    cudaSetDevice(1);

    const size_t cpsi_bytes =
        m_data_b_1.cpsi_gpu.size() * sizeof(cuDoubleComplex);

    cudaError_t error = cudaMemcpyPeerAsync(m_data_b_0.cpsi_gpu.data(),
                                            0,
                                            m_data_b_1.cpsi_gpu.data(),
                                            1,
                                            cpsi_bytes,
                                            copy_stream_0);

    if (error != cudaSuccess) {
        throw std::runtime_error(
            std::string("Failed to copy b.cpsi from GPU 1 to GPU 0: ") +
            cudaGetErrorString(error));
    }

    cudaSetDevice(0);
}

void GpuITBinaryGrossPitaevskiSolver::copy_data_from_gpu0() {
    cudaSetDevice(0);

    const size_t cpsi_bytes =
        m_data_b_0.cpsi_gpu.size() * sizeof(cuDoubleComplex);

    cudaError_t error = cudaMemcpyPeer(m_data_b_1.cpsi_gpu.data(),
                                       1,
                                       m_data_b_0.cpsi_gpu.data(),
                                       0,
                                       cpsi_bytes);

    if (error != cudaSuccess) {
        throw std::runtime_error(
            std::string("Failed to copy b.cpsi from GPU 1 to GPU 0: ") +
            cudaGetErrorString(error));
    }
}
