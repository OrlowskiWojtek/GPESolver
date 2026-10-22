#include "solver/cuda_solver/gpu_it_binary_solver.hpp"
#include "solver/cuda_solver/cuda_binary_kernels.cuh"
#include "solver/cuda_solver/cuda_kernels.cuh"
#include "utils/output/output.hpp"

GpuITBinaryGrossPitaevskiSolver::GpuITBinaryGrossPitaevskiSolver(
    AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator)
    , p_mix(BinaryMixtureParameters::getInstance()) {

    cudaMalloc(&d_norm_a, sizeof(double));
    cudaMalloc(&d_norm_b, sizeof(double));

    if (!initialize_gauss()) {
        throw std::runtime_error("Can't copy static arrays onto GPU");
    }
}

GpuITBinaryGrossPitaevskiSolver::~GpuITBinaryGrossPitaevskiSolver() {
    cudaFree(d_norm_a);
    cudaFree(d_norm_b);
}

void GpuITBinaryGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data.allocate(nx, ny, nz);

    flhy_a = GpuArray<double>(nx, ny, nz);
    flhy_b = GpuArray<double>(nx, ny, nz);
}

void GpuITBinaryGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver_a->execute();
    poisson_solver_b->execute();
}

void GpuITBinaryGrossPitaevskiSolver::calc_norm() {
    const int N = params->nx * params->ny * params->nz;

    // TODO: move into common kernel of two
    launch_kernel_calc_norm(
        m_data.cpsi_gpu_a.data(), d_norm_a, params->get_dxdydz(), N);
    launch_kernel_calc_norm(
        m_data.cpsi_gpu_b.data(), d_norm_b, params->get_dxdydz(), N);
}

void GpuITBinaryGrossPitaevskiSolver::normalize() {
    const int N = params->nx * params->ny * params->nz;

    // TODO: move into common kernel of two
    launch_kernel_normalize(m_data.cpsi_gpu_a.data(), d_norm_a, N);
    launch_kernel_normalize(m_data.cpsi_gpu_b.data(), d_norm_b, N);
}

void GpuITBinaryGrossPitaevskiSolver::calc_energy() {
    nvtxRangePushA("CalcEnergyMix");

    launch_kernel_calc_enes(m_data.cpsi_gpu_a.data(),
                            m_data.cpsi_gpu_b.data(),
                            m_data.pote_gpu.data(),
                            m_data.fi3d_gpu_a.data(),
                            m_data.fi3d_gpu_b.data(),
                            flhy_a.data(),
                            flhy_b.data(),
                            d_energy_a.get_data(),
                            d_energy_b.get_data(),
                            p_mix->spec_a.mass,
                            p_mix->spec_b.mass,
                            p_mix->spec_a.n_atoms,
                            p_mix->spec_b.n_atoms,
                            p_mix->g11,
                            p_mix->g12,
                            p_mix->g22,
                            p_mix->cdd11,
                            p_mix->cdd12,
                            p_mix->cdd22,
                            params->nx,
                            params->ny,
                            params->nz,
                            params->dx,
                            params->dy,
                            params->dz);

    // TODO: alloc on CPU (?)
    double energies[5];

    for(int i = 0; i < 2; i++){
        auto& target = i ? d_energy_a : d_energy_b;
        cudaError_t err = cudaMemcpy(energies,
                                     target.get_data(),
                                     5 * sizeof(double),
                                     cudaMemcpyDeviceToHost);

        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("Can't copy") +
                                     std::string(cudaGetErrorString(err)));
        }

        ene.e_kin = energies[0];
        ene.e_pot = energies[1];
        ene.e_int = energies[2];
        ene.e_ext = energies[3];
        ene.e_bmf = energies[4];
        ene.sum();

        enes.push_back(ene);
    }

    nvtxRangePop();
}

void GpuITBinaryGrossPitaevskiSolver::prepare_fft() {
    // ok, for now, normalization is broken, however assuming, that all number
    // of atoms = 40000
    poisson_solver_a = std::make_unique<CUFFTPoissonSolver>(&m_data.cpsi_gpu_a,
                                                            &m_data.fi3d_gpu_a);
    poisson_solver_b = std::make_unique<CUFFTPoissonSolver>(&m_data.cpsi_gpu_b,
                                                            &m_data.fi3d_gpu_b);
};

void GpuITBinaryGrossPitaevskiSolver::import_pote() {
    if (m_data.pote_gpu.size() != buf_data->pote.size())
        throw std::runtime_error("bad potential import");

    m_data.pote_gpu.copy_from_host(buf_data->pote.get_data());
};

void GpuITBinaryGrossPitaevskiSolver::import_data() {
}

void GpuITBinaryGrossPitaevskiSolver::import_binary_data() {
    if (m_data.cpsi_gpu_a.size() != binary_buf_data.spec_a.size() ||
        m_data.cpsi_gpu_b.size() != binary_buf_data.spec_b.size()) {
        throw std::runtime_error("Binary wavefunction dimensions do not match solver arrays.");
    }

    m_data.cpsi_gpu_a.copy_from_host(reinterpret_cast<const cuDoubleComplex*>(
        binary_buf_data.spec_a.get_data()));
    m_data.cpsi_gpu_b.copy_from_host(reinterpret_cast<const cuDoubleComplex*>(
        binary_buf_data.spec_b.get_data()));
}

void GpuITBinaryGrossPitaevskiSolver::export_data() {
    export_binary_data();
}

void GpuITBinaryGrossPitaevskiSolver::export_binary_data() {
    m_data.cpsi_gpu_a.copy_to_host(reinterpret_cast<cuDoubleComplex*>(
        binary_buf_data.spec_a.get_data()));
    m_data.cpsi_gpu_b.copy_to_host(reinterpret_cast<cuDoubleComplex*>(
        binary_buf_data.spec_b.get_data()));
}

const int GpuITBinaryGrossPitaevskiSolver::iter_per_summary() const {
    return 100;
}

void GpuITBinaryGrossPitaevskiSolver::iterate() {
    MEASURE_NVTX(calc_fi3d);
    MEASURE_NVTX(calc_lhy);
    MEASURE_NVTX(imag_iter_full_step);
    MEASURE_NVTX(calc_norm);
    MEASURE_NVTX(normalize);
}

void GpuITBinaryGrossPitaevskiSolver::imag_iter_full_step() {
    launch_kernel_imag_time_iteration(m_data.cpsi_gpu_a.data(),
                                      m_data.cpsi_gpu_b.data(),
                                      m_data.pote_gpu.data(),
                                      m_data.fi3d_gpu_a.data(),
                                      m_data.fi3d_gpu_b.data(),
                                      flhy_a.data(),
                                      flhy_b.data(),
                                      p_mix->spec_a.mass,
                                      p_mix->spec_b.mass,
                                      p_mix->spec_a.n_atoms,
                                      p_mix->spec_b.n_atoms,
                                      p_mix->spec_a.n_atoms * 1.5,
                                      p_mix->spec_b.n_atoms * 1.5,
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
    export_binary_data();
    p_mediator->save_binary_initial_state(binary_buf_data);
}

void GpuITBinaryGrossPitaevskiSolver::calc_lhy() {
    launch_kernel_calc_lhy(m_data.cpsi_gpu_a.data(),
                           m_data.cpsi_gpu_b.data(),
                           flhy_a.data(),
                           flhy_b.data(),
                           p_mix->spec_a.mass,
                           p_mix->spec_b.mass,
                           p_mix->spec_a.n_atoms,
                           p_mix->spec_b.n_atoms,
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
