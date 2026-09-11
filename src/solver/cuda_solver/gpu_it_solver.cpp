#include "solver/cuda_solver/gpu_it_solver.hpp"
#include "solver/cuda_solver/cuda_kernels.hpp"

GpuITGrossPitaevskiSolver::GpuITGrossPitaevskiSolver(AbstractSimulationMediator *mediator)
    : AbstractGrossPitaevskiSolver(mediator) {

    cudaMalloc(&d_norm, sizeof(double));

    cudaMalloc(&d_kin_dev, sizeof(double));
    cudaMalloc(&d_pot_dev, sizeof(double));
    cudaMalloc(&d_int_dev, sizeof(double));
    cudaMalloc(&d_ext_dev, sizeof(double));
    cudaMalloc(&d_bmf_dev, sizeof(double));
}

GpuITGrossPitaevskiSolver::~GpuITGrossPitaevskiSolver() {
    cudaFree(d_norm);

    cudaFree(d_kin_dev);
    cudaFree(d_pot_dev);
    cudaFree(d_int_dev);
    cudaFree(d_ext_dev);
    cudaFree(d_bmf_dev);
}

void GpuITGrossPitaevskiSolver::init_containers() {
    size_t nx = params->nx;
    size_t ny = params->ny;
    size_t nz = params->nz;

    m_data.allocate(nx, ny, nz);
}

void GpuITGrossPitaevskiSolver::prepare_fft() {
    poisson_solver = std::make_unique<CUFFTPoissonSolver>(&m_data.cpsi_gpu, &m_data.fi3d_gpu);
}

void GpuITGrossPitaevskiSolver::calc_fi3d() {
    poisson_solver->execute();
}

void GpuITGrossPitaevskiSolver::calc_norm() {
    int N  = params->nx * params->ny * params->nz;
    xnorma = launch_kernel_calc_norm(m_data.cpsi_gpu.data(), d_norm, N) * params->get_dxdydz();
}

void GpuITGrossPitaevskiSolver::normalize() {
    int N = params->nx * params->ny * params->nz;
    launch_kernel_normalize(m_data.cpsi_gpu.data(), N, xnorma);
}

void GpuITGrossPitaevskiSolver::imag_iteration_full() {
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
                                      params->time_dt,
                                      params->cdd,
                                      params->ggp11,
                                      params->gamma,
                                      params->n_atoms,
                                      params->w_15);

    m_data.cpsi_gpu.swap(m_data.cpsii_gpu);
}

void GpuITGrossPitaevskiSolver::calc_energy() {
    ene.e_kin = 0.;
    ene.e_pot = 0.;
    ene.e_int = 0.;
    ene.e_ext = 0.;
    ene.e_bmf = 0.;

    launch_kernel_calc_energies(m_data.cpsi_gpu.data(),
                                m_data.pote_gpu.data(),
                                m_data.fi3d_gpu.data(),
                                ene,
                                d_kin_dev,
                                d_pot_dev,
                                d_int_dev,
                                d_ext_dev,
                                d_bmf_dev,
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

void GpuITGrossPitaevskiSolver::export_data() {
    cudaDeviceSynchronize();
    m_data.cpsi_gpu.copy_to_host(reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
}

void GpuITGrossPitaevskiSolver::import_data() {
    if (m_data.cpsi_gpu.size() != buf_data->cpsi.size() ||
        m_data.cpsii_gpu.size() != buf_data->cpsii.size())
        throw std::runtime_error("bad wavefunction import");

    m_data.cpsi_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
    m_data.cpsii_gpu.copy_from_host(reinterpret_cast<cuDoubleComplex *>(buf_data->cpsi.get_data()));
}

void GpuITGrossPitaevskiSolver::import_pote() {
    if (m_data.pote_gpu.size() != buf_data->pote.size())
        throw std::runtime_error("bad potential import");

    m_data.pote_gpu.copy_from_host(buf_data->pote.get_data());
}

void GpuITGrossPitaevskiSolver::iterate() {
    calc_fi3d();
    imag_iteration_full();
    calc_norm();
    normalize();
}

void GpuITGrossPitaevskiSolver::adjust(int iter) {

}

void GpuITGrossPitaevskiSolver::finish() {
    export_data();
    p_mediator->save_initial_state(buf_data->cpsi);
}

const int GpuITGrossPitaevskiSolver::iter_per_summary() const {
    return 1000;
}
