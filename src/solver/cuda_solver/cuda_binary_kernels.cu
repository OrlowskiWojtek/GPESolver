#include "solver/cuda_solver/cuda_binary_kernels.cuh"
#include <stdio.h>


__constant__ double c_ug[8];
__constant__ double c_wg[8];

bool initialize_gauss(){
    constexpr double ug[8] = {
        0.019855071751231884,
        0.101666761293186630,
        0.237233795041835507,
        0.408282678752175098,
        0.591717321247824902,
        0.762766204958164493,
        0.898333238706813370,
        0.980144928248768116,
    };

    constexpr double wg[8] = {
        0.050614268145188129,
        0.111190517226687235,
        0.156853322938943644,
        0.181341891689180991,
        0.181341891689180991,
        0.156853322938943644,
        0.111190517226687235,
        0.050614268145188129,
    };

    cudaError_t error = cudaMemcpyToSymbol(
        c_ug,
        ug,
        sizeof(ug),
        0,
        cudaMemcpyHostToDevice);

    if (error != cudaSuccess) {
        return false;
    }

    error = cudaMemcpyToSymbol(
        c_wg,
        wg,
        sizeof(wg),
        0,
        cudaMemcpyHostToDevice);

    if (error != cudaSuccess) {
        return false;
    }

    return true;
}

__global__
void kernel_imag_time_iteration(
    cuDoubleComplex *__restrict__ cpsi_a,
    cuDoubleComplex *__restrict__ cpsi_b,
    const double *__restrict__ pote,
    const double *__restrict__ fi3d_a,
    const double *__restrict__ fi3d_b,
    const double *__restrict__ flhy_a,
    const double *__restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double n_atoms_a_15,
    const double n_atoms_b_15,
    const double cdd_11,
    const double cdd_12,
    const double cdd_22,
    const double ggp_11,
    const double ggp_12,
    const double ggp_22,
    const int nx,
    const int ny,
    const int nz,
    const double dx,
    const double dy,
    const double dz,
    const double dt
){
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int i = blockIdx.z * blockDim.z + threadIdx.z;

    // (warunki brzegowe Dirichleta)
    if (i < 1 || i >= nx - 1 || j < 1 || j >= ny - 1 || k < 1 || k >= nz - 1) {
        return;
    }

    const int idx    = i * ny * nz + j * nz + k;
    const int idx_ip1 = (i + 1) * ny * nz + j * nz + k;
    const int idx_im1 = (i - 1) * ny * nz + j * nz + k;
    const int idx_jp1 = i * ny * nz + (j + 1) * nz + k;
    const int idx_jm1 = i * ny * nz + (j - 1) * nz + k;
    const int idx_kp1 = i * ny * nz + j * nz + k + 1;
    const int idx_km1 = i * ny * nz + j * nz + k - 1;

    const cuDoubleComplex psi_a     = __ldg(&cpsi_a[idx]);
    const cuDoubleComplex psi_im1_a = __ldg(&cpsi_a[idx_im1]);
    const cuDoubleComplex psi_ip1_a = __ldg(&cpsi_a[idx_ip1]);
    const cuDoubleComplex psi_jm1_a = __ldg(&cpsi_a[idx_jm1]);
    const cuDoubleComplex psi_jp1_a = __ldg(&cpsi_a[idx_jp1]);
    const cuDoubleComplex psi_km1_a = __ldg(&cpsi_a[idx_km1]);
    const cuDoubleComplex psi_kp1_a = __ldg(&cpsi_a[idx_kp1]);

    const cuDoubleComplex psi_b     = __ldg(&cpsi_b[idx]);
    const cuDoubleComplex psi_im1_b = __ldg(&cpsi_b[idx_im1]);
    const cuDoubleComplex psi_ip1_b = __ldg(&cpsi_b[idx_ip1]);
    const cuDoubleComplex psi_jm1_b = __ldg(&cpsi_b[idx_jm1]);
    const cuDoubleComplex psi_jp1_b = __ldg(&cpsi_b[idx_jp1]);
    const cuDoubleComplex psi_km1_b = __ldg(&cpsi_b[idx_km1]);
    const cuDoubleComplex psi_kp1_b = __ldg(&cpsi_b[idx_kp1]);

    const double v_pote = pote[idx];
    const double v_fi3d_a = fi3d_a[idx];
    const double v_fi3d_b = fi3d_b[idx];

    const double v_a = v_pote + cdd_11 * v_fi3d_a +
        cdd_12 * v_fi3d_b;
    const double v_b = v_pote + cdd_12 * v_fi3d_b +
        cdd_22 * v_fi3d_a;

    double coef_x = -0.5 / (m_a * dx * dx);
    double coef_y = -0.5 / (m_a * dy * dy);
    double coef_z = -0.5 / (m_a * dz * dz);

    // ==== part a ====
    cuDoubleComplex laplacian_a;
    laplacian_a.x = coef_x * (psi_im1_a.x + psi_ip1_a.x - 2.0 * psi_a.x) +
                    coef_y * (psi_jm1_a.x + psi_jp1_a.x - 2.0 * psi_a.x) +
                    coef_z * (psi_km1_a.x + psi_kp1_a.x - 2.0 * psi_a.x);
    laplacian_a.y = coef_x * (psi_im1_a.y + psi_ip1_a.y - 2.0 * psi_a.y) +
                    coef_y * (psi_jm1_a.y + psi_jp1_a.y - 2.0 * psi_a.y) +
                    coef_z * (psi_km1_a.y + psi_kp1_a.y - 2.0 * psi_a.y);

    coef_x = -0.5 / (m_b * dx * dx);
    coef_y = -0.5 / (m_b * dy * dy);
    coef_z = -0.5 / (m_b * dz * dz);
    // ==== part b ====
    cuDoubleComplex laplacian_b;
    laplacian_b.x = coef_x * (psi_im1_b.x + psi_ip1_b.x - 2.0 * psi_b.x) +
                    coef_y * (psi_jm1_b.x + psi_jp1_b.x - 2.0 * psi_b.x) +
                    coef_z * (psi_km1_b.x + psi_kp1_b.x - 2.0 * psi_b.x);
    laplacian_b.y = coef_x * (psi_im1_b.y + psi_ip1_b.y - 2.0 * psi_b.y) +
                    coef_y * (psi_jm1_b.y + psi_jp1_b.y - 2.0 * psi_b.y) +
                    coef_z * (psi_km1_b.y + psi_kp1_b.y - 2.0 * psi_b.y);

    cuDoubleComplex linear_a;
    linear_a.x = laplacian_a.x + (v_a) * psi_a.x;
    linear_a.y = laplacian_a.y + (v_a) * psi_a.y;

    cuDoubleComplex linear_b;
    linear_b.x = laplacian_b.x + (v_b) * psi_b.x;
    linear_b.y = laplacian_b.y + (v_b) * psi_b.y;

    // ======== nonlinear part ============

    const double density_a = psi_a.x * psi_a.x + psi_a.y * psi_a.y;
    const double density_b = psi_b.x * psi_b.x + psi_b.y * psi_b.y;

    const double V_a_contact =
        ggp_11 * n_atoms_a * density_a + ggp_12 * n_atoms_b * density_b;
    const double V_b_contact =
        ggp_22 * n_atoms_b * density_b + ggp_12 * n_atoms_a * density_a;
    const double V_a_lhy = flhy_a[idx];
    const double V_b_lhy = flhy_b[idx];
    
    cuDoubleComplex nonlinear_a;
    nonlinear_a.x = psi_a.x * (V_a_contact + V_a_lhy);
    nonlinear_a.y = psi_a.y * (V_a_contact + V_a_lhy);
    cuDoubleComplex nonlinear_b;
    nonlinear_b.x = psi_b.x * (V_b_contact + V_b_lhy);
    nonlinear_b.y = psi_b.y * (V_b_contact + V_b_lhy);

    // ============ FINALIZATION ==============
    cpsi_a[idx].x = psi_a.x - dt * (linear_a.x + nonlinear_a.x);
    cpsi_a[idx].y = psi_a.y - dt * (linear_a.y + nonlinear_a.y);

    cpsi_b[idx].x = psi_b.x - dt * (linear_b.x + nonlinear_b.x);
    cpsi_b[idx].y = psi_b.y - dt * (linear_b.y + nonlinear_b.y);
}

void launch_kernel_imag_time_iteration(
    cuDoubleComplex *__restrict__ cpsi_a,
    cuDoubleComplex *__restrict__ cpsi_b,
    const double *__restrict__ pote,
    const double *__restrict__ fi3d_a,
    const double *__restrict__ fi3d_b,
    const double *__restrict__ flhy_a,
    const double *__restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double n_atoms_a_15,
    const double n_atoms_b_15,
    const double cdd_11,
    const double cdd_12,
    const double cdd_22,
    const double ggp_11,
    const double ggp_12,
    const double ggp_22,
    const int nx,
    const int ny,
    const int nz,
    const double dx,
    const double dy,
    const double dz,
    const double dt
) {
    dim3 block(8, 8, 8);
    dim3 grid((nz + 7) / 8, (ny + 7) / 8, (nx + 7) / 8);

    kernel_imag_time_iteration<<<grid, block>>>(
        cpsi_a, cpsi_b, pote, fi3d_a, fi3d_b, flhy_a, flhy_b,
        m_a, m_b, n_atoms_a, n_atoms_b, n_atoms_a_15, n_atoms_b_15,
        cdd_11, cdd_12, cdd_22, ggp_11, ggp_12, ggp_22,
        nx, ny, nz, dx, dy, dz, dt
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after kernel_imag_time_iteration: %s\n", cudaGetErrorString(err));
    }
}

__global__
void kernel_calc_lhy(
    const cuDoubleComplex* __restrict__ psi_a,
    const cuDoubleComplex* __restrict__ psi_b,
    double* __restrict__ flhy_a,
    double* __restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double g11,
    const double g12,
    const double g22,
    const double cdd11,
    const double cdd12,
    const double cdd22,
    const int nx,
    const int ny,
    const int nz,
    const double four_third_pi2,
    const double factor_m1,
    const double factor_m2
) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int i = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < 1 || i >= nx - 1 || j < 1 || j >= ny - 1 || k < 1 || k >= nz - 1) {
        return;
    }

    const int idx = i * ny * nz + j * nz + k;

    // Densities: |psi|^2 * n_atoms
    const double n1 = (psi_a[idx].x * psi_a[idx].x + psi_a[idx].y * psi_a[idx].y) * n_atoms_a;
    const double n2 = (psi_b[idx].x * psi_b[idx].x + psi_b[idx].y * psi_b[idx].y) * n_atoms_b;

    double mu1 = 0.0;
    double mu2 = 0.0;

    #pragma unroll
    for (int q = 0; q < 8; ++q) {
        const double u       = c_ug[q];
        const double wg     = c_wg[q];
        const double angular = 3.0 * u * u - 1.0;

        const double g11_eff = g11 + cdd11 * angular / 3.0;
        const double g22_eff = g22 + cdd22 * angular / 3.0;
        const double g12_eff = g12 + cdd12 * angular / 3.0;

        const double diff = g11_eff * n1 - g22_eff * n2;
        const double D    = sqrt(fmax(diff * diff + 4.0 * g12_eff * g12_eff * n1 * n2, 1e-30));
        const double lambda_p = 0.5 * (g11_eff * n1 + g22_eff * n2 + D);
        const double lambda_m = 0.5 * (g11_eff * n1 + g22_eff * n2 - D);

        // change to val * sqrt(val)
        const double s_p      = pow(fmax(lambda_p, 0.0), 1.5);
        const double s_m      = pow(fmax(lambda_m, 0.0), 1.5);

        const double dep_p_n1 = 0.5 * (g11_eff + (g11_eff * diff + 2.0 * g12_eff * g12_eff * n2) / D);
        const double dep_m_n1 = 0.5 * (g11_eff - (g11_eff * diff + 2.0 * g12_eff * g12_eff * n2) / D);
        const double dep_p_n2 = 0.5 * (g22_eff + (-g22_eff * diff + 2.0 * g12_eff * g12_eff * n1) / D);
        const double dep_m_n2 = 0.5 * (g22_eff - (-g22_eff * diff + 2.0 * g12_eff * g12_eff * n1) / D);

        mu1 += wg * (factor_m1 * s_m * dep_m_n1 + factor_m2 * s_p * dep_p_n1);
        mu2 += wg * (factor_m2 * s_m * dep_m_n2 + factor_m1 * s_p * dep_p_n2);
    }

    flhy_a[idx] = mu1;
    flhy_b[idx] = mu2;
}

void launch_kernel_calc_lhy(
    const cuDoubleComplex* __restrict__ psi_a,
    const cuDoubleComplex* __restrict__ psi_b,
    double* __restrict__ flhy_a,
    double* __restrict__ flhy_b,
    const double m_a,
    const double m_b,
    const double n_atoms_a,
    const double n_atoms_b,
    const double g11,
    const double g12,
    const double g22,
    const double cdd11,
    const double cdd12,
    const double cdd22,
    const int nx,
    const int ny,
    const int nz,
    const double four_third_pi2,
    const double factor_m1,
    const double factor_m2
) {
    dim3 block(8, 8, 8);
    dim3 grid((nz + 7) / 8, (ny + 7) / 8, (nx + 7) / 8);

    kernel_calc_lhy<<<grid, block>>>(
        psi_a, psi_b, flhy_a, flhy_b,
        m_a, m_b, n_atoms_a, n_atoms_b,
        g11, g12, g22, cdd11, cdd12, cdd22,
        nx, ny, nz,
        four_third_pi2, factor_m1, factor_m2
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after calc_lhy kernel: %s\n", cudaGetErrorString(err));
    }
}
