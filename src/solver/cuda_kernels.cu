#include "solver/cuda_kernels.hpp"

__global__
void imag_time_iteration_kernel(
    const cuDoubleComplex *__restrict__ cpsi, 
    const double *__restrict__ pote,          
    const double *__restrict__ fi3d,          
    cuDoubleComplex *__restrict__ cpsii,      
    int nx,
    int ny,
    int nz,
    double dx,
    double dy,
    double dz,
    double m,
    double imag_dt,
    double cdd,
    double ggp11,
    double gamma,
    double w,   // n_atoms
    double w_15 // w_15
) {
    // Indeks wątku w przestrzeni 3D
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // Pomijamy brzegi (warunki brzegowe Dirichlet)
    if (i < 1 || i >= nx - 1 || j < 1 || j >= ny - 1 || k < 1 || k >= nz - 1) {
        return;
    }

    // Linear idx (row-major: i*ny*nz + j*nz + k)
    int idx     = (i * ny + j) * nz + k;
    int idx_im1 = ((i - 1) * ny + j) * nz + k;
    int idx_ip1 = ((i + 1) * ny + j) * nz + k;
    int idx_jm1 = (i * ny + (j - 1)) * nz + k;
    int idx_jp1 = (i * ny + (j + 1)) * nz + k;
    int idx_km1 = (i * ny + j) * nz + (k - 1);
    int idx_kp1 = (i * ny + j) * nz + (k + 1);

    // Odczyt z __ldg() dla danych tylko do odczytu
    cuDoubleComplex psi     = __ldg(&cpsi[idx]);
    cuDoubleComplex psi_im1 = __ldg(&cpsi[idx_im1]);
    cuDoubleComplex psi_ip1 = __ldg(&cpsi[idx_ip1]);
    cuDoubleComplex psi_jm1 = __ldg(&cpsi[idx_jm1]);
    cuDoubleComplex psi_jp1 = __ldg(&cpsi[idx_jp1]);
    cuDoubleComplex psi_km1 = __ldg(&cpsi[idx_km1]);
    cuDoubleComplex psi_kp1 = __ldg(&cpsi[idx_kp1]);

    double v  = __ldg(&pote[idx]);
    double fi = __ldg(&fi3d[idx]);

    // Współczynniki dla Laplace'a
    double coef_x = -0.5 / (m * dx * dx);
    double coef_y = -0.5 / (m * dy * dy);
    double coef_z = -0.5 / (m * dz * dz);

    // Laplacian = psi(i+1) + psi(i-1) - 2*psi + ... (6 sąsiadów)
    cuDoubleComplex laplacian;
    laplacian.x = coef_x * (psi_im1.x + psi_ip1.x - 2.0 * psi.x) +
                  coef_y * (psi_jm1.x + psi_jp1.x - 2.0 * psi.x) +
                  coef_z * (psi_km1.x + psi_kp1.x - 2.0 * psi.x);
    laplacian.y = coef_x * (psi_im1.y + psi_ip1.y - 2.0 * psi.y) +
                  coef_y * (psi_jm1.y + psi_jp1.y - 2.0 * psi.y) +
                  coef_z * (psi_km1.y + psi_kp1.y - 2.0 * psi.y);

    // |psi|^2 = psi.x^2 + psi.y^2
    double norm_psi = psi.x * psi.x + psi.y * psi.y;
    double abs_psi  = sqrt(norm_psi); // lub cuCabs(psi)

    // Część liniowa: L(psi) = -1/(2m)*∇²psi + (v + cdd*fi3d)*psi
    cuDoubleComplex linear;
    linear.x = laplacian.x + (v + cdd * fi) * psi.x;
    linear.y = laplacian.y + (v + cdd * fi) * psi.y;

    // Część nieliniowa: N(psi) = [(g - cdd/3)*|psi|² + γ*|psi|³] * psi * w
    double nonlinear_factor = (ggp11 - cdd / 3.0) * norm_psi * w + gamma * pow(abs_psi, 3.0) * w_15;

    cuDoubleComplex nonlinear;
    nonlinear.x = nonlinear_factor * psi.x;
    nonlinear.y = nonlinear_factor * psi.y;

    // Suma i aktualizacja: psi_nowa = psi - dt * (linear + nonlinear)
    double dt    = imag_dt;
    cpsii[idx].x = psi.x - dt * (linear.x + nonlinear.x);
    cpsii[idx].y = psi.y - dt * (linear.y + nonlinear.y);
}

void launch_kernel_imag_time_iteration(
    const cuDoubleComplex* d_cpsi,
    const double* d_pote,
    const double* d_fi3d,
    cuDoubleComplex* d_cpsii,
    int nx, int ny, int nz,
    double dx, double dy, double dz,
    double m,
    double imag_dt,
    double cdd,
    double ggp11,
    double gamma,
    double w,
    double w_15)
{
    dim3 block(8, 8, 8);
    dim3 grid((nx + 7) / 8, (ny + 7) / 8, (nz + 7) / 8);

    imag_time_iteration_kernel<<<grid, block>>>(
        d_cpsi, d_pote, d_fi3d, d_cpsii,
        nx, ny, nz,
        dx, dy, dz, m, imag_dt, cdd, ggp11, gamma, w, w_15
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
}

__global__ 
void kernel_normalize(
    cuDoubleComplex* __restrict__ data,
    int N,
    double norm_factor
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        data[idx].x /= sqrt(norm_factor);
        data[idx].y /= sqrt(norm_factor);
    }
}

__global__ void kernel_calc_norm(
    const cuDoubleComplex* __restrict__ data,
    double* __restrict__ result,
    int N
) {
    __shared__ double sdata[256];
    
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    double sum = 0.0;
    if (idx < N) {
        sum = data[idx].x * data[idx].x + data[idx].y * data[idx].y;
    }
    sdata[tid] = sum;
    __syncthreads();
    
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s && idx + s < N) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        atomicAdd(result, sdata[0]);
    }
}

double launch_kernel_calc_norm(
    const cuDoubleComplex* data,
    int N
) {
    double h_result = 0.0;
    double* d_result;
    cudaMalloc(&d_result, sizeof(double));
    cudaMemset(d_result, 0, sizeof(double));
    
    int block = 256;
    int grid  = (N + block - 1) / block;
    kernel_calc_norm<<<grid, block>>>(data, d_result, N);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
    
    cudaMemcpy(&h_result, d_result, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_result);
    
    return h_result;
}

void launch_kernel_normalize(
    cuDoubleComplex* data,
    int N,
    double norm_factor
) {
    int block = 256;
    int grid  = (N + block - 1) / block;
    kernel_normalize<<<grid, block>>>(data, N, norm_factor);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
}

__global__ 
void kernel_calc_energies(
    const cuDoubleComplex* __restrict__ psi,
    const double* __restrict__ pote,
    const double* __restrict__ fi3d,
    double* __restrict__ e_kin,
    double* __restrict__ e_pot,
    double* __restrict__ e_int,
    double* __restrict__ e_ext,
    double* __restrict__ e_bmf,
    int nx, int ny, int nz,
    double dx, double dy, double dz,
    double m, double ggp11, double cdd, double gamma,
    double n_atoms, double w_15
) {
    __shared__ double s_kin[256];
    __shared__ double s_pot[256];
    __shared__ double s_int[256];
    __shared__ double s_ext[256];
    __shared__ double s_bmf[256];
    
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    double kin = 0.0, pot = 0.0, eint = 0.0, ext = 0.0, bmf = 0.0;
    
    if (idx < nx * ny * nz) {
        int i = idx / (ny * nz);
        int j = (idx / nz) % ny;
        int k = idx % nz;
        
        // Boundary check - skip edges for derivatives
        if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1 && k > 0 && k < nz - 1) {
            double psi_re = psi[idx].x;
            double psi_im = psi[idx].y;
            double psi_norm2 = psi_re * psi_re + psi_im * psi_im;
            
            // Kinetic energy
            double d2x = (psi[idx - ny * nz].x + psi[idx + ny * nz].x - 2.0 * psi_re) / (dx * dx);
            double d2y = (psi[idx - nz].x + psi[idx + nz].x - 2.0 * psi_re) / (dy * dy);
            double d2z = (psi[idx - 1].x + psi[idx + 1].x - 2.0 * psi_re) / (dz * dz);
            
            kin = -(d2x + d2y + d2z) * psi_re;
            
            // Potential energy
            pot = pote[idx] * psi_norm2;
            
            // Interaction energy
            eint = 0.5 * ggp11 * psi_norm2 * psi_norm2;
            
            // Dipole-dipole energy
            ext = 0.5 * cdd * fi3d[idx] * psi_norm2 * n_atoms 
                  - cdd / 3.0 * psi_norm2 * psi_norm2 / 2.0 * n_atoms * n_atoms;
            
            // Beyond mean-field
            bmf = 2.0 / 5.0 * gamma * pow(psi_norm2, 2.5);
        }
    }
    
    s_kin[tid] = kin;
    s_pot[tid] = pot;
    s_int[tid] = eint;
    s_ext[tid] = ext;
    s_bmf[tid] = bmf;
    __syncthreads();
    
    // Reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            s_kin[tid] += s_kin[tid + s];
            s_pot[tid] += s_pot[tid + s];
            s_int[tid] += s_int[tid + s];
            s_ext[tid] += s_ext[tid + s];
            s_bmf[tid] += s_bmf[tid + s];
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        atomicAdd(e_kin, s_kin[0]);
        atomicAdd(e_pot, s_pot[0]);
        atomicAdd(e_int, s_int[0]);
        atomicAdd(e_ext, s_ext[0]);
        atomicAdd(e_bmf, s_bmf[0]);
    }
}

void launch_kernel_calc_energies(
    const cuDoubleComplex* psi,
    const double* pote,
    const double* fi3d,
    energies_t& ene,
    int nx, int ny, int nz,
    double dx, double dy, double dz,
    double m, double ggp11, double cdd, double gamma,
    double n_atoms, double w_15
) {
    double *d_kin_dev, *d_pot_dev, *d_int_dev, *d_ext_dev, *d_bmf_dev;
    
    cudaMalloc(&d_kin_dev, sizeof(double));
    cudaMalloc(&d_pot_dev, sizeof(double));
    cudaMalloc(&d_int_dev, sizeof(double));
    cudaMalloc(&d_ext_dev, sizeof(double));
    cudaMalloc(&d_bmf_dev, sizeof(double));
     
    int N = nx * ny * nz;
    int block = 256;
    int grid = (N + block - 1) / block;
    
    kernel_calc_energies<<<grid, block>>>(
        psi, pote, fi3d, d_kin_dev, d_pot_dev, d_int_dev, d_ext_dev, d_bmf_dev,
        nx, ny, nz, dx, dy, dz, m, ggp11, cdd, gamma, n_atoms, w_15
    );
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
    
    cudaMemcpy(&ene.e_kin, d_kin_dev, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ene.e_pot, d_pot_dev, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ene.e_int, d_int_dev, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ene.e_ext, d_ext_dev, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ene.e_bmf, d_bmf_dev, sizeof(double), cudaMemcpyDeviceToHost);
    
    cudaFree(d_kin_dev);

    ene.sum();
}

__global__ 
void kernel_potential_half_step_inplace(
    cuDoubleComplex* __restrict__ psi,
    const double* __restrict__ V,
    const double* __restrict__ fi3d,
    double dt_factor,
    double cdd,
    double ggp11,
    double gamma,
    double n_atoms,
    double w_15,
    int nx, int ny, int nz
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (i >= nx || j >= ny || k >= nz) return;
    
    int idx = (i * ny + j) * nz + k;
    
    double v_ext = __ldg(&V[idx]);
    double fi = __ldg(&fi3d[idx]);
    
    double psi_re = psi[idx].x;
    double psi_im = psi[idx].y;
    double norm_psi = psi_re * psi_re + psi_im * psi_im;
    double abs_psi = sqrt(norm_psi);
    
    double v_int = (ggp11 - cdd / 3.0) * norm_psi * n_atoms 
                   + gamma * pow(abs_psi, 3.0) * w_15;
    
    double total_potential = v_ext + cdd * fi + v_int;
    
    double ph = -dt_factor * total_potential;
    double c = cos(ph);
    double s = sin(ph);
    
    psi[idx].x = psi_re * c - psi_im * s;
    psi[idx].y = psi_re * s + psi_im * c;
}

void launch_kernel_potential_half_step(
    cuDoubleComplex* psi,
    const double* V,
    const double* fi3d,
    double dt,
    double cdd,
    double ggp11,
    double gamma,
    double n_atoms,
    double w_15,
    int nx, int ny, int nz
) {
    dim3 block(8, 8, 8);
    dim3 grid((nx + 7) / 8, (ny + 7) / 8, (nz + 7) / 8);
    
    double dt_factor = dt / 2.0;
    
    kernel_potential_half_step_inplace<<<grid, block>>>(
        psi, V, fi3d, dt_factor, cdd, ggp11, gamma, n_atoms, w_15, nx, ny, nz
    );
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error after imag_Time_Iteration kernel: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
}

