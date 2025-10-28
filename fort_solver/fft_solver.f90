module FFTSolver
    use iso_c_binding
    implicit none

    logical :: initialized = .false.

    integer :: nx2, ny2, nz2
    real(8) :: pi

    complex(8), allocatable :: psi_r(:,:,:), psi_k(:,:,:), Vdip_k(:,:,:)
    real(8), allocatable :: fi3(:,:,:) ! copy for fi3do
    type(C_PTR) :: plan_fwd, plan_bwd   ! FFTW plans stay here

contains

    subroutine init_fft_memory(nx, ny, nz)
        include 'fftw3.f03'
        integer, intent(in) :: nx, ny, nz

        if (initialized) return 

        nx2 = 2*nx  + 1
        ny2 = 2*ny + 1
        nz2 = 2*nz + 1

        pi=4*atan(1.0)

        allocate(psi_r(nx2,ny2,nz2))
        allocate(psi_k(nx2,ny2,nz2))
        allocate(Vdip_k(nx2,ny2,nz2))
        allocate(fi3(-nx:nx, -ny:ny, -nz:nz))

        ! Tworzymy plany FFT tylko raz
        call fftw_plan_with_nthreads(8)
        plan_fwd = fftw_plan_dft_3d(nx2,ny2,nz2, psi_r, psi_k, FFTW_FORWARD, FFTW_MEASURE)
        plan_bwd = fftw_plan_dft_3d(nx2,ny2,nz2, psi_k, psi_r, FFTW_BACKWARD, FFTW_MEASURE)

        initialized = .true.
    end subroutine init_fft_memory

    subroutine free_fft_memory()
        include 'fftw3.f03'
        if (.not. initialized) return
        call fftw_destroy_plan(plan_fwd)
        call fftw_destroy_plan(plan_bwd)
        deallocate(psi_r, psi_k, Vdip_k)
        initialized = .false.
    end subroutine free_fft_memory

end module FFTSolver
