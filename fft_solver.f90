module dip_fft_mem
    use iso_c_binding
    implicit none

    logical :: initialized = .false.

    integer :: nx, ny2, nz2

    complex(8), allocatable :: psi_r(:,:,:), psi_k(:,:,:), Vdip_k(:,:,:)
    type(C_PTR) :: plan_fwd, plan_bwd   ! FFTW plans stay here

contains

    subroutine init_fft_memory(n, ny, nz)
        include 'fftw3.f03'
        integer, intent(in) :: n, ny, nz

        if (initialized) return 

        nx  = 2*n  + 1
        ny2 = 2*ny + 1
        nz2 = 2*nz + 1

        allocate(psi_r(nx,ny2,nz2))
        allocate(psi_k(nx,ny2,nz2))
        allocate(Vdip_k(nx,ny2,nz2))

        ! Tworzymy plany FFT tylko raz
        call fftw_plan_with_nthreads(8)
        plan_fwd = fftw_plan_dft_3d(nx,ny2,nz2, psi_r, psi_k, FFTW_FORWARD, FFTW_MEASURE)
        plan_bwd = fftw_plan_dft_3d(nx,ny2,nz2, psi_k, psi_r, FFTW_BACKWARD, FFTW_MEASURE)

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

end module dip_fft_mem
