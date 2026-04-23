add_subdirectory(utils)

add_library(SolverManager
    filemanager/file_manager.cpp
    initializer/initializer.cpp
    parameters/parameters.cpp
    manager/sim_manager.cpp
    context/context.cpp
    solver/fft_solvers/fft_context.cpp
    solver/solver.cpp
)

target_link_libraries(SolverManager PUBLIC
    Output
    Types
    Json
    UnitsConverter
    compiler_flags        
)

target_include_directories(SolverManager PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    
    set(CMAKE_CUDA_ARCHITECTURES 60)
    target_link_libraries(SolverManager PUBLIC CUDA::cufftw CUDA::cufft CUDA::cudart)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")
endif()

if(USE_CUDA)
    enable_language(CUDA)

    target_sources(SolverManager PRIVATE
        solver/fft_solvers/gpu/cuda_fft_kernels.cu
        solver/fft_solvers/gpu/cufft_poisson_solver.cpp
        solver/fft_solvers/gpu/cufft_split_solver.cpp
    )
    
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fast-math -O3")
    
    target_link_libraries(SolverManager PUBLIC
        CUDA::cudart
        CUDA::cufft
    )
else()
	find_package(FFTW3 REQUIRED COMPONENTS double threads)
	message("Found libraries: " ${FFTW3_LIBRARIES})

	message("Linking fftw3 libraries")

    target_sources(SolverManager PRIVATE
        solver/fft_solvers/cpu/fftw_poisson_solver.cpp
        solver/fft_solvers/cpu/fftw_split_solver.cpp
    )

    target_link_libraries(SolverManager PUBLIC
		${FFTW3_LIBRARIES})
endif()
