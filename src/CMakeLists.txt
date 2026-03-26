add_subdirectory(utils)

add_library(SolverManager
    filemanager/file_manager.cpp
    initializer/initializer.cpp
    parameters/parameters.cpp
    manager/sim_manager.cpp
    context/context.cpp
    solver/fft_context.cpp
    solver/fft_poisson_solver.cpp
    solver/fft_rt_split_solver.cpp
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
        solver/cuda_kernels.cu
    )
    
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fast-math -O3")
    
    target_link_libraries(SolverManager PUBLIC
        CUDA::cudart
        CUDA::cufft
    )

    #target_link_libraries(SolverManager PUBLIC SolverKernels)
endif()
