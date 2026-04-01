cmake_minimum_required(VERSION 3.20)

project(GrossPitaevskiSolver VERSION 1.0)

add_library(compiler_flags INTERFACE)

set(CMAKE_CXX_STANDARD 14)

# Dla kompilatorów intela
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel|IntelLLVM")
    # IntelLLVM = icpx, Intel = icpc
    target_compile_options(compiler_flags INTERFACE
        -O3
        -ffast-math
        -axCORE-AVX2,CORE-AVX512
        #-axarch
        #-xHost 
        #-fp-model 
        #fast=2
        #-ipo                          # link-time optimization
    )
else()
    # GCC / Clang
    # TODO: add Wall
    target_compile_options(compiler_flags INTERFACE
        -Wall
        -O3 
        -ffast-math 
        -funroll-loops
    )
endif()

option(USE_CUDA "Build using CUDA" FALSE)

# add CUDA if option is on
if(USE_CUDA)
    message("Trying to compile using CUDA")
    find_package(CUDA REQUIRED)
    #find_package(CUDAToolkit REQUIRED)

    message("Found CUDA toolkit, ver. " ${CUDAToolkit_VERSION})
    set(CMAKE_CUDA_ARCHITECTURES 60 86)  # GTX 1050 = Pascal, compute 6.1, 86 RTX 3050
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    add_definitions(-DUSE_CUDA)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

add_subdirectory(src)
add_subdirectory(lib)

add_executable(${PROJECT_NAME} main.cpp)


target_link_libraries(${PROJECT_NAME}
                        SolverManager 
                        compiler_flags)

if(NOT USE_CUDA)
	find_package(FFTW3 REQUIRED COMPONENTS double threads)
	message("Found libraries: " ${FFTW3_LIBRARIES})

	message("Linking fftw3 libraries")
	target_link_libraries(${PROJECT_NAME}
		${FFTW3_LIBRARIES})
endif()
