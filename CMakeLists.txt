cmake_minimum_required(VERSION 3.20)

project(GrossPitaevskiSolver VERSION 1.0)

add_library(compiler_flags INTERFACE)

set(CMAKE_CXX_STANDARD 14)

option(USE_CUDA "Build using CUDA" FALSE)
option(USE_OPENMP "Build using openMP for CPU parallelization" FALSE)
option(BUILD_TESTS "Build unit tests" FALSE)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel|IntelLLVM")
    message("Compiling with intel")
    target_compile_options(compiler_flags INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-O3>
        $<$<COMPILE_LANGUAGE:CXX>:-axCORE-AVX2,CORE-AVX512>
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI|NVIDIA")
    message("Compiling with NVHPC")
    target_compile_options(compiler_flags INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-O3>
    )
else()
    message("Compiling with default")
    target_compile_options(compiler_flags INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall>
        $<$<COMPILE_LANGUAGE:CXX>:-O3>
        $<$<COMPILE_LANGUAGE:CXX>:-march=native>
    )
endif()

# add CUDA if option is on
if(USE_CUDA)
    message("Trying to compile using CUDA")
    find_package(CUDA REQUIRED)
    #find_package(CUDAToolkit REQUIRED)

    message("Found CUDA toolkit, ver. " ${CUDAToolkit_VERSION})
    set(CMAKE_CUDA_ARCHITECTURES 90)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    add_definitions(-DUSE_CUDA)
endif()

if(USE_OPENMP)
    add_definitions(-DUSE_OPENMP)
endif()

if(BUILD_TESTS)
    add_definitions(-DBUILD_TESTS)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

add_subdirectory(src)
add_subdirectory(lib)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
                        SolverManager 
                        compiler_flags)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()
