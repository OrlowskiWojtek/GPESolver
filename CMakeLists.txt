cmake_minimum_required(VERSION 3.20)

project(GrossPitaevskiSolver VERSION 1.0)

add_library(compiler_flags INTERFACE)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3 -ffast-math -funroll-loops")

target_compile_options(compiler_flags INTERFACE
        -march=native -O3 -ffast-math -funroll-loops
)
# Dla kompilatorów intela
# a więc u mnie nie ma wielkiej różnicy w czasie  ->
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel|IntelLLVM")
    # IntelLLVM = icpx, Intel = icpc
    target_compile_options(compiler_flags INTERFACE
        -O3 -xHost -fp-model fast=2
        -qopt-zmm-usage=high
        #-ipo                          # link-time optimization
        -qopenmp
    )
else()
    # GCC / Clang
    target_compile_options(compiler_flags INTERFACE
        -O3 -march=native -ffast-math -funroll-loops
    )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(FFTW3 REQUIRED COMPONENTS double threads)

add_subdirectory(src)
add_subdirectory(lib)

add_executable(${PROJECT_NAME} main.cpp)

message("Found libraries: " ${FFTW3_LIBRARIES})

target_link_libraries(${PROJECT_NAME} SolverManager ${FFTW3_LIBRARIES} compiler_flags)
