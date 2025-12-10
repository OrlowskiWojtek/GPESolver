cmake_minimum_required(VERSION 3.07...3.31)

project(GrossPitaevskiSolver VERSION 1.0)

add_library(compiler_flags INTERFACE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3 -ffast-math -funroll-loops")

# Opcje dla konkretnych kompilatorów
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX2 -qopt-zmm-usage=high -ipo")
endif()

add_subdirectory(utils)
add_subdirectory(solver)
add_subdirectory(lib)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(FFTW3 REQUIRED COMPONENTS double threads)

add_executable(${PROJECT_NAME} main.cpp)

message("Found libraries: " ${FFTW3_LIBRARIES})

target_link_libraries(${PROJECT_NAME} GPSolver ${FFTW3_LIBRARIES} compiler_flags)
