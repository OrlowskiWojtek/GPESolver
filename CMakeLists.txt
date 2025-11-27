cmake_minimum_required(VERSION 3.07...3.31)

project(GrossPitaevskiSolver VERSION 1.0)

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_17)

add_subdirectory(utils)
add_subdirectory(solver)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(FFTW3 REQUIRED COMPONENTS double threads)

add_executable(${PROJECT_NAME} main.cpp)

message("Found libraries: " ${FFTW3_LIBRARIES})

target_link_libraries(${PROJECT_NAME} GPSolver ${FFTW3_LIBRARIES} compiler_flags)
