cmake_minimum_required(VERSION 3.10)
project(GPESolver LANGUAGES Fortran)

#include(cmake/FindFFTW.cmake)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# Znajdź bibliotekę FFTW
# Jeśli masz ją w niestandardowym miejscu, podaj: -DFFTW_ROOT=/ścieżka
find_package(FFTW3 REQUIRED COMPONENTS double)

enable_language(Fortran)

# Dodaj pliki źródłowe
add_executable(GPESolver
    main.f90
)

message("Found libraries:" ${FFTW3_LIBRARIES})

# Linkowanie z FFTW
target_include_directories(GPESolver PUBLIC ${FFTW3_INCLUDE_DIR})
target_link_libraries(GPESolver PRIVATE ${FFTW3_LIBRARIES})
