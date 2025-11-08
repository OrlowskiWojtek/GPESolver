cmake_minimum_required(VERSION 3.30)

project(GrossPitaevskiSolver VERSION 1.0)

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_17)

add_subdirectory(utils)
add_subdirectory(solver)

find_package(Armadillo REQUIRED)

if(Armadillo_FOUND AND NOT TARGET Armadillo::Armadillo)
  add_library(Armadillo::Armadillo INTERFACE IMPORTED)
  set_target_properties(
    Armadillo::Armadillo
    PROPERTIES
      INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}"
  )
endif()

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} Armadillo::Armadillo GPSolver)
