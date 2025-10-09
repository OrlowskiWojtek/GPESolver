cmake_minimum_required(VERSION 3.30)

project(GrossPitaevskiSolver VERSION 1.0)

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_17)

add_subdirectory(solver)
add_subdirectory(utils)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME})
