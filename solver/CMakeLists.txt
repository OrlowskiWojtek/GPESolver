add_library(GPSolver
            src/params.cpp
            src/poisson_solver.cpp
            src/solver.cpp)

target_compile_options(GPSolver PRIVATE "-O3;-Wall")
find_package(OpenMP REQUIRED CXX)

target_link_libraries(GPSolver UnitsConverter Types OpenMP::OpenMP_CXX)
target_include_directories(GPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
