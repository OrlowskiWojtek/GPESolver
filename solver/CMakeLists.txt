add_library(GPSolver
            src/params.cpp
            src/solver.cpp)

target_compile_options(GPSolver PRIVATE "-O3;-Wall")

target_link_libraries(GPSolver UnitsConverter Types)
target_include_directories(GPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
