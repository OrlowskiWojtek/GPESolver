add_library(GPSolver
            src/params.cpp
            src/file_manager.cpp
            src/fft_context.cpp
            src/fft_poisson_solver.cpp
            src/fft_rt_split_solver.cpp
            src/solver.cpp)

target_link_libraries(GPSolver UnitsConverter Types json)
target_include_directories(GPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
