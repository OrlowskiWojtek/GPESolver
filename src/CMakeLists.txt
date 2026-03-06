add_subdirectory(utils)

add_library(SolverManager
    filemanager/file_manager.cpp
    initializer/initializer.cpp
    parameters/parameters.cpp
    manager/sim_manager.cpp
    context/context.cpp
    solver/fft_context.cpp
    solver/fft_poisson_solver.cpp
    solver/fft_rt_split_solver.cpp
    solver/solver.cpp
)

target_link_libraries(SolverManager PUBLIC
    Output
    Types
    Json
    UnitsConverter
    compiler_flags        
)

target_include_directories(SolverManager PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
