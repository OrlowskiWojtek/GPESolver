add_library(TimeDependentSolver
            src/params.cpp
            src/solver.cpp
            src/units.cpp)

target_include_directories(TimeDependentSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
