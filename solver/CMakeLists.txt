add_library(TimeDependentSolver
            src/params.cpp
            src/solver.cpp)


target_link_libraries(TimeDependentSolver UnitsConverter)
target_include_directories(TimeDependentSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
