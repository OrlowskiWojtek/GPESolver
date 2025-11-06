add_library(GPSolver
            src/params.cpp
            src/solver.cpp)

target_link_libraries(GPSolver UnitsConverter)
target_include_directories(GPSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
