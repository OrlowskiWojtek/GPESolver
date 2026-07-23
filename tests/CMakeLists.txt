enable_testing()

find_package(GTest REQUIRED)
find_package(CUDAToolkit REQUIRED)

# === GPU Array Test ===
add_executable(gpu_3d_array_test 
    utils/types/gpu_3d_array_test.cpp
)

target_compile_definitions(gpu_3d_array_test PRIVATE USE_CUDA)
target_include_directories(gpu_3d_array_test PRIVATE 
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(gpu_3d_array_test PRIVATE
    GTest::gtest_main
    pthread
    CUDA::cudart
)

add_test(NAME GpuArrayTest COMMAND gpu_3d_array_test)

# === FileManager Test ===
add_executable(file_manager_test 
    filemanager/file_manager_test.cpp
)

# Define BUILD_TESTS for friend class access
target_compile_definitions(file_manager_test PRIVATE BUILD_TESTS)
target_include_directories(file_manager_test PRIVATE 
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(file_manager_test PRIVATE
    GTest::gtest_main
    pthread
    SolverManager   
)

add_test(NAME FileManagerTest COMMAND file_manager_test)
