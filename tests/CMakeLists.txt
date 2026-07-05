enable_testing()

find_package(GTest REQUIRED)
find_package(CUDAToolkit REQUIRED)

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
