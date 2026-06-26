#ifdef USE_CUDA

#include <gtest/gtest.h>
#include "utils/types/mat3d/gpu_3d_array.hpp"

class GpuArrayTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {
        // Cleanup after each test
    }
};

TEST_F(GpuArrayTest, DefaultConstructor) {
    GpuArray<double> arr;

    EXPECT_EQ(arr.size(), 0);
    EXPECT_EQ(arr.nx(), 0);
    EXPECT_EQ(arr.ny(), 0);
    EXPECT_EQ(arr.nz(), 0);
    EXPECT_TRUE(arr.empty());
}

TEST_F(GpuArrayTest, DimsConstructor) {
    GpuArray<double> arr(2, 3, 4);
    EXPECT_EQ(arr.nx(), 2);
    EXPECT_EQ(arr.ny(), 3);
    EXPECT_EQ(arr.nz(), 4);
    EXPECT_EQ(arr.size(), 2 * 3 * 4);
}

TEST_F(GpuArrayTest, CopyFromHost) {
    GpuArray<int> arr(3);
    std::vector<int> host = {1, 2, 3};
    arr.copy_from_host(host);
    
    std::vector<int> result = arr.copy_to_host();
    EXPECT_EQ(result, host);
}

TEST_F(GpuArrayTest, MoveConstructor) {
    GpuArray<double> arr1(10);
    double* original_data = arr1.data();
    
    GpuArray<double> arr2(std::move(arr1));
    
    EXPECT_EQ(arr2.size(), 10);
    EXPECT_EQ(arr2.data(), original_data);
    EXPECT_TRUE(arr1.empty());
}

TEST_F(GpuArrayTest, Reset) {
    GpuArray<double> arr(5);
    ASSERT_FALSE(arr.empty());
    
    arr.reset();
    
    EXPECT_TRUE(arr.empty());
    EXPECT_EQ(arr.size(), 0);
}

TEST_F(GpuArrayTest, Index3D) {
    // nx=2, ny=3, nz=4 -> size=24
    GpuArray<double> arr(2, 3, 4);
    
    // Sprawdź kilka punktów
    EXPECT_EQ(arr.idx(0, 0, 0), 0);
    EXPECT_EQ(arr.idx(0, 0, 1), 1);
    EXPECT_EQ(arr.idx(0, 1, 0), 4);   // i=0, j=1, k=0 -> (0*3+1)*4+0 = 4
    EXPECT_EQ(arr.idx(1, 0, 0), 12);  // i=1, j=0, k=0 -> (1*3+0)*4+0 = 12
}

TEST_F(GpuArrayTest, CopyVector) {
    GpuArray<double> arr(2, 2, 1);
    std::vector<double> host = {1.0, 2.0, 3.0, 4.0};
    
    arr.copy_from_host(host);
    
    std::vector<double> result = arr.copy_to_host();
    EXPECT_EQ(result.size(), 4);
    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[3], 4.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif
