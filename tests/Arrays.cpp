#include "gtest/gtest.h"
#include "Arrays.h"


// Check array shape is set correctly
TEST(ArrayTest2D, GetShape){

    unsigned int height = 3;
    unsigned int width = 7;

    Array2D test_array(extents[height][width]);

    ASSERT_EQ(height,test_array.shape()[0]);
    ASSERT_EQ(width,test_array.shape()[1]);
    
}

// Run all tests
int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}