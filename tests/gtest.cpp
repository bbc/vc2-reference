#include "gtest/gtest.h"

TEST(CategoryTest,SpecificTest){
    // Trivial test for testing the tests
    ASSERT_EQ(0,0);
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}