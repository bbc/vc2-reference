#include "gtest/gtest.h"
#include "Quantisation.h"

TEST(QuantisationTest, QuantFactorTestZero){
    
    ASSERT_EQ(quant(12,0),12);
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}