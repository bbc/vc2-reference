#include "gtest/gtest.h"
#include "Quantisation.h"


// Check correct exception is thrown when a quantisation factor is too large
TEST(QuantisationFactorTest, ThrowErrorIfTooLarge){
    try{
        quant(12,130);
    } catch (std::logic_error& exception){
        EXPECT_EQ(std::string(exception.what()), "quantization index exceeds maximum implemented value.");
    }
}


// Parameterised test to check quantising a few values
struct quantised_values{
    int value;
    int q;
    int result;
};

struct QuantisationFactorTest : testing::Test, testing::WithParamInterface<quantised_values>{};

TEST_P(QuantisationFactorTest, CheckNormalOperation){
    auto as = GetParam();
    EXPECT_EQ(quant(as.value, as.q), as.result);

};

INSTANTIATE_TEST_SUITE_P(Default, QuantisationFactorTest,
    testing::Values(
        quantised_values{12,0,12},
        quantised_values{12,2,8},
        quantised_values{-12,2,-8},
        quantised_values{-12,-2,-12}
    ));


// Run all tests
int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}