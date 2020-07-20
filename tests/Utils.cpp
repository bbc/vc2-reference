#include "gtest/gtest.h"
#include "Utils.h"

// Parameterised test to check getting the picture number
struct picture_values{
    int fieldNumber;
    unsigned long long frameNumber;
    const int fieldsPerFrame;
    unsigned long result;
};

struct GetPictureNumberTest : testing::Test, testing::WithParamInterface<picture_values>{};

TEST_P(GetPictureNumberTest, CheckNormalOperation){
    auto as = GetParam();
    EXPECT_EQ(utils::getPictureNumber(as.fieldNumber, as.frameNumber, as.fieldsPerFrame), as.result);
};

INSTANTIATE_TEST_SUITE_P(Default, GetPictureNumberTest,
    testing::Values(
        picture_values{0,0,1,0},
        picture_values{1,0,1,1},
        picture_values{2,0,2,2},
        picture_values{1,1,1,2},
        picture_values{2,1,2,4},
        picture_values{1,2,2,5},
        picture_values{0,((1ULL)<<32)-1 ,1,((1ULL)<<32)-1},  // 2^32-1
        picture_values{0,(1ULL)<<32,1,0}         
    ));
    
// Parameterised test to check getting the picture number exceptions
struct picture_exceptions{
    int fieldNumber;
    unsigned long long frameNumber;
    const int fieldsPerFrame;
    std::string exceptionString;
};

struct GetPictureNumberExceptionTest : testing::Test, testing::WithParamInterface<picture_exceptions>{};

TEST_P(GetPictureNumberExceptionTest, CheckExeptionsCorrectlyThrown){
    auto as = GetParam();
        try{
        utils::getPictureNumber(as.fieldNumber, as.frameNumber, as.fieldsPerFrame);
        FAIL(); // Fail if an exception is not thrown
    } catch (std::logic_error& exception){
        EXPECT_EQ(std::string(exception.what()), as.exceptionString);
    }
};

INSTANTIATE_TEST_SUITE_P(Default, GetPictureNumberExceptionTest,
    testing::Values(
        picture_exceptions{-5,0,1,"field number should be positive"},
        picture_exceptions{2,0,1,"field number exceeds number of fields per frame"},
        picture_exceptions{0,0,3,"number of fields per frame should be 1 (progressive) or 2 (interlaced)"}
    ));


// Run all tests
int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
