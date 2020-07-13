#include "gtest/gtest.h"
#include "DataUnit.h"


// Parameterised test to check valid header parsing
// for operator >> (std::istream& stream, DataUnit &d)
struct header_contents{
    unsigned int type;
    unsigned int next_parse_offset;
    unsigned int prev_parse_offset;
    unsigned int type_result;
};

struct HeaderParseTest : testing::Test, testing::WithParamInterface<header_contents>{
    DataUnit du;
};

TEST_P(HeaderParseTest, CheckValidHeaders){
    
    auto as = GetParam();

    Bytes next(4,as.next_parse_offset); 
    Bytes prev(4,as.prev_parse_offset);

    std::stringstream ss;
    ss  << 'B' //0x41
        << 'B' //0x41
        << 'C' //0x42
        << 'D' //0x43
        << (unsigned char)as.type 
        << next
        << prev;

    ss >> du;

    EXPECT_EQ(du.type, as.type_result);
    EXPECT_EQ((unsigned int)du.next_parse_offset, as.next_parse_offset);
    EXPECT_EQ((unsigned int)du.prev_parse_offset, as.prev_parse_offset);
    EXPECT_EQ(du.length(), (int)as.next_parse_offset - 13);

};

INSTANTIATE_TEST_SUITE_P(Default, HeaderParseTest,
    testing::Values(
        header_contents{0x00,0,0,SEQUENCE_HEADER},
        header_contents{0x10,20,20,END_OF_SEQUENCE},
        header_contents{0x20,128,0,AUXILIARY_DATA},
        header_contents{0x30,0,3245,PADDING_DATA},
        header_contents{0xC8,23,0,LD_PICTURE},
        header_contents{0xE8,13,13,HQ_PICTURE},
        header_contents{0xCC,13,13,LD_FRAGMENT},
        header_contents{0xEC,13,0,HQ_FRAGMENT}  
    ));

// Check correct exception is thrown when the prefix is incorrect
TEST(HeaderParseTest, ThrowErrorIfPrefixIncorrect){
    DataUnit du;
    std::stringstream ss;
    ss  << 'A' 
        << 'B' 
        << 'C' 
        << 'D';

    try{
        ss >> du;
    } catch (std::logic_error& exception){
        EXPECT_EQ(std::string(exception.what()), 
        "Read bytes do not match expected parse_info_header.");
    }
}

// Check correct exception is thrown when the type is incorrect
TEST(HeaderParseTest, ThrowErrorIfTypeIncorrect){
    DataUnit du;
    std::stringstream ss;
    ss  << 'B' 
        << 'B' 
        << 'C' 
        << 'D'
        << (unsigned char) 0xFF;
        
    try{
        

        ss >> du;

    } catch (std::logic_error& exception){
        EXPECT_EQ(std::string(exception.what()), 
        "Stream Error: Unknown data unit type.");
    }
}



// Run all tests
int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}