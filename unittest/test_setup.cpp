#include <sstream>

#include <gtest/gtest.h>
#include "../src/setup.hpp"

TEST(SetupTest, ReadGoodConfig) {
    std::istringstream stream("good_param_1 100\ngood_param_2 200");

    int status;
    parse_config(stream, status);

    EXPECT_EQ(status, 0);
}

TEST(SetupTest, ReadBadConfig) {
    std::istringstream stream("bad_param");
    
    int status;
    parse_config(stream, status);

    EXPECT_EQ(status, 1);
}


// It's not opening the files for some reason