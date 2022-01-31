#include <sstream>

#include <gtest/gtest.h>
#include "../src/setup.hpp"

TEST(SetupTest, ReadGoodConfig) {
    // Simple example of proper config: does it go to expected map?
    std::istringstream stream("good_param_1 100\ngood_param_2 200");
    ConfigMap expected_map = {
        {"good_param_1", "100"},
        {"good_param_2", "200"}
    };

    int status;
    ConfigMap actual_map = parse_config(stream, status);

    EXPECT_EQ(status, 0);
    EXPECT_EQ(expected_map, actual_map);
}

TEST(SetupTest, ReadDupedConfig) {
    // Should not replace value
    std::istringstream stream("good_param_1 100\ngood_param_1 200");
    ConfigMap expected_map = {
        {"good_param_1", "100"}
    };

    int status;
    ConfigMap actual_map = parse_config(stream, status);

    EXPECT_EQ(status, 0);
    EXPECT_EQ(expected_map, actual_map);
}


TEST(SetupTest, ReadBadConfig) {
    // Should give status 1 as bad_param has no value
    std::istringstream stream("bad_param");
    
    int status;
    parse_config(stream, status);

    EXPECT_EQ(status, 1);
}