#include <sstream>

#include <gtest/gtest.h>
#include "../sph/setup.hpp"

TEST(ConfigMapTest, ReadGoodConfig) {
    // Simple example of proper config: does it go to expected map?
    std::istringstream stream("good_param_1 100\ngood_param_2 200");
    ConfigMap expected_map = {
        {"good_param_1", "100"},
        {"good_param_2", "200"}
    };

    ConfigMap actual_map = parse_config(stream);

    EXPECT_EQ(expected_map, actual_map);
}

TEST(ConfigMapTest, ReadDupedConfig) {
    // Should not replace value
    std::istringstream stream("good_param_1 100\ngood_param_1 200");
    ConfigMap expected_map = {
        {"good_param_1", "100"}
    };

    ConfigMap actual_map = parse_config(stream);

    EXPECT_EQ(expected_map, actual_map);
}

TEST(ConfigMapTest, ReadBadConfig) {
    // Should give status 1 as bad_param has no value
    std::istringstream stream("bad_param");
    EXPECT_EXIT(
        parse_config(stream), testing::ExitedWithCode(1),
        "Parsing error on line [0-9]* of config file."
    );
}

TEST(ConfigClassTest, ReadConfigMap) {
    ConfigMap config_map = {
        {"n_part", "3"},
        {"d_unit", "4"}
    };

    Config config(config_map);

    EXPECT_EQ(config.n_part, 3);
    EXPECT_EQ(config.d_unit, 4);
}