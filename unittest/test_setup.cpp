#include <sstream>

#include <gtest/gtest.h>
#include "../sph/setup.hpp"

TEST(ConfigClassTest, ReadConfigStream) {
    std::istringstream stream("d_unit 2\nn_part 4");

    Config config(stream);

    EXPECT_EQ(config.d_unit, 2);
    EXPECT_EQ(config.n_part, 4);
}

TEST(ConfigClassTest, ReadBadConfig) {
    // Should give status 1 as bad_param has no value
    std::istringstream stream("bad_param");
    EXPECT_EXIT(
        Config config(stream), testing::ExitedWithCode(1),
        "Parsing error on line 1 of config file."
    );
}

TEST(ConfigClassTest, ReadDupedConfig) {
    // Duplicated value
    std::istringstream stream("d_unit 2\nd_unit 3000\nn_part 4");
    
    Config config(stream);

    EXPECT_EQ(config.d_unit, 2);
    EXPECT_EQ(config.n_part, 4);
}

TEST(ConfigClassTest, ReadMissingConfig) {
    // Missing value
    std::istringstream stream("d_unit 2");
    
    EXPECT_EXIT(
        Config config(stream), testing::ExitedWithCode(1),
        "Failed to find a definition for property 'n_part'"
    );
}