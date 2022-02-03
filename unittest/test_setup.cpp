#include <sstream>

#include <gtest/gtest.h>
#include "../sph/setup.hpp"


// Helper method to build config file stringstreams using a vector. A bit more readable than
// a string with a bunch of \ns in it
std::istringstream build_config_stream(std::vector<std::string> line_vector) {
    // istringstream is read-only, so can't accumulate onto that
    std::string accum;
    
    for (std::string s : line_vector) {
        accum += s;
        accum += "\n";
    }

    return std::istringstream(accum);
}

TEST(ConfigClassTest, ReadConfigStream) {
    std::istringstream stream = build_config_stream({
        "d_unit 2",
        "n_part 4"
    });

    Config config(stream);

    EXPECT_EQ(config.d_unit, 2);
    EXPECT_EQ(config.n_part, 4);
}

TEST(ConfigClassTest, ReadBadConfig) {
    // Should give status 1 as bad_param has no value
    std::istringstream stream = build_config_stream({
        "bad_param",
        "d_unit 3000",
        "n_part 4"
    });

    EXPECT_EXIT(
        Config config(stream), testing::ExitedWithCode(1),
        "Parsing error on line 1 of config file."
    );
}

TEST(ConfigClassTest, ReadDupedConfig) {
    // Duplicated value
    std::istringstream stream = build_config_stream({
        "d_unit 2",
        "d_unit 3000",
        "n_part 4"
    });
    
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