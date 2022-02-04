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
        "n_part 4",
        "d_unit 3",
        "t_unit 2",
        "limit 1"
    });

    Config config(stream);

    EXPECT_EQ(config.n_part, 4);
    EXPECT_EQ(config.d_unit, 3);
    EXPECT_EQ(config.t_unit, 2);
    EXPECT_EQ(config.limit, 1);
}

TEST(ConfigClassTest, ReadBadConfig) {
    // Should give status 1 as bad_param has no value
    std::istringstream stream = build_config_stream({
        "bad_param",
        "d_unit 3000",
        "n_part 4",
        "t_unit 2",
        "limit 1"
    });

    EXPECT_EXIT(
        Config config(stream), testing::ExitedWithCode(1),
        "Parsing error on line 1 of config file."
    );
}

TEST(ConfigClassTest, ReadDupedConfig) {
    // Duplicated value
   std::istringstream stream = build_config_stream({
        "n_part 4",
        "d_unit 3",
        "t_unit 2",
        "limit 1",
        "d_unit 300000"
    });
    
    Config config(stream);

    EXPECT_EQ(config.n_part, 4);
    EXPECT_EQ(config.d_unit, 3);
}

TEST(ConfigClassTest, ReadMissingConfig) {
    // Missing value(s)
    std::istringstream stream("d_unit 2");
    
    EXPECT_EXIT(
        Config config(stream), testing::ExitedWithCode(1),
        "Failed to find a definition for property 'n_part'"
    );
}