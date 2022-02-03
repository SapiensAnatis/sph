#include <sstream>

#include <gtest/gtest.h>
#include "../sph/setup.hpp"

TEST(ConfigClassTest, ReadConfigStream) {
    std::istringstream stream("d_unit 2\nn_part 4");

    Config config(stream);

    EXPECT_EQ(config.d_unit, 2);
    EXPECT_EQ(config.n_part, 4);
}