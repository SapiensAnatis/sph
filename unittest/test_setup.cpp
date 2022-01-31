#include <fstream>

#include <gtest/gtest.h>
#include "../src/setup.hpp"

TEST(SetupTest, ReadGoodConfig) {
    std::ifstream config_stream;
    config_stream.open("config_good.txt");

    int cfg_read_status;
    parse_config(config_stream, cfg_read_status);

    EXPECT_EQ(cfg_read_status, 0);
}

TEST(SetupTest, ReadBadConfig) {
    std::ifstream config_stream;
    config_stream.open("config_bad.txt");

    int cfg_read_status = 4;
    parse_config(config_stream, cfg_read_status);

    EXPECT_EQ(cfg_read_status, 1);
}