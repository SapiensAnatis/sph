/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * test_setup.cpp contains unit tests relating to the initial setup of the program. This includes
 * checking that configuration files are read correctly (and that the expected behaviour occurs
 * when configuration files are malformed), and that the Particle vector is setup correctly.
 */

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
        "mass 10",
        "limit 1",
        "v_0 8",
        "smoothing_length 1",
        "pressure_calc 0"
    });

    auto config_reader = ConfigReader(stream);
    Config config = config_reader.config;

    EXPECT_EQ(config.n_part, 4);
    EXPECT_EQ(config.d_unit, 3);
    EXPECT_EQ(config.t_unit, 2);
    EXPECT_EQ(config.mass, 10);
    EXPECT_EQ(config.limit, 1);
    EXPECT_EQ(config.v_0, 8);
    EXPECT_EQ(config.smoothing_length, 1);
    EXPECT_EQ(config.pressure_calc, Isothermal);
}

TEST(ConfigClassTest, ReadBadConfig) {
    // Should give status 1 as bad_param has no value
    std::istringstream stream = build_config_stream({
        "bad_param",
        "d_unit 3000",
        "n_part 4",
        "mass 1"
        "t_unit 2",
        "limit 1",
        "v_0 10",
        "smoothing_length 1",
        "pressure_calc 0"
    });

    EXPECT_EXIT(
        ConfigReader config_reader(stream), testing::ExitedWithCode(1),
        "Parsing error on line 1 of config file."
    );
}

TEST(ConfigClassTest, ReadDupedConfig) {
    // Duplicated value
    std::istringstream stream = build_config_stream({
        "n_part 4",
        "d_unit 3",
        "t_unit 2",
        "mass 10",
        "limit 1",
        "d_unit 300000",
        "v_0 12",
        "v_0 19",
        "smoothing_length 1",
        "pressure_calc 0"
    });
    
    auto config_reader = ConfigReader(stream);
    Config config = config_reader.config;

    EXPECT_EQ(config.n_part, 4);
    EXPECT_EQ(config.d_unit, 3);
    EXPECT_EQ(config.v_0, 12);
}

TEST(ConfigClassTest, ReadMissingConfig) {
    // Missing value(s)
    std::istringstream stream("d_unit 2");
    
    EXPECT_EXIT(
        ConfigReader config_reader(stream), testing::ExitedWithCode(1),
        "Failed to find a definition for property 'n_part'"
    );
}

/*
Since switching to array pointers, sizeof(p_arr) is always 8 -- I think it decays to a regular
pointer and renders size info inaccessible.

TEST(ParticleSetup, CorrectNParticles) {
    // Check the right number of particles are made
    std::istringstream stream = build_config_stream({
        "n_part 20",
        "d_unit 1",
        "t_unit 1",
        "mass 10",
        "limit 1",
        "v_0 10",
        "smoothing_length 1",
        "pressure_calc 0"
    });

    auto config = Config(stream);
    ParticleArrayPtr p_arr(new Particle[config.n_part]);

    EXPECT_EQ(sizeof(p_arr), 20 * sizeof(Particle));
} 
*/

TEST(ParticleSetup, CorrectMass) {
    std::istringstream stream = build_config_stream({
        "n_part 25",
        "d_unit 1",
        "t_unit 1",
        "mass 15",
        "limit 1",
        "v_0 12",
        "smoothing_length 1",
        "pressure_calc 0"
    });

    auto config_reader = ConfigReader(stream);
    Config config = config_reader.config;

    ParticleArrayPtr p_arr(new Particle[config.n_part]);
    init_particles(config, p_arr);

    for (int i = 0; i < config.n_part; i++) {
        EXPECT_EQ(p_arr[i].mass, 15);
    }
}


TEST(ParticleSetup, CorrectVZero) {
    // Check that all particles are given either plus or minus v_0. Velocity should be positive if
    // position is negative, but velocity should be negative if position is positive.
    std::istringstream stream = build_config_stream({
        "n_part 20",
        "d_unit 1",
        "t_unit 1",
        "mass 10",
        "limit 1",
        "v_0 12",
        "smoothing_length 1",
        "pressure_calc 0"
    });

    auto config_reader = ConfigReader(stream);
    Config config = config_reader.config;
    
    ParticleArrayPtr p_arr(new Particle[config.n_part]);
    init_particles(config, p_arr);

    for (int i = 0; i < config.n_part; i++) {
        // For the edge case p.pos == 0 (exceedingly unlikely because RNG), see line 122-ish of
        // setup.cpp init_particles(); v will be negative
        if (p_arr[i].pos < 0) {
            EXPECT_EQ(p_arr[i].vel, 12);
        } else {
            EXPECT_EQ(p_arr[i].vel, -12);
        }
    }
}

/*
// The random nature of the particle distribution means that if this test *does* fail, it will
probably be sporadically and very rarely. So it's pretty badly designed.
// Removed for this reason
TEST(ParticleSetup, WithinBounds) {
    // Check that all particles are within the specified boundary
    std::istringstream stream = build_config_stream({
        "n_part 20",
        "d_unit 2",
        "t_unit 1",
        "limit 2"
    });

    auto config = Config(stream);
    ParticleVector pv = init_particles(config);

    double max_expected = config.limit * config.d_unit;

    for (Particle p : pv) {
        EXPECT_LE(std::abs(p.pos), max_expected);
    }
}
*/