/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * test_calculator.cpp defines unit tests for the calculation methods, such as those for density
 * and acceleration.
 */

#include <gtest/gtest.h>

#include "../sph/particle.hpp"
#include "../sph/calculators.hpp"
#include "../sph/setup.hpp"

// Shared objects between test suites
class CalcTestFixture : public ::testing::Test {
    protected:
        ParticleVector p_vec;
        Config config;

        CalcTestFixture() {
            // Particle vector with 3 static particles: one at -0.5, one at 0, and one at 0.5
            this->p_vec = ParticleVector {
                Particle(-0.5, 0, 1),
                Particle(0, 0, 1),
                Particle(0.5, 0, 1)
            };

            // Create config
            this->config = Config(); 

            this->config.n_part = 3;
            this->config.d_unit = 1;
            this->config.t_unit = 1;
            this->config.mass = 1;
            this->config.limit = 1;
            this->config.v_0 = 10;
            this->config.smoothing_length = 1;
        }
};

TEST_F(CalcTestFixture, DensityCalc) {
    // Compare against hand-calculated values
    auto dc = DensityCalculator(config);

    for (Particle &p : p_vec) {
        dc(p, p_vec);
    }

    EXPECT_FLOAT_EQ(p_vec[0].density, 31.0/48.0);
    EXPECT_FLOAT_EQ(p_vec[1].density, 23.0/24.0);
    EXPECT_FLOAT_EQ(p_vec[2].density, 31.0/48.0);
}

// TODO: Acceleration calculation test