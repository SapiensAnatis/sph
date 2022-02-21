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
            this->config.pressure_calc = Isothermal;
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

/*
TEST_F(CalcTestFixture, AccelerationCalc) {
    auto ac = AccelerationCalculator(config);
    // Note that the test fixture is re-constructed for each test; there is no persistency between
    // test suites. So we need to set the densities to start with. We use fixed values rather than
    // DensityCalculator, so this one won't fail if the density calculation starts to fail.
    p_vec[0].density = 31.0/48.0;
    p_vec[1].density = 23.0/24.0;
    p_vec[2].density = 31.0/48.0;

    // Hand-calculated isothermal pressures for c_s = 10
    double P_0 = 100 * p_vec[0].density;
    double P_1 = 100 * p_vec[1].density;
    double P_2 = 100 * p_vec[2].density;

    // grad_W(0.5)
    double g_W_0p5 = -0.625;
    // grad_W(1)
    double g_W_1 = -0.5;

    // Artificial viscosity
    // _1: Between central particle and either edge particle
    // _2: Between both edge particles

    
    for (Particle &p : p_vec) {
        ac(p, p_vec);
    }

    
    EXPECT_FLOAT_EQ(p_vec[0].vel, );
    EXPECT_FLOAT_EQ(p_vec[1].vel, );
    EXPECT_FLOAT_EQ(p_vec[2].vel, );
    
}

This was starting to be a gigantic pain so I decided to just do the density one and leave the rest
of the code to be tested by reproducing the analytical solution

*/
