/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.cpp implements the class methods defined in calculators.hpp.
 */

#include <cmath>
#include <exception>

#include "calculators.hpp"
#include "kernel.hpp"

#pragma region DensityCalculator

void DensityCalculator::operator()(Particle &p, const ParticleVector &p_vec) {
    double h = Kernel::smoothing_length(this->config);
    double density = 0;

    for (Particle p_i : p_vec) {
        if (p_i != p) { // i != j
            double q = std::abs(p_i.pos - p.pos) / h;
            double w = Kernel::kernel(q);
            
            density += p_i.mass * (w / h);
        }
    }

    p.density = density;
}

#pragma endregion
#pragma region MomentumCalculator

void AccelerationCalculator::operator()(Particle &p, const ParticleVector &p_vec) {
    throw std::logic_error("Not implemented yet!");
}

#pragma endregion