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

void AccelerationCalculator::operator()(Particle &p_i, const ParticleVector &p_vec) {
    // I have tried to use variable names that correspond to how this equation is typeset in the
    // Bate thesis. Pr = pressure, p = particle, rho = density, W = weighting
    double Pr_i = this->pressure_isothermal(p_i);
    double Pr_rho_i = Pr_i / std::pow(p_i.density, 2);

    double acc = 0;

    for (Particle p_j : p_vec) {
        if (p_j != p_i) {
            // It probably isn't efficient to declare so many variables, but it makes the code more
            // readable, and they probably get optimized out by the compiler anyway(?)
            double r_ij = p_j.pos - p_i.pos;
            double q = r_ij / Kernel::smoothing_length(this->config);
            double grad_W = Kernel::d_kernel(q);

            double Pr_j = this->pressure_isothermal(p_j);
            double Pr_rho_j = Pr_j / std::pow(p_j.density, 2);
            
            acc += p_j.mass * (Pr_rho_i + Pr_rho_j) * grad_W;
        }
    }

    p_i.acc = acc;
}

double AccelerationCalculator::pressure_isothermal(const Particle &p) {
    double c_s = this->sound_speed();
    return std::pow(c_s, 2) * p.density;
}

double AccelerationCalculator::sound_speed() {
    // 10 m.s.^-1 in code units
    return 10 / (this->config.d_unit / this->config.t_unit);
}

#pragma endregion