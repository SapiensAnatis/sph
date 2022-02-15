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
    double c_s = this->sound_speed();
    double Pr_i = this->pressure_isothermal(p_i, c_s);
    double Pr_rho_i = Pr_i / std::pow(p_i.density, 2);
    double h = Kernel::smoothing_length(this->config);

    double acc = 0;

    for (Particle p_j : p_vec) {
        if (p_j != p_i) {
            // It probably isn't efficient to declare so many variables, but it makes the code more
            // readable, and they probably get optimized out by the compiler anyway(?)
            double r_ij = std::abs(p_i.pos - p_j.pos);
            double q = r_ij / h;
            double grad_W = Kernel::d_kernel(q);

            double Pr_j = this->pressure_isothermal(p_j, c_s);
            double Pr_rho_j = Pr_j / std::pow(p_j.density, 2);

            double visc_ij = artificial_viscosity(p_i, p_j, r_ij, h, c_s);
            
            acc += p_j.mass * (Pr_rho_i + Pr_rho_j + visc_ij) * grad_W;
        }
    }

    p_i.acc = acc;
}

double AccelerationCalculator::pressure_isothermal(const Particle &p, double c_s) {
    return std::pow(c_s, 2) * p.density;
}

double AccelerationCalculator::sound_speed() {
    // 10 m.s.^-1 in code units
    return 10 / (this->config.d_unit / this->config.t_unit);
}

double AccelerationCalculator::artificial_viscosity(
    const Particle &p_i, 
    const Particle &p_j, 
    double r_ij,
    double h,
    double c_s
) {
    double v_ij = p_i.vel - p_j.vel;
    double dot = v_ij * r_ij;
    
    if (dot > 0) 
        return 0;
    
    double eta_sq = this->eta_coeff * std::pow(h, 2);
    double mu_ij = (h * dot)/(std::pow(r_ij, 2) + eta_sq);
    double rho_ij = (p_i.density + p_j.density) / 2;

    double result = 0;
    result += -this->alpha * c_s * mu_ij;
    result += this->beta * std::pow(mu_ij, 2);
    result /= rho_ij;

    return result;
}

#pragma endregion