/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.cpp implements the class methods defined in calculators.hpp.
 */

#include <cmath>
#include <exception>
#include <iostream>

#include "calculators.hpp"

#pragma region Kernel

double Kernel::kernel(double q) {
    double w;
    double sigma = 2.0/3.0;

    // Third-order Schoenberg B-spline (M_4), as detailed in Price (2010) pg. 3.
    if (q < 1) {
        w = 0.25 * std::pow(2 - q, 3) - std::pow(1 - q, 3);
        return sigma * w;
    } else if (q >= 1 && q < 2) {
        w = 0.25 * std::pow(2 - q, 3);
        return sigma * w;
    } else {
        // > 2
        return 0;
    }
}

double Kernel::d_kernel(double q) {
    double w;
    double sigma = 2.0/3.0;

    // Derivative of kernel(), analytically
    if (q < 1) {
        w = 0.25 * -3 * std::pow(2 - q, 2) + 3 * std::pow(1 - q, 2);
        return sigma * w;
    } else if (q >= 1 && q < 2) {
        w = 0.25 * -3 * std::pow(2 - q, 2);
        return sigma * w;
    } else {
        // > 2
        return 0;
    }
}

// This returns a constant for now, but may later use variable smoothing lengths
double Kernel::smoothing_length(const Config &c) {
    return c.smoothing_length * c.d_unit;
}

#pragma endregion
#pragma region DensityCalculator

void DensityCalculator::operator()(Particle &p_i) {
    double h = Kernel::smoothing_length(config);
    double density = 0;

    for (int i = 0; i < config.n_part; i++) {
        Particle &p_j = p_all[i];
        if (p_i != p_j) { // i != j
            double q = std::abs(p_i.pos - p_j.pos) / h;
            double w = Kernel::kernel(q);
            
            density += p_i.mass * (w / h);
        }
    }

    p_i.density = density;
}

#pragma endregion
#pragma region AccelerationCalculator

void AccelerationCalculator::operator()(Particle &p_i) {
    // I have tried to use variable names that correspond to how this equation is typeset in the
    // Bate thesis. Pr = pressure, p = particle, rho = density, W = weight function
    double c_s = sound_speed();
    double Pr_i = pressure_isothermal(p_i, c_s);
    double Pr_rho_i = Pr_i / std::pow(p_i.density, 2);
    double h = Kernel::smoothing_length(config);

    double acc = 0;

    // Density = 0 will cause div by zero and screw everything up. Should never really happen
    ensure_nonzero_density(p_i);

    for (int i = 0; i < config.n_part; i++) {
        Particle &p_j = p_all[i];
        ensure_nonzero_density(p_j);

        if (p_j != p_i) {
            // It probably isn't efficient to declare so many variables, but it makes the code more
            // readable, and they probably get optimized out by the compiler anyway(?)
            double r_ij = std::abs(p_i.pos - p_j.pos);
            double q = r_ij / h;
            double grad_W = Kernel::d_kernel(q);

            double Pr_j;
            if (config.pressure_calc == Isothermal)
                Pr_j = pressure_isothermal(p_j, c_s);
            else
                throw std::invalid_argument("Adiabatic EoS not yet implemented!");

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
    // 10 m.s^-1 in code units. Normally we multiply by the unit to use code units,
    // but this was originally given in m.s^-1 so the process is reversed
    return 10 / (config.d_unit / config.t_unit);
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
    
    double eta_sq = eta_coeff * std::pow(h, 2);
    double mu_ij = (h * dot)/(std::pow(r_ij, 2) + eta_sq);
    double rho_ij = (p_i.density + p_j.density) / 2;

    double result = 0;
    result += -alpha * c_s * mu_ij;
    result += beta * std::pow(mu_ij, 2);
    result /= rho_ij;

    return result;
}

void AccelerationCalculator::ensure_nonzero_density(const Particle &p) {
    if (p.density < CALC_EPSILON) {
        std::cerr << "[ERROR] Particle had density less than epsilon " << CALC_EPSILON << std::endl;
        std::cerr << "[ERROR] Particle id: " << p.id << " has density: " << p.density << std::endl;
        // std::cerr << "[ERROR] Current time " << No way of accessing that, oops!

        throw new std::logic_error("Particle had density less than epsilon!");
    }
}

#pragma endregion