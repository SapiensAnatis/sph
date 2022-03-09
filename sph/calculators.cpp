/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.cpp implements the class methods defined in calculators.hpp, and also has a few
 * free-floating non-exposed functions that are used by the GSL setup.
 */


#include <cmath>
#include <exception>
#include <iostream>

#include "calculators.hpp"
#include "kernel.hpp"
#include "smoothing_length.hpp"

#pragma region DensityCalculator

void DensityCalculator::operator()(Particle &p) {
    std::pair<double, double> root_result = rootfind_h(p, p_all, config);

    // Retrieve smoothing length (with sanity check)
    double h = root_result.first;
    if (h < 0) {
        std::cerr << "[ERROR] Smoothing length root-finding for particle id: " << p.id << 
                  " returned negative smoothing length: " << h << std::endl;
        throw new std::logic_error("Root-finding returned negative smoothing length");
    }

    p.h = h;

    double density = root_result.second;
    // Sometimes it finds a solution where density < 0, which is quite bizarre.
    if (density > 0) {
        p.density = density;
    } else {
        // We can probably resolve this by recalculating
        std::cout << "[WARN] Smoothing length root-finding for particle id: " << p.id << 
                  " returned negative density: " << density << std::endl;

        double new_density = calc_density(p, p.h, p_all.get(), config.n_part);

        std::cout << "[WARN] Recalculated above negative density as " << new_density << std::endl;
        p.density = new_density;
    }
}

#pragma endregion
#pragma region AccelerationCalculator

void AccelerationCalculator::operator()(Particle &p_i) {
    if (p_i.type == Ghost)
        return;

    // I have tried to use variable names that correspond to how this equation is typeset in the
    // Bate thesis. Pr = pressure, p = particle, rho = density, W = weight function
    double c_s = sound_speed();
    double Pr_i = pressure_isothermal(p_i, c_s);
    p_i.pressure = Pr_i;
    double Pr_rho_i = Pr_i / std::pow(p_i.density, 2);

    double acc = 0;

    // Density = 0 will cause div by zero and screw everything up. Should never really happen
    ensure_nonzero_density(p_i);

    for (int i = 0; i < config.n_part; i++) {
        Particle &p_j = p_all[i];
        ensure_nonzero_density(p_j);

        if (p_j != p_i) {
            double r_ij = p_i.pos - p_j.pos;

            // The unit vector is +-1, depending on the sign of the vector, because we are in 1D
            double r_ij_unit = (r_ij > 0) ? 1 : -1;

            double q = std::abs(r_ij) / p_i.h;
            double grad_W = d_kernel(q) * r_ij_unit; // Rosswog 2009 eq. 25

            double Pr_j;
            if (config.pressure_calc == Isothermal)
                Pr_j = pressure_isothermal(p_j, c_s);
            else
                throw std::invalid_argument("Adiabatic EoS not yet implemented!");

            double Pr_rho_j = Pr_j / std::pow(p_j.density, 2);

            double visc_ij = artificial_viscosity(p_i, p_j, r_ij, p_i.h, c_s);
            double to_add = -p_j.mass * (Pr_rho_i + Pr_rho_j + visc_ij) * grad_W;
            acc += to_add;
        }
    }

    p_i.acc = acc;
}

double AccelerationCalculator::pressure_isothermal(const Particle &p, double c_s) {
    return std::pow(c_s, 2) * p.density;
}

double AccelerationCalculator::sound_speed() {
    // 1 m.s^-1
    return 1;
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