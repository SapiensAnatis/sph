/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.cpp implements the class methods defined in calculators.hpp, and also has a few
 * free-floating non-exposed functions that are used by the GSL setup.
 */


#include <cmath>
#include <exception>
#include <iostream>

#include "define.hpp"
#include "calculators.hpp"
#include "kernel.hpp"
#include "smoothing_length.hpp"

double Calculator::grad_W(const Particle &p_i, const Particle &p_j, double h) {
    double r_ij = p_i.pos - p_j.pos;
    // The unit vector is +-1, depending on the sign of the vector, because we are in 1D
    double r_ij_unit = (r_ij > 0) ? 1 : -1;

    double q = std::abs(r_ij) / h;
    double grad_W = dkernel_dq(q) * r_ij_unit; // Rosswog 2009 eq. 25

    return grad_W;
}

#pragma region DensityCalculator

#ifdef USE_VARIABLE_H
// Variable smoothing length implementation
void DensityCalculator::operator()(Particle &p) {
    std::pair<double, double> root_result = rootfind_h(p, p_arr, config);

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

        double new_density = calc_density(p, p.h, p_arr.get(), config.n_part);

        std::cout << "[WARN] Recalculated above negative density as " << new_density << std::endl;
        p.density = new_density;
    }
}
#endif

#ifndef USE_VARIABLE_H
void DensityCalculator::operator()(Particle &p) {
    double d_sum = 0;

    for (int i = 0; i < config.n_part; i++) {
        Particle p_j = p_arr[i];
        double q = std::abs(p.pos - p_j.pos) / CONSTANT_H;
        double w = kernel(q);

        d_sum += p.mass * (w / CONSTANT_H);
    }

    p.density = d_sum;
}

#endif

#pragma endregion
#pragma region AccelerationCalculator

// Version of acceleration calculation that accounts for variable smoothing length, by adding in
// 'omega terms' (Rosswog eqns. 118-121). If USE_VARIABLE_H isn't defined then calc_omega() just
// returns 1, simplifying it to the standard SPH expression.
void AccelerationCalculator::operator()(Particle &p_i) {
    if (p_i.type == Ghost)
        return;

    // I have tried to use variable names that correspond to how this equation is typeset in the
    // Bate thesis. Pr = pressure, p = particle, rho = density, W = weight function
    double c_s;
    double Pr_i;

    // Annoyingly, in the isothermal case pressure is dependent on sound speed, but in the adiabatic
    // case, sound speed is dependent on pressure. So the order switches based on which one is used.
    if (config.pressure_calc == Isothermal) {
        c_s = sound_speed(p_i);
        Pr_i = pressure_isothermal(p_i, c_s);
    } else if (config.pressure_calc == Adiabatic) {
        Pr_i = pressure_adiabatic(p_i);
        c_s = sound_speed(p_i);
    } else {
        throw std::logic_error("Unknown pressure calculation mode!");
    }

    // Keep track of pressures as they can be used to verify the analytical solution
    p_i.pressure = Pr_i;

    double Pr_rho_i = Pr_i / std::pow(p_i.density, 2) / calc_omega(p_i, p_arr, config);

    double acc = 0;

    // Density = 0 will cause div by zero and screw everything up. Should never really happen
    ensure_nonzero_density(p_i);

    for (int i = 0; i < config.n_part; i++) {
        Particle &p_j = p_arr[i];
        ensure_nonzero_density(p_j);

        if (p_j != p_i) {
            double r_ij = p_i.pos - p_j.pos;

            double grad_W_i = grad_W(p_i, p_j, p_i.h);
            // Different smoothing length of particle j. Gradient still w.r.t. i
            double grad_W_j = grad_W(p_i, p_j, p_j.h);

            double Pr_j;
            if (config.pressure_calc == Isothermal)
                Pr_j = pressure_isothermal(p_j, c_s);
            else if (config.pressure_calc == Adiabatic)
                Pr_j = pressure_adiabatic(p_j);
            else
                throw std::logic_error("Unknown pressure calculation mode!");

            double Pr_rho_j = Pr_j / std::pow(p_j.density, 2) / calc_omega(p_j, p_arr, config);

            // TODO: Figure out how artificial viscosity fits into this equation!
            double visc_ij = artificial_viscosity(p_i, p_j, r_ij, p_i.h, c_s);

            // Rosswog 2009 eqn 120 (plus viscosity?) I know it's horrible, I'm sorry
            double to_add = -p_j.mass * ((grad_W_i * Pr_rho_i) + (grad_W_i * visc_ij) + (grad_W_j * Pr_rho_j));
            acc += to_add;
        }
    }

    p_i.acc = acc;
}

double AccelerationCalculator::pressure_isothermal(const Particle &p, double c_s) {
    return std::pow(c_s, 2) * p.density;
}

double AccelerationCalculator::pressure_adiabatic(const Particle &p) {
    return (GAMMA - 1) * p.u * p.density;
}

double AccelerationCalculator::sound_speed(Particle p) {
    if (config.pressure_calc == Isothermal)
        return 1;
    else if (config.pressure_calc == Adiabatic)
        return std::sqrt(GAMMA * p.pressure / p.density);
    else
        throw std::logic_error("Unknown pressure calculation mode!");
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

#pragma region EnergyCalculator

void EnergyCalculator::operator()(Particle &p) {
    // Price 2012 eqn. 35
    double omega = calc_omega(p, p_arr, config);
    double Pr_rho = p.pressure / (omega * std::pow(p.density, 2));

    double sum = 0;
    for (int i = 0; i < config.n_part; i++) {
        Particle p_j = p_arr[i];

        double r_ij = p.pos - p_j.pos;
        double v_ij = p.vel - p_j.vel;
        double c_s = sound_speed(p);
        
        double visc = Pr_rho + 0.5 * artificial_viscosity(p, p_j, r_ij, p.h, c_s);

        // Citation: Private correspondence between Loren-Aguilar and Seal (??)
        sum += p_j.mass * visc * v_ij * grad_W(p, p_j, p.h);
    }

    p.du_dt = sum;

    if (p.id == 42) {
        std::cout << "du_dt 42 " << p.du_dt << std::endl;
    }
}

#pragma endregion