/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * calculators.cpp implements the class methods defined in calculators.hpp, and also has a few
 * free-floating non-exposed functions that are used by the GSL setup.
 */

// Maximum number of iterations for root-finding of smoothing length
#define H_MAX_ITER 1000
// Epsilon for root-finding of smoothing length -- used to decide when to declare success
#define H_EPSILON 1e-3
// Show root-finding debug messages
#define H_DEBUG

#include <cmath>
#include <exception>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

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
    return c.smoothing_length;
}

#pragma endregion
#pragma region DensityCalculator
// The code for root-finding the variable smoothing length is largely based on the Rosenbrock
// example in the GSL documentation:
// https://www.gnu.org/software/gsl/doc/html/multiroots.html#examples

// Params for root-finding method
struct params
{
    Particle* p; // Particle in question
    Particle* p_arr; // Pointer to array of particles
    int n_part; // Length of above array
    double eta; // Smoothing length parameter; see Price 2010 eq. 10
};

// This method and smoothing_f are not in the headers, since they aren't public interfaces. They are
// solely used by DensityCalculator::operator().
double calc_density(const Particle &p, const double h, const Particle* p_arr, const int n_part) {
    double d_sum = 0;
    for (int i = 0; i < n_part; i++) {
        Particle p_j = *(p_arr + i);
        double q = std::abs(p.pos - p_j.pos) / h;
        double w = Kernel::kernel(q);

        d_sum += p.mass * (w / h);
    }

    return d_sum;
}


// Method defining the system of equations. Can't be a member of DensityCalculator because it
// changes the type of the function making it unacceptable to the GSL setup
int smoothing_f(const gsl_vector* x, void* params, gsl_vector* f) {
    // Get parameters
    Particle p = *((struct params*)params)->p;
    Particle* p_arr = ((struct params*)params)->p_arr;
    int n_part = ((struct params*)params)->n_part;
    double eta = ((struct params*)params)->eta;

    // Smoothing length
    const double x0 = gsl_vector_get(x, 0);

    // Density
    const double x1 = gsl_vector_get(x, 1);

    // Calculate density via sum over other particles
    double density = 0;
    for (int i = 0; i < n_part; i++) {
        Particle p_j = *(p_arr + i);
        double q = std::abs(p.pos - p_j.pos) / x0;
        double w = Kernel::kernel(q);

        density += p.mass * (w / x0);
    }

    // Smoothing length equation: h - eta(m/rho) = 0
    const double f0 = x0 - eta*(p.mass / x1);
    // Density equation: rho - sum = 0
    const double f1 = x1 - density;

    gsl_vector_set(f, 0, f0);
    gsl_vector_set(f, 1, f1);

    return GSL_SUCCESS;
}

void print_state (size_t iter, gsl_multiroot_fsolver* s)
{
  printf("iter = %3lu x = % .3f % .3f f(x) = % .3e % .3e\n",
         iter,
         gsl_vector_get (s->x, 0),
         gsl_vector_get (s->x, 1),
         gsl_vector_get (s->f, 0),
         gsl_vector_get (s->f, 1));
}

// operator() sets up a GSL solver to employ a root-finding method and vary the smoothing length.
// I'm using a root-finding algorithm that doesn't require derivatives because taking the derivative
// of the density sounds like a nightmare
void DensityCalculator::operator()(Particle &p_i) {
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    struct params param = {
        &p_i,
        p_all.get(),
        config.n_part,
        config.smoothing_length
    };

    gsl_multiroot_function f = {&smoothing_f, 2, &param};

    // Initial guesses for {smoothing length, density}
    // Mean particle spacing * neighbour parameter
    double h_guess = ((config.limit * 2) / config.n_part) * config.smoothing_length; 
    // Number density * mass
    double rho_guess = (config.n_part / (config.limit * 2)) * config.mass;

    double x_init[2] = {h_guess, rho_guess};
    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, 2);
    gsl_multiroot_fsolver_set(s, &f, x);

    // print_state(iter, s);

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        if (status)
            break;

        // print_state(iter, s);
        
        status = gsl_multiroot_test_residual(s->f, H_EPSILON);
    } while (status == GSL_CONTINUE && iter < H_MAX_ITER);

    #ifdef H_DEBUG
    if (status != GSL_SUCCESS) {
        std::cout << "[WARN] Smoothing length root-finding failed for particle id " << p_i.id << 
                    " with status '" << gsl_strerror(status) << "'" << std::endl;
    }
    #endif
    
    
    // Retrieve smoothing length (with sanity check)
    double h = gsl_vector_get(s->x, 0);
    if (h < 0) {
        std::cerr << "[ERROR] Smoothing length root-finding for particle id: " << p_i.id << 
                  " returned negative smoothing length: " << h << std::endl;
        throw new std::logic_error("Root-finding returned negative smoothing length");
    }
    p_i.h = gsl_vector_get(s->x, 0);

    double density = gsl_vector_get(s->x, 1);

    // Sometimes it finds a solution where density < 0, which is quite bizarre.
    if (density > 0) {
        p_i.density = density;
    } else {
        // We can probably resolve this by recalculating
        std::cout << "[WARN] Smoothing length root-finding for particle id: " << p_i.id << 
                  " returned negative density: " << density << std::endl;

        double new_density = calc_density(p_i, p_i.h, p_all.get(), config.n_part);

        std::cout << "[WARN] Recalculated above negative density as " << new_density << std::endl;
        p_i.density = new_density;
    }

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
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
    double h = Kernel::smoothing_length(config);

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

            double q = std::abs(r_ij) / h;
            double grad_W = Kernel::d_kernel(q) * r_ij_unit; // Rosswog 2009 eq. 25

            double Pr_j;
            if (config.pressure_calc == Isothermal)
                Pr_j = pressure_isothermal(p_j, c_s);
            else
                throw std::invalid_argument("Adiabatic EoS not yet implemented!");

            double Pr_rho_j = Pr_j / std::pow(p_j.density, 2);

            double visc_ij = artificial_viscosity(p_i, p_j, r_ij, h, c_s);
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