#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <iostream>

#include "smoothing_length.hpp"
#include "define.hpp"
#include "kernel.hpp"

// Params for root-finding method
// In hindsight, I should've used a ParticleArrayPtr in this params struct, but I suppose I had an
// unconscious bias against using 'fancy' C++ stuff as the GSL documentation and examples that I
// based this code off of is designed for C and is quite spartan.
// I don't want to change it as this stage as things could probably go wrong and I just want to
// submit!!!
struct params
{
    const Particle* p; // Particle in question
    const Particle* p_arr; // Pointer to array of particles
    int n_part; // Length of above array
    double h_fact; // Smoothing length parameter; see Price 2012 eq. 10
};



// Calculate the derivative of the weighting function with respect to h

double calc_dW_dh(const Particle &p_1, const Particle &p_2, double h) {
    /*
    // This is my attempt at calculating it, but it breaks the program even though I'm sure it's
    // analytically equivalent to the below
    // If W(r, h) = 1/h w(q) then dW(r, h)/dh = -w(q)/h^2 + 1/h * dw(q)/dh by product rule
    // dw(q)/dh = dw(q)/dq * dq/dh

    double r_ij = std::abs(p_1.pos - p_2.pos);
    double q = r_ij / h;
    double dq_dh = -r_ij / std::pow(h, 2);
    double dw_dh = dkernel_dq(q) * dq_dh;

    double dcapitalW_dh = -kernel(q)/std::pow(h, 2) + dw_dh/h;

    return dcapitalW_dh;
    */

    double r_ij = p_1.pos - p_2.pos;
    double q = std::abs(r_ij) / h;

    // PHANTOM paper eq. 16, adapted for 1D where W(r,h) = 1/h w(q) instead of 1/h^3 w(q)
    return (kernel(q) + q*dkernel_dq(q)) / (-std::pow(h, 2));
}

// Calculate the derivative of the summation with respect to h
double calc_density_dh(const Particle &p, double h, const Particle* p_arr, int n_part) {
    double d_sum = 0;
    for (int i = 0; i < n_part; i++) {
        Particle p_j = *(p_arr + i);
        double dW_dh = calc_dW_dh(p, p_j, h);
        d_sum += p_j.mass * dW_dh;
    }

    return d_sum;
}

double calc_omega(const Particle &p, ParticleArrayPtr p_arr, Config c) {
    #ifdef USE_VARIABLE_H
    // Price 2012 eq. 27
    double o_sum = calc_density_dh(p, p.h, p_arr.get(), c.n_part);

    double dh_drho = -p.h / p.density;
    o_sum *= dh_drho;
    return 1 - o_sum;
    #endif
    
    #ifndef USE_VARIABLE_H
    return 1;
    #endif
}

// Summation density calculation
double calc_density(const Particle &p, double h, const Particle* p_arr, int n_part) {
    double d_sum = 0;
    for (int i = 0; i < n_part; i++) {
        Particle p_j = *(p_arr + i);
        double q = std::abs(p.pos - p_j.pos) / h;
        double w = kernel(q);

        d_sum += p.mass * (w / h);
    }

    return d_sum;
}

// Method defining the system of density and smoothing length equations.
double smoothing_f(double x, void* params) {
    // Get parameters
    const Particle p = *((struct params*)params)->p;
    const Particle* p_arr = ((struct params*)params)->p_arr;
    int n_part = ((struct params*)params)->n_part;
    double h_fact = ((struct params*)params)->h_fact;

    // Calculate density via sum over other particles
    double density_sum = calc_density(p, x, p_arr, n_part);
    // Calculate density via expression (Price 2018 eq. 10)
    double density_exp = p.mass * h_fact / x;

    // Smoothing length equation (Price 2018 eq. 9). The aim is to optimize this to be 0.
    return density_sum - density_exp;
    
}

// Get derivative of smoothing length equation with respect to h
double smoothing_df(double x, void *params) {
    // Get parameters
    const Particle p = *((struct params*)params)->p;
    const Particle* p_arr = ((struct params*)params)->p_arr;
    int n_part = ((struct params*)params)->n_part;
    double h_fact = ((struct params*)params)->h_fact;

    // Price 2018 eq. 12
    double drho_dh_sum = calc_density_dh(p, x, p_arr, n_part);
    double drho_dh_exp = -p.mass * h_fact / std::pow(x, 2);

    return drho_dh_sum - drho_dh_exp;
}

// Function returning value and derivative. GSL algorithms want this apparently
void smoothing_fdf(double x, void *params, double *y, double *dy) {
    *y = smoothing_f(x, params);
    *dy = smoothing_df(x, params);
}

// Fallback bisection method. Not in header file since it's only called into by rootfind_h in case
// Newton's method fails
double rootfind_h_fallback(
    const Particle &p, 
    const ParticleArrayPtr p_arr,
    const Config c
) {
    int status;
    int iter = 0;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    // Since it's a fallback, and is unlikely to be used that much, set the intervals really wide
    double x = CALC_EPSILON;
    double x_lo = CALC_EPSILON; 
    double x_hi = 2*c.limit;

    struct params param = {
        &p,
        p_arr.get(),
        c.n_part,
        c.h_factor
    };

    gsl_function f = {
        &smoothing_f,
        &param
    };
    
    // Could probably use Brent to be honest, but again -- it's not going to be used unless Newton's
    // method fails, so it's probably wise to keep it simple and extremely reliable
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc(T);
    status = gsl_root_fsolver_set(s, &f, x_lo, x_hi);

    // Solver loop
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        x = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);

        status = gsl_root_test_interval(x_lo, x_hi, 0, H_EPSILON);
    } while (status == GSL_CONTINUE && iter < H_MAX_ITER_BS);

    // If this fails, then we're probably in trouble!
    if (status != GSL_SUCCESS) {
        std::cout << "[WARN] Fallback smoothing length root-finding failed for particle id " << p.id
                  << " with status '" << gsl_strerror(status) << "'" << std::endl;
    }

    gsl_root_fsolver_free(s);
    return x;

}

// Main Newton's method solver
double rootfind_h(
    const Particle &p, // p probably doesn't need to be passed by reference...oops
    const ParticleArrayPtr p_arr,
    const Config c
) {
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;

    int status;
    size_t iter = 0;

    // Initialize solver parameters
    struct params param = {
        &p,
        p_arr.get(),
        c.n_part,
        c.h_factor
    };

    gsl_function_fdf f = {
        &smoothing_f,
        &smoothing_df,
        &smoothing_fdf,
        &param
    };

    double x0, x;

    // Use previous smoothing length as a first guess for the iteration
    // That is, if it's not zero due to us currently setting up the initial smoothing lengths!
    // This generally makes Newton's method converge extremely quickly, within just a handful of
    // iterations.
    if (p.h > CALC_EPSILON) {
        x0 = p.h;
        x = p.h;
    } else {
        double mean_p_spacing = 2*c.limit / (c.n_part-1);
        x0 = c.h_factor * mean_p_spacing;
        x = c.h_factor * mean_p_spacing;
    }

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc(T);
    gsl_root_fdfsolver_set(s, &f, x);
    
    // Main solver loop
    do {
        iter++;
        status = gsl_root_fdfsolver_iterate(s);
        x0 = x;
        x = gsl_root_fdfsolver_root(s);
        status = gsl_root_test_delta(x, x0, 0, H_EPSILON);

    } while (status == GSL_CONTINUE && iter < H_MAX_ITER_NR);

    if (status != GSL_SUCCESS) {
        // Fallback to bisection
        #ifdef H_WARNINGS
        std::cout << "[WARN] Smoothing length root-finding failed for particle id " << p.id
                  << " with status '" << gsl_strerror(status) << "'" << std::endl;
        
        std::cout << "[WARN] Repeating root-finding process using bisection." << std::endl;
        #endif
        x = rootfind_h_fallback(p, p_arr, c);
    }

    gsl_root_fdfsolver_free(s);
    return x;
}

