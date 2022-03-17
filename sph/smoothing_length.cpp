#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <iostream>

#include "smoothing_length.hpp"
#include "define.hpp"
#include "kernel.hpp"

// The code for root-finding the variable smoothing length is largely based on the Rosenbrock
// example in the GSL documentation:
// https://www.gnu.org/software/gsl/doc/html/multiroots.html#examples

// Params for root-finding method
struct params
{
    const Particle* p; // Particle in question
    const Particle* p_arr; // Pointer to array of particles
    int n_part; // Length of above array
    double eta; // Smoothing length parameter; see Price 2012 eq. 10
};



// Calculate the derivative of the weighting function with respect to h
// If W(r, h) = 1/h w(q) then dW(r, h)/dh = -1/h^2 * w(q) + 1/h * dw(q)/dh by product rule
// dw(q)/dh = dw(q)/dq * dq/dh
double calc_dw_dh(const Particle &p_1, const Particle &p_2, double h) {
    double r_ij = std::abs(p_1.pos - p_2.pos);
    double q = r_ij / h;
    double dq_dh = -r_ij / std::pow(h, 2);
    double dw_dh = dkernel_dq(q) * dq_dh;

    double dcapitalW_dh = -1/std::pow(h, 2) * kernel(q) + 1/h * dw_dh;

    return dcapitalW_dh;
}

// Calculate the derivative of the summation with respect to h
double calc_density_dh(const Particle &p, double h, const Particle* p_arr, int n_part) {
    double d_sum = 0;
    for (int i = 0; i < n_part; i++) {
        Particle p_j = *(p_arr + i);
        double dw_dh = calc_dw_dh(p, p_j, h);
        d_sum -= p_j.mass * dw_dh / std::pow(h, 2);
    }

    return d_sum;
}

double calc_omega(const Particle &p, ParticleArrayPtr p_arr, Config c) {
    #ifdef USE_VARIABLE_H
    double o_sum = 0;
    for (int i = 0; i < c.n_part; i++) {
        Particle p_j = p_arr[i];
        double dw_dh = calc_dw_dh(p, p_j, p.h);
        o_sum += p_j.mass * dw_dh;
    }

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
    double eta = ((struct params*)params)->eta;

    // Calculate density via sum over other particles
    double density = calc_density(p, x, p_arr, n_part);

    // Smoothing length equation: h - eta(m/rho) = 0
    return x - eta*(p.mass / density);
}

// Get derivative of smoothing length equation with respect to h
double smoothing_df(double x, void *params) {
    // Get parameters
    const Particle p = *((struct params*)params)->p;
    const Particle* p_arr = ((struct params*)params)->p_arr;
    int n_part = ((struct params*)params)->n_part;
    double eta = ((struct params*)params)->eta;

    double drho_dh = calc_density_dh(p, x, p_arr, n_part);

    // d[h - eta(m/rho)]/d[h]
    return 1 + eta * p.mass / std::pow(drho_dh, 2);
}

void smoothing_fdf(double x, void *params, double *y, double *dy) {
    *y = smoothing_f(x, params);
    *dy = smoothing_df(x, params);
}

double rootfind_h_fallback(
    const Particle &p, 
    const ParticleArrayPtr p_arr,
    const Config c
) {
    int status;
    int iter = 0;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double x = 0;
    double x_lo = CALC_EPSILON; 
    double x_hi = c.limit;

    struct params param = {
        &p,
        p_arr.get(),
        c.n_part,
        c.smoothing_length
    };

    gsl_function f = {
        &smoothing_f,
        &param
    };

    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &f, x_lo, x_hi);

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        x = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);

        status = gsl_root_test_interval(x_lo, x_hi, 0, H_EPSILON);
    } while (status == GSL_CONTINUE && iter < H_MAX_ITER);

    if (status != GSL_SUCCESS) {
        std::cout << "[WARN] Fallback smoothing length root-finding failed for particle id " << p.id
                  << " with status '" << gsl_strerror(status) << "'" << std::endl;
    }

    gsl_root_fsolver_free(s);
    return x;

}

double rootfind_h(
    const Particle &p, 
    const ParticleArrayPtr p_arr,
    const Config c
) {
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;

    int status;
    size_t iter = 0;

    struct params param = {
        &p,
        p_arr.get(),
        c.n_part,
        c.smoothing_length
    };

    gsl_function_fdf f = {
        &smoothing_f,
        &smoothing_df,
        &smoothing_fdf,
        &param
    };

    double x0, x = 0.1;
    
    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc(T);
    gsl_root_fdfsolver_set(s, &f, x);
    
    do {
        iter++;
        status = gsl_root_fdfsolver_iterate(s);
        x0 = x;
        x = gsl_root_fdfsolver_root(s);
        status = gsl_root_test_delta(x, x0, 0, H_EPSILON);

    } while (status == GSL_CONTINUE && iter < H_MAX_ITER);

    if (status != GSL_SUCCESS) {
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

