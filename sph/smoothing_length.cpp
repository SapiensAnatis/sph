#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
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
    double eta; // Smoothing length parameter; see Price 2010 eq. 10
};

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
int smoothing_f(const gsl_vector* x, void* params, gsl_vector* f) {
    // Get parameters
    const Particle p = *((struct params*)params)->p;
    const Particle* p_arr = ((struct params*)params)->p_arr;
    int n_part = ((struct params*)params)->n_part;
    double eta = ((struct params*)params)->eta;

    // Smoothing length
    const double x0 = gsl_vector_get(x, 0);
    // Density
    const double x1 = gsl_vector_get(x, 1);

    // Calculate density via sum over other particles
    double density = calc_density(p, x0, p_arr, n_part);

    // Smoothing length equation: h - eta(m/rho) = 0
    const double f0 = x0 - eta*(p.mass / x1);
    // Density equation: rho - summation = 0
    const double f1 = x1 - density;

    gsl_vector_set(f, 0, f0);
    gsl_vector_set(f, 1, f1);

    return GSL_SUCCESS;
}

// Calculate the derivative of the weighting function with respect to h
double calc_dw_dh(const Particle &p_1, const Particle &p_2, double h) {
    double r_ij = std::abs(p_1.pos - p_2.pos);
    double q = r_ij / h;
    double dq_dh = -r_ij / std::pow(h, 2);
    double dw_dh = dkernel_dq(q) * dq_dh;

    return dw_dh;
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

// Interface is a bit different since this isn't called in GSL-based methods but rather in a context
// with the shared pointer and Config objects etc
double calc_omega(const Particle &p, ParticleArrayPtr p_arr, Config c) {
    double o_sum = 0;
    for (int i = 0; i < c.n_part; i++) {
        Particle p_j = p_arr[i];
        double dw_dh = calc_dw_dh(p, p_j, p.h);
        o_sum += p_j.mass * dw_dh;
    }

    return 1 - o_sum;
}

#ifdef USE_DERIVATIVE
// Print the current state of the solver. Useful when wanting to see the step-by-step in case it
// produces silly values. Only used if GSL_DEBUG is defined.
void print_state (size_t iter, gsl_multiroot_fdfsolver* s)
{
    printf("iter = %3lu x = % .3f % .3f f(x) = % .3e % .3e\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1));
}

// Get Jacobian of the system of equations
int smoothing_df(const gsl_vector* x, void *params, gsl_matrix* J)
{
    // Get parameters
    const Particle p = *((struct params*)params)->p;
    const Particle* p_arr = ((struct params*)params)->p_arr;
    int n_part = ((struct params*)params)->n_part;
    double eta = ((struct params*)params)->eta;

    // Smoothing length
    const double x0 = gsl_vector_get(x, 0);
    // Density
    const double x1 = gsl_vector_get(x, 1);

    // d[h - eta*m/rho]/d[h]
    const double df00 = 1;
    // d[h - eta*m/rho]/d[rho]
    const double df01 = -(eta*p.mass)/std::pow(x1, 2);

    // d[rho - summation]/d[h]
    const double df10 = -calc_density_dh(p, x0, p_arr, n_part);
    // d[rho - summation]/d[rho]
    const double df11 = 1;

    gsl_matrix_set (J, 0, 0, df00);
    gsl_matrix_set (J, 0, 1, df01);
    gsl_matrix_set (J, 1, 0, df10);
    gsl_matrix_set (J, 1, 1, df11);

    return GSL_SUCCESS;
}

// Get function and derivative simultaneously, because GSL docs says this is good??
int smoothing_fdf(const gsl_vector* x, void *params, gsl_vector* f, gsl_matrix* J)
{
    smoothing_f(x, params, f);
    smoothing_df(x, params, J);

    return GSL_SUCCESS;
}

std::pair<double, double> rootfind_h(
    const Particle &p, 
    const ParticleArrayPtr p_arr,
    const Config c
) {
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;

    int status;
    size_t iter = 0;

    struct params param = {
        &p,
        p_arr.get(),
        c.n_part,
        c.smoothing_length
    };

    gsl_multiroot_function_fdf f = {
        &smoothing_f,
        &smoothing_df,
        &smoothing_fdf, 
        2, 
        &param
    };

    // Initial guesses for {smoothing length, density}
    // Mean particle spacing * neighbour parameter
    double h_guess = ((c.limit * 2) / c.n_part) * c.smoothing_length; 
    // Number density * mass
    double rho_guess = (c.n_part / (c.limit * 2)) * c.mass;

    // std::cout << "Initial guesses: h " << h_guess << " rho " << rho_guess << std::endl;

    double x_init[2] = {h_guess, rho_guess};
    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc(T, 2);
    gsl_multiroot_fdfsolver_set(s, &f, x);

    #ifdef H_DEBUG 
    print_state(iter, s);
    #endif

    do {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate(s);

        if (status)
            break;

        #ifdef H_DEBUG 
        print_state(iter, s); 
        #endif
        
        status = gsl_multiroot_test_residual(s->f, H_EPSILON);
    } while (status == GSL_CONTINUE && iter < H_MAX_ITER);

    #ifdef H_ERRORS
    if (status != GSL_SUCCESS) {
        std::cout << "[WARN] Smoothing length root-finding failed for particle id " << p.id << 
                     " with status '" << gsl_strerror(status) << "'" << std::endl;
    }
    #endif

    double h = gsl_vector_get(s->x, 0);
    double rho = gsl_vector_get(s->x, 1);

    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x);

    return std::pair<double, double> {h, rho};
}
#endif

#ifndef USE_DERIVATIVE
// if we're not using derivatives, every mention of gsl_multiroot_fdfsolver has to be changed to
// gsl_multiroot_fdfsolver etc. so there's a lot of code reuse involved

void print_state (size_t iter, gsl_multiroot_fsolver* s)
{
    printf("iter = %3lu x = % .3f % .3f f(x) = % .3e % .3e\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1));
}

std::pair<double, double> rootfind_h(
    const Particle &p, 
    const ParticleArrayPtr p_arr,
    const Config c
) {
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    struct params param = {
        &p,
        p_arr.get(),
        c.n_part,
        c.smoothing_length
    };

    gsl_multiroot_function f = {&smoothing_f, 2, &param};

    // Initial guesses for {smoothing length, density}
    // Mean particle spacing * neighbour parameter
    double h_guess = ((c.limit * 2) / c.n_part) * c.smoothing_length; 
    // Number density * mass
    double rho_guess = (c.n_part / (c.limit * 2)) * c.mass;

    double x_init[2] = {0.1, 1};
    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fsolver_dnewton;
    s = gsl_multiroot_fsolver_alloc(T, 2);
    gsl_multiroot_fsolver_set(s, &f, x);

    #ifdef GSL_DEBUG 
    print_state(iter, s);
    #endif

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        if (status)
            break;

        #ifdef GSL_DEBUG 
        print_state(iter, s); 
        #endif
        
        status = gsl_multiroot_test_residual(s->f, H_EPSILON);
    } while (status == GSL_CONTINUE && iter < H_MAX_ITER);

    #ifdef H_DEBUG
    if (status != GSL_SUCCESS) {
        std::cout << "[WARN] Smoothing length root-finding failed for particle id " << p.id << 
                     " with status '" << gsl_strerror(status) << "'" << std::endl;
    }
    #endif

    double h = gsl_vector_get(s->x, 0);
    double rho = gsl_vector_get(s->x, 1);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return std::pair<double, double> {h, rho};
}
#endif
