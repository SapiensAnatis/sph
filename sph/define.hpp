/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * define.hpp contains preprocessor directives and const variable declarations that are used as
 * compile-time configuration parameters for the program. These were previously scattered closer to
 * where they were used, but it makes more sense to have them in one centralized location.
 *
 * These are loosely split into sections by the filename that they're most relevant to, but one
 * directive may be relevant in more than one file, so it's only a guideline.
 * 
 * In general, toggles for behaviour are set using #ifdef macros in the relevant file, so to disable
 * behaviour you comment the #define out. e.g. to enable constant smoothing length, comment the line
 * containing #define USE_VARIABLE_H
 * 
 * WARNING: Always run make clean after changing anything in this file, because otherwise make will
 * reuse binaries that were compiled using old settings, which can lead to all sorts of weird stuff!
 */

#ifndef define_hpp // Include guard
#define define_hpp

// === calculators.cpp ===

// Epsilon value -- when checking if a floating point is 0, check if it's less than this instead
const double CALC_EPSILON = 1e-8;
// Use variable smoothing length
#define USE_VARIABLE_H
// If not using variable smoothing length, constant value to use
const double CONSTANT_H = 0.2;
// Polytropic index for adiabatic equation of state. Used to calculate pressures and set initial
// particle energy (Bate eq. 2.23)
const double GAMMA = 5.0/3.0;

// === smoothing_length.cpp ===

// Maximum number of iterations for root-finding of smoothing length for Newton-Raphsen
// It should converge really quickly, so this is set rather low to detect when it's going in circles
// as something is probably wrong if it takes more than 5 iterations.
#define H_MAX_ITER_NR 10
// Maximum number of iterations for root-finding of smoothing length for fallback bisection
#define H_MAX_ITER_BS 1000
// Epsilon for root-finding of smoothing length -- used to decide when to declare success
#define H_EPSILON 1e-4

// Show root-finding warnings (i.e. when fallback bisection method is used)
#define H_WARNINGS

// === sph.cpp ===

// Don't start the evolution and only generate initial conditions. Useful when debugging setup or
// root-finding.
// #define SETUP_ONLY 

#endif