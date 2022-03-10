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
 */

#ifndef define_hpp // Include guard
#define define_hpp

// === smoothing_length.cpp ===

// Maximum number of iterations for root-finding of smoothing length
#define H_MAX_ITER 1000
// Epsilon for root-finding of smoothing length -- used to decide when to declare success
#define H_EPSILON 1e-4
// Show individual root-finding steps
// #define H_DEBUG
// Show root-finding errors
#define H_ERRORS
// Use a derivative-based root-finding algorithm (HYBRIDSJ)
#define USE_DERIVATIVE

// === calculators.cpp ===

// Epsilon value -- when checking if a floating point is 0, check if it's less than this instead
const double CALC_EPSILON =  1e-8;
// Use variable smoothing length
#define USE_VARIABLE_H
// If not using variable smoothing length, constant value to use
const double CONSTANT_H = 0.3;

// === setup.cpp ===
// Use a uniform, evenly-spaced distribution of particles. If commented out, the distribution is
// randomly generated
#define UNIFORM_DIST

// === sph.cpp ===
// Don't start the evolution and only generate initial conditions. Useful when debugging setup or
// root-finding.
// #define SETUP_ONLY 

#endif