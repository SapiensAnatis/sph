/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * smoothing_length.hpp defines methods that interface with the GSL to use root-finding algorithms
 * to determine the appropriate varied smoothing length. I moved these out of calculators.cpp as
 * they were starting to represent a lot of code in an already dense file.
 */

// Include guard
#ifndef smoothing_length_hpp
#define smoothing_length_hpp

#include <utility>

#include "basictypes.hpp"

// Maximum number of iterations for root-finding of smoothing length
#define H_MAX_ITER 1000
// Epsilon for root-finding of smoothing length -- used to decide when to declare success
#define H_EPSILON 1e-4
// Show root-finding steps
// #define H_DEBUG
// Show root-finding errors
#define H_ERRORS
// Use a derivative-based root-finding algorithm (HYBRIDSJ)
#define USE_DERIVATIVE


// Given a particle, the particle array, and initial guesses for smoothing length and density, use
// the GSL rootfinding capabilities to find a solution to the system of equations and return a pair
// of <smoothing length, density>.
std::pair<double, double> rootfind_h(
    const Particle &p, 
    const ParticleArrayPtr p_arr,
    const Config c
);

#endif