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

// Calculate 'omega' parameter from Rosswog 2009 eq. 111
// Incorporation of this quantity into the momentum equation is required when using variable
// smoothing lengths.
double calc_omega(const Particle &p, ParticleArrayPtr p_arr, Config c);

// Given a particle, the particle array, and initial guesses for smoothing length and density, use
// the GSL rootfinding capabilities to find a solution to the system of equations and return a pair
// of <smoothing length, density>.
std::pair<double, double> rootfind_h(
    const Particle &p, 
    const ParticleArrayPtr p_arr,
    const Config c
);

#endif