/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * ghost_particles.hpp defines a function to initialize the ghost particles of the simulation. This
 * is separate from setup.hpp as it is necessary to redefine the ghost particles after each timestep
 * in sph.cpp as well.
 */

#ifndef ghost_particles_hpp
#define ghost_particles_hpp

#include "basictypes.hpp"

// Setup ghost particles. To be done for initial conditions and after each timestep. Takes particle
// array pointer and Config by reference, as the array must be reallocated and n_part must be
// changed.
// This assumes that smoothing lengths for each particle has already been set.
void setup_ghost_particles(ParticleArrayPtr &p_arr, Config &config);

#endif