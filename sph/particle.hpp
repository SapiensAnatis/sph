/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * particle.hpp defines the basic particle unit. These are set up in an array according to n_part.
 */

#ifndef particle_hpp // Include guard
#define particle_hpp

#include <vector>


enum ParticleType {
    Alive,
    Dead,
    Ghost
};

struct Particle {
    /* 1-dimensional values */
    double pos;
    double vel;

    double mass;
    double u; // energy

    // Ctor needs to take 0 arguments so we can initialize a vector of n Particles
    Particle()
        : pos(0), vel(0), mass(0), u(0)
    {}

    ParticleType type;
};

#endif