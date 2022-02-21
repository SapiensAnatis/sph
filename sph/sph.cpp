/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * sph.cpp implements the functions defined and explained in sph.hpp.
 */

#include "sph.hpp"
#include <gsl/gsl_odeiv2.h>

void SPHSimulation::start(double end_time) {
    // First calculate density and acceleration for all particles at t = 0
    for (int i = 0; i < this->config.n_part; i++) {
        Particle& p = this->p_arr[i];

        // These are Calculator objects which operate on a particle and set some specific property.
        // See calculators.hpp/cpp for more info.
        this->dc(p, p_arr);
        this->ac(p, p_arr);
    }

    while (this->current_time < end_time) {
        this->step_forward();
    }
}

void SPHSimulation::step_forward() {
    // Call into the integrator. For now it's a simple velocity verlet one because I remember
    // how to write that from the nbody assignment, and the GSL documentation scares me
    for (int i = 0; i < this->config.n_part; i++) {
        Particle& p = this->p_arr[i];

        // Half-step velocity
        p.vel += p.acc * (this->timestep / 2);
        // Position
        p.pos += p.vel * (this->timestep);

        // Recalculate density and acceleration, as position has changed
        this->dc(p, p_arr);
        this->ac(p, p_arr);

        // Remaining half-step velocity
        p.vel += p.acc * (this->timestep);
    }

    this->current_time += timestep;
    this->step_counter++;
}