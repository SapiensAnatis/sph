/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * ghost_particles.cpp implements the function from ghost_particles.hpp.
 */

#include <vector>
#include <memory>
#include <iostream>

#include "ghost_particles.hpp"
#include "define.hpp"
#include "calculators.hpp"
#include "kernel.hpp"

// How much room is there at the end of the array for ghost particles?
static int max_n_ghost = 0;

void setup_ghost_particles(ParticleArrayPtr &p_arr, Config &config) {
    // Collect particles near the left and right boundary
    std::vector<Particle> ghost_particles;

    // For current quartic kernel this should be 2.5 smoothing lengths. Defined in kernel.hpp.
    
    // In theory you should only check for particles within half the radius of the boundary, since
    // mirroring a particle that is farther away than that will not affect said particle. However,
    // I found that particles near the edge needed their neigbours as well as themselves replicated
    // to preserve the boundary conditions.
    double radius = KERNEL_RADIUS;
    
    for (int i = 0; i < config.n_part; i++) {
        Particle p = p_arr[i];
        if (p.type == Ghost)
            continue;

        // Sanity check; smoothing length may be uninitialized
        if (p.h < CALC_EPSILON) {
            std::cerr << "[ERROR] Error in ghost particle initialization: Particle id " << p.id 
                      << " has smoothing length " << p.h << std::endl;
            throw std::logic_error(
                "Ghost particle initialization: Particle had zero smoothing length!"
            );
        }

        // If a particle's distance to either boundary is greater than half the truncation radius of
        // the kernel, then it would not be affected by a mirrored ghost particle if were one to be
        // created.
        double max_pos = config.limit - (radius * p.h);
        double min_pos = -config.limit + (radius * p.h);

        if (p.pos > max_pos || p.pos < min_pos) {
            // Add copy of p to vector
            ghost_particles.push_back(p);
        }
    }

    // Change copies so that they have inverted velocities, mirrored positions, and correct type
    for (Particle &p : ghost_particles) {
        double vec;
        if (p.pos < 0) {
            // Left border
            vec = -config.limit - p.pos;
        } else {
            // Right border
            vec = config.limit - p.pos;
        }

        // Mirror around boundary by adding 2*(vector joining particle and boundary) to its position
        p.pos += 2 * vec;
        p.vel *= -1;
        p.type = Ghost;
    }

    int new_n_ghost = ghost_particles.size();
    int n_alive = config.n_part - config.n_ghost;

    if (new_n_ghost > max_n_ghost) {
        /* 
         * If we need more room in the array. There isn't really any point shrinking the array if we
         * need less, as we may end up needing more in a future timestep... so this way we avoid
         * having to carry out this procedure on most timesteps.
         */
        Particle* old_ptr = p_arr.get();
        Particle* new_ptr;

        try {
            new_ptr = new Particle[n_alive + new_n_ghost];
        } catch (std::bad_alloc &e) {
            size_t bytes = (config.n_part + new_n_ghost) * sizeof(Particle);
            std::cerr << "[ERROR] Failed to reallocate array for creation of ghost particles!" << std::endl;
            std::cerr << "[ERROR] Attempted to allocate " << bytes << " bytes for " << n_alive
                      << " particles and " << new_n_ghost << " ghost particles" << std::endl;
            exit(1);
        }

        // Copy over existing alive particles
        std::copy(old_ptr, old_ptr + n_alive, new_ptr);

        // Reinitialize shared ptr to use new_ptr; old_ptr is obsolete after copying the alive
        // particles. This will free the memory taken by old_ptr as the smart pointer detects it is
        // no longer in use
        p_arr.reset(new_ptr);

        std::cout << "[INFO] Array reallocated to resize ghost partition from " << max_n_ghost 
                  << " to " << new_n_ghost << " particles." << std::endl;
        
        max_n_ghost = new_n_ghost;
    }

    // Copy over new ghost particles
    std::copy(ghost_particles.begin(), ghost_particles.end(), p_arr.get() + n_alive);

    // Update config values
    config.n_ghost = new_n_ghost;
    config.n_part = n_alive + new_n_ghost;
}