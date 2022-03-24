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

    // For current quartic kernel this should be 1.25 smoothing lengths, but other kernels may have
    // this as 1 smoothing length
    double half_k_radius = KERNEL_RADIUS;

    for (int i = 0; i < config.n_part; i++) {
        Particle p = p_arr[i];
        if (p.type == Ghost)
            continue;

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
        double max_pos = config.limit - (half_k_radius * p.h);
        double min_pos = -config.limit + (half_k_radius * p.h);

        if (p.pos > max_pos || p.pos < min_pos) {
            // Add copy of p to vector
            ghost_particles.push_back(p);
        }
    }

    // Change copies so that they have inverted velocities and mirrored positions, also correct type
    for (Particle &p : ghost_particles) {
        // Add twice the vector joining the particle and boundary to the new particle, so it is
        // mirrored around the boundary.
        double vec, vel;
        if (p.pos < 0) {
            // Left border
            vec = -config.limit - p.pos;
            vel = -1;
        } else {
            // Right border
            vec = config.limit - p.pos;
            vel = 1;
        }

        p.pos += 2 * vec;
        p.vel = vel;
        p.type = Ghost;
    }

    // Reallocate particle array
    // It's a bit easier to work with raw pointers when copying with data, so we switch back to
    // shared pointers later.

    int old_n_ghost = config.n_ghost;
    int new_n_ghost = ghost_particles.size();
    int n_alive = config.n_part - config.n_ghost;

    Particle* old_ptr = p_arr.get();
    Particle* new_ptr;

    // If we need more room in the array. There isn't really any point shrinking the array if we
    // need less, as we may end up needing more in a future timestep... so this way we avoid having
    // to carry out this procedure on most timesteps.
    // The fact that we end up with empty/old memory at the end of the array doesn't matter, because
    // we loop over it using i=0; i<n_part; i++ (rather than looking for a terminator etc) so these
    // regions will be ignored.

    if (new_n_ghost > max_n_ghost) {
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
        // no longer in use (or more precisely will do so later once the Calculators no longer have
        // references to it).
        p_arr.reset(new_ptr);

        std::cout << "[INFO] Reallocated particle array to resize ghost partition from "
                  << max_n_ghost << " particles to " << new_n_ghost << " particles" << std::endl;
        
        max_n_ghost = new_n_ghost;
    }

    // Copy over new ghost particles
    std::copy(ghost_particles.begin(), ghost_particles.end(), p_arr.get() + n_alive);

    // Update config values
    config.n_ghost = new_n_ghost;
    config.n_part = n_alive + new_n_ghost;
}