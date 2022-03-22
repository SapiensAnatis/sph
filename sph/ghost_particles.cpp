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
    int n_ghost = ghost_particles.size();
    int n_alive = config.n_part - config.n_ghost;

    Particle* old_ptr = p_arr.get();
    Particle* new_ptr;

    // Reset particle counter before we allocate new ones, so we don't end up with constantly
    // increasing particle ids
    _particle_counter = 0;

    try {
        new_ptr = new Particle[n_alive + n_ghost];
    } catch (std::bad_alloc &e) {
        size_t bytes = (config.n_part + n_ghost) * sizeof(Particle);
        std::cerr << "[ERROR] Failed to reallocate array for creation of ghost particles!" << std::endl;
        std::cerr << "[ERROR] Attempted to allocate " << bytes << " bytes for " << n_alive
                  << " particles and " << n_ghost << " ghost particles" << std::endl;
        exit(1);
    }

    std::copy(old_ptr, old_ptr + n_alive, new_ptr);

    // Reinitialize shared ptr to use new_ptr; old_ptr is obsolete after copying the alive
    // particles. This will free the memory taken by old_ptr as the smart pointer detects it is no
    // longer in use.
    p_arr.reset(new_ptr);
    // Copy over new ghost particles
    std::copy(ghost_particles.begin(), ghost_particles.end(), p_arr.get() + n_alive);

    // Update config values
    config.n_ghost = n_ghost;
    config.n_part = n_alive + n_ghost;
}