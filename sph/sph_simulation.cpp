/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * sph.cpp implements the functions defined and explained in sph.hpp.
 */

#include <cstdio>
#include <iostream>
#include <filesystem> // Support for this is a bit questionable, but should work with recent g++
#include <system_error>

#include "sph_simulation.hpp"
#include "ghost_particles.hpp"

void SPHSimulation::start(double end_time) {
    std::cout << "[INFO] Simulation time: " << current_time << " / " << end_time << std::endl;

    if (!std::filesystem::exists("dumps")) {
        std::error_code dir_ec;
        std::filesystem::create_directory("dumps", dir_ec);
        
        if (dir_ec.value() != 0) {
            std::cerr << "[ERROR] Failed to make directory ./dumps/ to store dump files." << std::endl;
            std::cerr << "[ERROR] Error code " << dir_ec.value() << " with message " 
                      << dir_ec.message()  << std::endl;
            std::cerr << "[HINT] You can probably get around this by just making the directory ./dumps"
                      << "manually..." << std::endl;
            exit(1);
        }
    }

    // Initial densities/acceleration/pressure etc was handled in setup.cpp
    file_write();

    // And so it begins. Note that `while(current_time < end_time)` produces
    
    // [INFO] Simulation time: 0.9 / 1
    // [INFO] Simulation time: 1 / 1
    // [INFO] Simulation time: 1.1 / 1
    
    // probably due to rounding error!
    
    #ifndef SETUP_ONLY
    while (current_time < (end_time - CALC_EPSILON)) {
        current_time += timestep;
        // These print statements help to identify where the program has had an error, should one
        // occur.
        std::cout << "[INFO] Simulation time: " << current_time << " / " << end_time << std::endl;
        step_forward();
        file_write();
    }
    #endif
}

void SPHSimulation::step_forward() {
    // Call into the integrator. For now it's a simple velocity verlet one because I remember
    // how to write that from the nbody assignment, and the GSL documentation scares me
    
    for (int i = 0; i < config.n_part; i++) {
        Particle& p = p_arr[i];

        // Don't evolve ghost particles
        if (p.type == Ghost) 
            continue;

        // Half-step velocity
        p.vel += p.acc * (timestep / 2);
        // Thermal energy
        p.u += p.du_dt * (timestep / 2);

        // Position
        p.pos += p.vel * (timestep);
    }

    /*
    // Artefacts of me debugging ghost-particle-induced segmentation faults
    std::cout << "p_arr pre-update: " << p_arr << std::endl;
    std::cout << "n_part pre-update: " << config.n_part << std::endl;
    */

    // Now that we've moved the particles, reinitialize ghost particles
    setup_ghost_particles(p_arr, config);
    // Update calculators with new n_part and possibly array pointer
    dc.update(config, p_arr);
    ac.update(config, p_arr);
    ec.update(config, p_arr);

    /*
    std::cout << "p_arr post-update: " << p_arr << std::endl;
    std::cout << "n_part post-update: " << config.n_part << std::endl;
    */

    // Perform the final half of the integration
    for (int i = 0; i < config.n_part; i++) {
        Particle& p = p_arr[i];
        if (p.type == Ghost) 
            continue;

        // Recalculate density and density-dependent quantities
        dc(p);
        ac(p);
        ec(p);

        // Remaining half-step velocity
        p.vel += p.acc * (timestep / 2);
        p.u += p.du_dt * (timestep / 2);
    }
}

void SPHSimulation::file_write() {
    // Directory should hopefully have been made in start()
    outstream.open("./dumps/" + std::to_string(dump_counter) + ".txt");

    // File header
    outstream << "# This file was dumped at t = " << current_time << std::endl;
    outstream << "# Column definitions:" << std::endl;
    outstream << "# Particle ID / Type / Smoothing length / Density / Pressure / Acceleration / Velocity / Position / Thermal energy" << std::endl;
    outstream << "# Aligned definition 'tags' for easier reading:" << std::endl;
    outstream << "# ID    TYPE     H          DENSITY  PRESS    ACCEL     VEL       POS       U" << std::endl;

    // Particle information
    for (int i = 0; i < config.n_part; i++) {
        Particle& p = p_arr[i];
        // Bit of C-style code here...
        // I want to format the strings so the floats use the same d.p. and it all lines up nicely.
        // But for some reason, no major compiler has an implementation of std::format from C++20
        // yet, and I didn't feel like adding an external dependency e.g.
        // https://github.com/fmtlib/fmt
        
        char buffer[256];
        sprintf(buffer, 
                "%4d    %s    %3.5f    %3.3f    %3.3f    %+3.3f    %+3.3f    %+3.3f    %3.3f\n", 
                p.id, ParticleTypeNames[p.type], p.h, p.density, p.pressure, p.acc, p.vel, p.pos, p.u);
        outstream << buffer;
    }

    outstream.close();

    dump_counter++;
}