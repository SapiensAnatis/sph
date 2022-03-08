/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * sph.cpp implements the functions defined and explained in sph.hpp.
 */

#include <cstdio>
#include <iostream>

#include "sph.hpp"

void SPHSimulation::start(double end_time) {
    // First calculate density and acceleration for all particles at t = 0
    // Consecutive for loops: acceleration calculation requires that density is defined for every
    // other particle (otherwise div by zero!)

    // These print statements help to identify where the program has had an error, if one occurs.
    std::cout << "[INFO] Simulation time: " << current_time << " / " << end_time << std::endl;
    
    for (int i = 0; i < config.n_part; i++) {
        dc(p_arr[i]);
    }
    
    for (int i = 0; i < config.n_part; i++) {
        ac(p_arr[i]);
    }

    file_write();

    // And so it begins. Note that `while(current_time < end_time)` produces
    
    // [INFO] Simulation time: 0.9 / 1
    // [INFO] Simulation time: 1 / 1
    // [INFO] Simulation time: 1.1 / 1
    
    // probably due to rounding error!
    
    while (current_time < (end_time - CALC_EPSILON)) {
        current_time += timestep;
        std::cout << "[INFO] Simulation time: " << current_time << " / " << end_time << std::endl;
        step_forward();
    }
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
        // Position
        p.pos += p.vel * (timestep);

        // Recalculate density and acceleration, as position has changed
        dc(p);
        ac(p);

        // Remaining half-step velocity
        p.vel += p.acc * (timestep / 2);
    }

    file_write();
}

void SPHSimulation::file_write() {
    // TODO: make this less hardcoded later
    outstream.open("/home/jay/Dropbox/University/Y4/PHYM004/sph/dumps/" + std::to_string(dump_counter) + ".txt");
    outstream << "# This file was dumped at t = " << current_time << std::endl;
    outstream << "# Columns definitions:" << std::endl;
    outstream << "#ID     TYPE     DENSITY  PRESS    ACCEL     VEL       POS" << std::endl;

    for (int i = 0; i < config.n_part; i++) {
        Particle& p = p_arr[i];
        // Bit of C-style code here...
        // I want to format the strings so the floats use the same d.p. and it all lines up nicely.
        // But for some reason, no major compiler has an implementation of std::format from C++20
        // yet, and I didn't feel like adding an external dependency e.g.
        // https://github.com/fmtlib/fmt
        
        char buffer[256];
        sprintf(buffer, "%4d    %s    %3.3f    %3.3f    %+3.3f    %+3.3f    %+3.3f\n", p.id, ParticleTypeNames[p.type], p.density, p.pressure, p.acc, p.vel, p.pos);
        outstream << buffer;
    }

    outstream.close();

    dump_counter++;
}

void SPHSimulation::dump_to_stdout() {
    for (int i = 0; i < config.n_part; i++) {
        Particle& p = p_arr[i];
        std::cout << p.id << "\t" << p.density << "\t" << p.acc << "\t" << p.vel << "\t" << p.pos << std::endl;
    }
}