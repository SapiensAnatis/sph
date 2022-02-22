/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * sph.cpp implements the functions defined and explained in sph.hpp.
 */

#include <gsl/gsl_odeiv2.h>

#include "sph.hpp"

void SPHSimulation::start(double end_time) {
    // First calculate density and acceleration for all particles at t = 0
    // Consecutive for loops: acceleration calculation requires that density is define for every
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

    // And so it begins
    while (current_time < end_time) {
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

        // Half-step velocity
        p.vel += p.acc * (timestep / 2);
        // Position
        p.pos += p.vel * (timestep);

        // Recalculate density and acceleration, as position has changed
        dc(p);
        ac(p);

        // Remaining half-step velocity
        p.vel += p.acc * (timestep);
    }

    file_write();
}

void SPHSimulation::file_write() {
    outstream.open("/home/jay/Dropbox/University/Y4/PHYM004/sph/dumps/" + std::to_string(dump_counter) + ".txt");
    outstream << "# Code units: distance = " << config.d_unit << ", time = " << config.t_unit << std::endl;
    outstream << "# This file was dumped at t = " << current_time << std::endl;
    outstream << "# Columns: Particle ID / Density / Particle acceleration / Particle velocity / Particle position" << std::endl;

    for (int i = 0; i < config.n_part; i++) {
        Particle& p = p_arr[i];
        outstream << p.id << "    " << p.density << "    " << p.acc << "    " << p.vel << "    " << p.pos << std::endl;
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