/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * sph.hpp defines the SPHSimulation object, which represents a run of the simulation by containing 
 * the current time, current timestep, Config, particle vector, etc. This interacts with the 
 * 'calculators' module to calculate the density and acceleration at each time-step and evolve the
 * latter by using an integrator.
 */

#ifndef sph_hpp
#define sph_hpp

#include <iostream>

#include "setup.hpp"
#include "calculators.hpp"

class SPHSimulation {
    public:
        // ctor
        SPHSimulation(const Config &c, ParticleArrayPtr p_arr) 
            : config(c), dc(c, p_arr), ac(c, p_arr), timestep(c.t_i)
        {
            // Copy pointer
            this->p_arr = ParticleArrayPtr(p_arr);
        }

        // Start the simulation (and block the thread until current_time reaches end_time)
        void start(double end_time);

    private:
        const Config &config;
        
        ParticleArrayPtr p_arr;

        DensityCalculator dc;
        AccelerationCalculator ac;
        
        double current_time = 0;
        double timestep;

        int dump_counter = 0;
        std::ofstream outstream;

        // Step the simulation forward
        void step_forward();

        // Write particle information to a file: "dumps/{dump_counter}.txt", and then increment
        // dump_counter
        void file_write();

        // Dump particle information to stdout (for debugging)
        void dump_to_stdout();

};

#endif