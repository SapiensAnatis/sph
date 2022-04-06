/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * sph.hpp defines the SPHSimulation object, which represents a run of the simulation by containing
 * the current time, current timestep, Config, particle vector, etc. This interacts with the
 * 'calculators' module to calculate the density and acceleration, etc. at each time-step and
 * evolves each particle's position, velocity and internal energy by using an integrator.
 */

#ifndef sph_simulation_hpp
#define sph_simulation_hpp

#include <fstream>

#include "define.hpp"
#include "basictypes.hpp"
#include "calculators.hpp"

class SPHSimulation {
    public:
        // ctor
        SPHSimulation(Config c, ParticleArrayPtr p_arr) 
            : config(c), p_arr(p_arr), dc(c, p_arr), ac(c, p_arr), ec(c, p_arr), timestep(c.t_i)
        {}

        // Start the simulation (and block the thread until current_time reaches end_time)
        void start(double end_time);

    private:
        Config config;
        
        ParticleArrayPtr p_arr;

        DensityCalculator dc;
        AccelerationCalculator ac;
        EnergyCalculator ec;
        
        double current_time = 0;
        double timestep;

        int dump_counter = 0;
        std::ofstream outstream;

        // Step the simulation forward
        void step_forward();

        // Write particle information to a file: "./dumps/{dump_counter}.txt", and then increment
        // dump_counter
        void file_write();

};

#endif