/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * main.cpp contains the entrypoint for the program. All of the options are given in the config
 * file, due to my own personal trauma relating to invoking sph-NG setup by piping in a 20+ line
 * long text file without comments or parameter names. As a result, this entrypoint is brief and
 * contains no user input prompts.
 *
 * It mainly serves to look for a command line argument asking for a different config file, and open
 * whichever config file to a readable stream that is given to parse_config.
 *
 * It has no corresponding header file as there is no need for the functions herein to be accessed
 * elsewhere.
 */

#include <string>
#include <fstream>
#include <iostream>
#include <memory>

#include "setup.hpp"
#include "basictypes.hpp"
#include "sph.hpp"


int main(int argc, char* argv[]) {
    std::string filename;

    // Check if alternative filename argument was given
    if (argc >= 2) {
        // Second argument; first is always program name
        filename = argv[1];
    } else {
        filename = "./config.txt";
    }
    
    std::cout << "[INFO] Using config file " << filename << std::endl;

    std::ifstream config_stream;
    config_stream.open(filename);

    if (!config_stream) {
        std::cerr << "[ERROR] Could not open config file " << filename << " for reading. Perhaps"
        " the file could not be found, or you do not have permission to read it." << std::endl;
        
        return 1;
    }

    // Read in config file, and only take actual values
    auto config_reader = ConfigReader(config_stream);
    Config config = config_reader.GetConfig();

    // Allocate memory for particle array. Using a vector would've been way easier but I thought an
    // array would be mOrE eFfIcIeNt and now I can't be bothered to change all the references to
    // this type
    ParticleArrayPtr p_arr;

    try {
        p_arr = boost::make_shared<Particle[]>(config.n_part);
    } catch (std::bad_alloc &e) {
        size_t bytes = config.n_part * sizeof(Particle);

        std::cerr << "[ERROR] Failed to allocate memory for particle array!" << std::endl;
        std::cerr << "[ERROR] Attempted to allocate " << bytes << " bytes for " << config.n_part
                  << " particles" << std::endl;
        exit(1);
    }

    // Initialize position, velocity, and mass values. p_arr will be reallocated to fit the ghost
    // particles in.
    std::cout << "[INFO] Initializing particle array..." << std::endl;
    init_particles(config, p_arr);

    // Create simulation object
    auto sim = SPHSimulation(config, p_arr);
    sim.start(1);

    return 0;
}