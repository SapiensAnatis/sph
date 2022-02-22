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
#include "particle.hpp"
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
    
    // std::cout << "[INFO] Using config file " << filename << std::endl;

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

    // Allocate memory for particle array. Using a vector probably would've been a whole lot easier,
    // but there's no need for all the features that vectors provide such as swapping and resizing.
    // Given what the code's actually doing, a static array is fine and is probably more efficient.

    // I'm also torturing myself by using 'smart pointers' in an effort to learn new things. These
    // will automatically be freed -- in the case of a shared pointer, once the last reference to it
    // goes out of scope. I have to use the boost version because while std::shared_ptr supports
    // array types, std::make_shared doesn't (?!). Apparently that changed in C++20, but I still
    // got compiler errors when using that standard.

    // Yes, I am now slightly wishing I would've kept everything simple and written it in C
    ParticleArrayPtr p_arr = boost::make_shared<Particle[]>(config.n_part);

    // Initialize position, velocity, and mass values
    init_particles(config, p_arr);

    // Create simulation object
    auto sim = SPHSimulation(config, p_arr);
    sim.start(1);

    return 0;
}