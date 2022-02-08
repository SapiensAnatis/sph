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
#include <chrono>

#include "setup.hpp"
#include "density.hpp"

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
        // The logging system is not very advanced, but hopefully using these simple tags should
        // enable users to perform ad-hoc filtering using grep or similar
        std::cerr << "[ERROR] Could not open config file " << filename << " for reading. Perhaps"
        " the file could not be found, or you do not have permission to read it." << std::endl;
        
        return 1;
    }

    auto config = Config(config_stream);

    ParticleVector pv = init_particles(config);

    auto d_calc = DensityCalculator(config);
    auto outstream = std::ofstream("/home/jay/Dropbox/University/Y4/PHYM004/sph/densities.txt");
    outstream << "# Code units are " << config.d_unit << " for distance and " << config.t_unit << " for time." << std::endl;
    outstream << "# Each line is one particle." << std::endl;
    outstream << "# Column defs: position (x), density" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (Particle p : pv) {
        // Calculate density at each particle
        d_calc(p, pv);

        outstream << p.pos << "\t" << p.density << std::endl;
    }

    // I'm usually against `using namespace`s but these type names are making me reconsider...
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << duration.count() << std::endl;

    return 0;
}