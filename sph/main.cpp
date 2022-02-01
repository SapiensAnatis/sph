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

#include "setup.hpp"

void load_properties(ConfigMap &config_map);

int main(int argc, char* argv[]) {
    std::string filename;

    // Check if alternative filename argument was given
    if (argc >= 2) {
        // Second argument; first is always program name
        filename = argv[1];
    } else {
        filename = "config.txt";
    }
    
    std::cout << "[INFO] Using config file " << filename << std::endl;

    std::ifstream config_stream;
    config_stream.open(filename);

    if (!config_stream) {
        // The logging system is not very advanced, but hopefully using these simple tags should
        // enable users to perform ad-hoc filtering using grep or similar
        std::cerr << "[ERROR] Could not open config file " << filename << " for reading. Perhaps"
        " the file could not be found, or you do not have permission to read it." << std::endl;
        
        return 1;
    }

    ConfigMap config_map = parse_config(config_stream);
    Config config(config_map);

    std::cout << config.n_part;

    return 0;
}

// This is in main.cpp rather than setup.cpp as static variables are still only file-scoped, so
// changes to Config:: members in setup.cpp would not be reflected here.