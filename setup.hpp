/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * setup.hpp defines functions implemented in read_config.cpp. As you may have
 * guessed, these relate to reading the configuration file and setting up the
 * program. It also defines a static Config class, through which these values
 * are accessed while the program is running.
 *
 * The configuration format is proprietary, and it does reinvent the wheel
 * somewhat as there are plenty of libraries out there that do this. That said,
 * it is simple enough code that I prefer to do it this way rather than adding
 * an external dependency.
 */

#ifndef setup_hpp // Include guard
#define setup_hpp

#include <fstream>

namespace Config {
    // Variable names correspond to parameter names in config file
    static double d_unit;
    static uint n_part; // Max: ~4 billion. No need for long in my opinion.
};

void parse_config(std::ifstream &cfg_filestream);

#endif