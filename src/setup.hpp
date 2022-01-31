/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * setup.hpp defines functions implemented in read_config.cpp. As you may have
 * guessed, these relate to reading the configuration file and setting up the
 * program. It also defines a static Config namespace, through which these
 * values are accessed while the program is running.
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

// parse_config: takes in an infile stream of the config file, and attempts to
// initialize the values in the config namespace. A status code is written to
// the second parameter -- it's 0 if everything went okay, and 1 if something
// went wrong.
void parse_config(std::istream &cfg_stream, int &status_code);

#endif