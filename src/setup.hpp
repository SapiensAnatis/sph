/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * setup.hpp defines functions implemented in read_config.cpp. As you may have guessed, these relate
 * to reading the configuration file and setting up the program. It also defines a Config struct
 * with shared static properties.
 *
 * The configuration format is proprietary, and it does reinvent the wheel somewhat as there are
 * plenty of libraries out there that do this. That said, it is simple enough code that I prefer to
 * do it this way rather than adding an external dependency.
 */

#ifndef setup_hpp // Include guard
#define setup_hpp

#include <fstream>
#include <map>
#include <string>

struct Config {
    // Variable names correspond to parameter names in config file
    static double d_unit;
    static uint n_part; // Max: ~4 billion. No need for long in my opinion.
};

typedef std::map<std::string, std::string> ConfigMap; 

// parse_config: takes in a stream of the config file, and creates a <string, string> map of
// <propertyname, propertyvalue> to be properly parsed and stored by load_properties. A status code
// is written to the second parameter -- it's 0 if everything went okay, and 1 if something went
// wrong (e.g. parse error).
ConfigMap parse_config(std::istream &cfg_stream, int &status_code);

// load_properties: takes in a ConfigMap (from above) and initializes the values of all the needed
// properties in the Config class. This is where datatype conversion from string will be performed.
// It uses the same status code convention as above.
void load_properties(ConfigMap &config_map, int &status_code);

// Overloads of set_property. These take in a particular type of Config member by reference as well
// as a string property value, and each overload has a different way of converting the property
// value based on the type of the Config member.
void set_property(uint &prop, const std::string &prop_name, const std::string &prop_value);
void set_property(double &prop, const std::string &prop_name, const std::string &prop_value);

#endif