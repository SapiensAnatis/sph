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

#include <string>
#include <fstream>
#include <map>
#include <vector>

#include "particle.hpp"

// Configuration properties that are compiled into the code
const double v_0 = 10;

typedef std::map<std::string, std::string> ConfigMap;

// Config class, used to store configuration properties. Has a constructor that takes in the
// ConfigMap and performs datatype conversipn.
class Config {
    public:
        // Variable names correspond to parameter names in config file
        int n_part;
        double d_unit;
        double t_unit;
        double mass;
        double limit;
        double v_0;
        double smoothing_length;

        Config(std::istream &config_stream);
    private:
        // parse_config: takes in a stream of the config file, and creates a <string, string> map of
        // <propertyname, propertyvalue> to be converted later in the Config constructor. 
        static ConfigMap parse_config(std::istream &cfg_stream);

        // Method to access ConfigMap and return an error if key not found. Helps to reduce code reuse in
        // set_property overloads.
        static std::string read_config_map(ConfigMap &config_map, const std::string &prop_name);

        // Overloads of set_property. These take in a particular type of Config member by reference as well
        // as a string property value, and each overload has a different way of converting the property
        // value based on the type of the Config member.
        static void set_property(int &prop, ConfigMap &config_map, const std::string &prop_name);
        static void set_property(double &prop, ConfigMap &config_map, const std::string &prop_name);
};

// Method to set up initial particle array
ParticleVector init_particles(Config c);

#endif