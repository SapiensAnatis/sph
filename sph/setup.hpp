/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * ConfigReader.hpp defines functions implemented in read_config.cpp. As you may have guessed, these relate
 * to reading the configuration file and setting up the program. It also defines a Config struct
 * with shared static properties.
 *
 * The configuration format is proprietary, and it does reinvent the wheel somewhat as there are
 * plenty of libraries out there that do this. That said, it is simple enough code that I prefer to
 * do it this way rather than adding an external dependency.
 */

#ifndef ConfigReader_hpp // Include guard
#define ConfigReader_hpp

#include <string>
#include <fstream>
#include <map>
#include <memory>
#include <boost/smart_ptr/make_shared.hpp>

#include "particle.hpp"

// propertyname, value map read in from file
typedef std::map<std::string, std::string> ConfigMap;
// Particle array
typedef boost::shared_ptr<Particle []> ParticleArrayPtr;

enum PressureCalc {
    Isothermal,
    Adiabatic
};

// Actual configuration values. When reading the config, the program initializes a ConfigReader,
// which is composed of a Config that it initializes the values of. The program can then grab the
// Config from its ConfigReader, and avoid passing around the complete ConfigReader class.
struct Config {
    int n_part;
    double d_unit;
    double t_unit;
    double mass;
    PressureCalc pressure_calc;
    double limit;
    double v_0;
    double smoothing_length;
};

// Config class, used to store configuration properties. Has a constructor that takes in the
// ConfigMap and performs datatype conversion.
class ConfigReader {
    public:
        // Ctor reads stream and initializes Config object
        ConfigReader(std::istream &config_stream);
        
        // Return determined Config object
        Config GetConfig();
    private:
        // parse_config: takes in a stream of the config file, and creates a <string, string> map of
        // <propertyname, propertyvalue> to be converted later in the Config constructor. 
        static ConfigMap parse_config(std::istream &cfg_stream);

        // Method to access ConfigMap and return an error if key not found. Helps to reduce code
        // reuse in set_property overloads.
        static std::string read_config_map(ConfigMap &config_map, const std::string &prop_name);

        // Overloads of set_property. These take in a particular type of Config member by reference
        // as well as a string property value, and each overload has a different way of converting
        // the property value based on the type of the Config member.
        static void set_property(int &prop, ConfigMap &config_map, const std::string &prop_name);
        static void set_property(double &prop, ConfigMap &config_map, const std::string &prop_name);

        // Data
        Config config;
};

// Take in a pointer to a particle array, and loop through it to properly initialize the particles.
void init_particles(const Config &c, ParticleArrayPtr p_arr_ptr);

#endif