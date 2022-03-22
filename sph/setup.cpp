/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * setup.cpp implements the functions from setup.hpp -- primarily the Config class methods and
 * constructor.
 */

#include <iostream>
#include <memory>
#include <algorithm>

#include "setup.hpp"
#include "define.hpp"
#include "basictypes.hpp"
#include "ghost_particles.hpp"

#pragma region ConfigParsing

ConfigReader::ConfigReader(std::istream &config_stream) {
    ConfigMap config_map = parse_config(config_stream);
    
    set_property(config.n_part, config_map, "n_part");
    set_property(config.mass, config_map, "mass");
    set_property(config.pressure_calc, config_map, "pressure_calc");
    set_property(config.limit, config_map, "limit");
    set_property(config.v_0, config_map, "v_0");
    set_property(config.h_factor, config_map, "h_factor");
    set_property(config.t_i, config_map, "t_i");

    // 'Runtime' properties
    config.n_ghost = 0;
}

Config ConfigReader::GetConfig() {
    return config;
}

ConfigMap ConfigReader::parse_config(std::istream &cfg_stream) {
    ConfigMap result_map;

    int current_line = 0;
    std::string l;

    while (std::getline(cfg_stream, l)) {    
        current_line++;

        if (l[0] != '#' && l.length() > 0) { // Ignore comments and empty lines
            // Split string on space
            size_t space = 0;
            std::string propname;
            std::string propvalue;

            space = l.find(' ');

            if (space == std::string::npos) {
                // Error if no space found
                std::cerr << "[ERROR] Parsing error on line " << current_line << " of config file.";
                exit(1);
            }

            propname = l.substr(0, space);
            propvalue = l.substr(space + 1, l.length());

            if (propname.length() == 0 || propvalue.length() == 0) {
                // Error if there was a space but nothing on one side of it
                std::cerr << "[ERROR] Parsing error on line " << current_line << " of config file.";
                exit(1);
            }

            // What if the parameter is already in the map? Not a fatal error, but the user should
            // be told
            if (result_map.count(propname)) {
                std::cout << "[WARN] Ignoring second definition of parameter '" << propname << 
                "' on line " << current_line << " of config file." << std::endl; 
            }
            
            // .insert() will do nothing if the value is already in the map, as per the above
            result_map.insert(std::pair<std::string, std::string>(propname, propvalue));
        }
    }

    return result_map;
}

std::string ConfigReader::read_config_map(ConfigMap &config_map, const std::string &prop_name) {
    ConfigMap::iterator it = config_map.find(prop_name);
    if (it != config_map.end()) {
        // Key exists
        return it->second;
    } else {
        // The problem with this is it means every parameter is required...it should be fine for a
        // simple SPH program without optional functionality
        std::cerr << "[ERROR] Failed to find a definition for property '" << prop_name << "' in the"
        " config file." << std::endl;
        exit(1);
    }
}

void ConfigReader::set_property(int &prop, ConfigMap &config_map, const std::string &prop_name) {
    std::string prop_value = read_config_map(config_map, prop_name);
    try {
        prop = std::stoi(prop_value);
    } catch (const std::invalid_argument &ia) {
        std::cerr << "[ERROR] Failed to parse value '" << prop_value << "' for property '"
                  << prop_name << "' of type 'int'." << std::endl;
        exit(1);
    } catch (const std::out_of_range &e) {
        std::cerr << "[ERROR] The value '" << prop_value << "' is out of range for property '"
                  << prop_name << "' of type 'int'" << std::endl;
        exit(1);

    }
}

void ConfigReader::set_property(double &prop, ConfigMap &config_map, const std::string &prop_name) {
    std::string prop_value = read_config_map(config_map, prop_name);
    try {
        prop = std::stod(prop_value);
    } catch (const std::invalid_argument &ia) {
        std::cerr << "[ERROR] Failed to parse value '" << prop_value << "' for property '"
                  << prop_name << "' of type 'double'." << std::endl;
        exit(1);
    }
}

void ConfigReader::set_property(PressureCalc &prop, ConfigMap &config_map, const std::string &prop_name) {
    // Just need to get string property as int, then pass to PressureCalc enum constructor. Use
    // existing int fetch method
    int tmp_prop;
    set_property(tmp_prop, config_map, prop_name);
    prop = (PressureCalc)tmp_prop;
}

#pragma endregion
#pragma region ParticleInitialization

void init_particles(Config &config, ParticleArrayPtr &p_arr)
{
    double max_x = config.limit;
    double min_x = -max_x;
    double spacing = (max_x - min_x) / (config.n_part-1);

    // Don't put particles directly on the boundaries, as this becomes problematic
    // when trying to reflect them around the boundary to create ghost particles.
    max_x -= spacing / 2;
    min_x += spacing / 2;

    for (int i = 0; i < config.n_part; i++) {
        Particle& p  = p_arr[i];

        double spacing = (max_x - min_x) / (config.n_part-1);
        double pos = min_x + spacing * (i);
        
        // +v_0 if pos negative, -v_0 otherwise
        double vel = (pos < 0) ? config.v_0 : -config.v_0;
        
        p.pos = pos;
        p.mass = config.mass;
        // p.vel will be overwritten later if the adiabatic option is enabled, as soon as the
        // acceleration is known, which defines the pressure as an intermediate and thus the sound
        // speed
        p.vel = vel;

        // Set initial adiabatic energy
        if (config.pressure_calc == Adiabatic)
            p.u = 1/(GAMMA - 1);

        
        #ifdef USE_VARIABLE_H
        // Set variable h to a guess. Don't actually do the rootfinding, because that leads to the
        // edge particles having higher smoothing lengths and influences the ghost particle setup.
        p.h = config.h_factor * spacing;
        #endif
        
        #ifndef USE_VARIABLE_H
        // Set constant h
        p.h = CONSTANT_H;
        #endif
    }

    // In the adiabatic case, we must first calculate accelerations so that we can set the
    // initial velocitites of particles to the adiabatic sound speed, which depends on pressure.
    if (config.pressure_calc == Adiabatic) {
        auto ac = AccelerationCalculator(config, p_arr);
        for (int i = 0; i < config.n_part; i++) {
            ac(p_arr[i]);
            double c_s = ac.sound_speed(p_arr[i]);
            p_arr[i].vel = (p_arr[i].pos < 0) ? c_s : -c_s;
        }
    }
    
    // Next step: setup ghost particles. 
    setup_ghost_particles(p_arr, config);

    std::cout << "[INFO] Initialized " << config.n_ghost << " ghost particles." << std::endl;
    std::cout << "[INFO] Initialized " << config.n_part << " total particles." << std::endl;
    std::cout << "[INFO] Calculating initial conditions..." << std::endl;

    // Calculate conditions at T = 0
    auto dc2 = DensityCalculator(config, p_arr);
    auto ac2 = AccelerationCalculator(config, p_arr);
    auto ec2 = EnergyCalculator(config, p_arr);

    for (int i = 0; i < config.n_part; i++) {
        dc2(p_arr[i]);
    }
    
    // Once density is defined for all particles, can calculate derived quantities
    for (int i = 0; i < config.n_part; i++) {
        ac2(p_arr[i]);
        ec2(p_arr[i]);
    }

}

#pragma endregion


