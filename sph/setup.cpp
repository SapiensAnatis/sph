/* 
 * PHYM004 Project 2 / Jay Malhotra
 *
 * setup.cpp implements the functions from setup.hpp -- primarily the Config class methods and
 * constructor.
 */

#include <iostream>
#include <random>
#include <memory>
#include <algorithm>
#include <string.h>

#include "define.hpp"
#include "setup.hpp"

#pragma region ConfigParsing

ConfigReader::ConfigReader(std::istream &config_stream) {
    ConfigMap config_map = parse_config(config_stream);
    
    set_property(config.n_part, config_map, "n_part");
    set_property(config.mass, config_map, "mass");
    set_property(config.pressure_calc, config_map, "pressure_calc");
    set_property(config.limit, config_map, "limit");
    set_property(config.v_0, config_map, "v_0");
    set_property(config.smoothing_length, config_map, "smoothing_length");
    set_property(config.t_i, config_map, "t_i");

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

    #ifndef UNIFORM_DIST
        // Ensure randomness of position generation
        auto eng = std::default_random_engine(std::random_device{}());
        auto rand = std::uniform_real_distribution<double>(min_x, max_x);
    #endif

    for (int i = 0; i < config.n_part; i++) {
        Particle& p  = p_arr[i];

        #ifndef UNIFORM_DIST
            double pos = rand(eng);
        #endif

        #ifdef UNIFORM_DIST
            double spacing = (max_x - min_x) / (config.n_part-1);
            double pos = min_x + spacing * i;
        #endif
        
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

        // Set constant h
        #ifndef USE_VARIABLE_H
        p.h = CONSTANT_H;
        #endif
    }

    
    init_ghost_particles(config, p_arr);

    // Calculate conditions at T = 0. We are doing this here rather than in SPHSimulation as the
    // adiabatic test requires initial velocities are set to sound speed, which requires knowing
    // the pressures and densities, which requires calculating initial acceleration.

    // Create a new DensityCalculator with new array reference from ghost particle reallocation
    auto dc = DensityCalculator(config, p_arr);
    auto ac = AccelerationCalculator(config, p_arr);
    auto ec = EnergyCalculator(config, p_arr);

    for (int i = 0; i < config.n_part; i++) {
        dc(p_arr[i]);
    }

    for (int i = 0; i < config.n_part; i++) {
        ac(p_arr[i]);
    }

    // Now that pressures are defined (by AccelerationCalculator) we can overwrite initial
    // velocities for the adiabatic test
    if (config.pressure_calc == Adiabatic) {
        for (int i = 0; i < config.n_part; i++) {
            double c_s = ac.sound_speed(p_arr[i]);
            p_arr[i].vel = (p_arr[i].pos < 0) ? c_s : -c_s;
        }
    }

    std::cout << "[INFO] Initialized " << config.n_ghost << " ghost particles." << std::endl;
    std::cout << "[INFO] Initialized " << config.n_part << " total particles." << std::endl;
}

void init_ghost_particles(Config &config, ParticleArrayPtr &p_arr) {
    // Initialize ghost particles. This is not the 'proper' way of doing it, which is based on
    // neighbour trees etc., but essentially the way it works is: For the leftmost and rightmost
    // particle in the array, collect all the particles within a smoothing length and duplicate
    // them, then rotate them 180 about the selected particle, placing them outside of the boundary.
    // They are then given opposite velocities to the selected particle to stop it from going out
    // of bounds.

    // Seperate function as the scope of init_particles is cluttered with variable names I want to
    // use and it's big enough already

    // Search for boundary particles. If UNIFORM_DIST is defined, we could just immediately pick
    // the indices 0 and n_part, but let's be universal for the sake of it!
    double min_x = 0;
    double min_x_idx;

    double max_x = 0;
    double max_x_idx;

    // Just going through keeping track of what the lowest/highest positions seen so far is, and
    // their corresponding indices
    for (int i = 0; i < config.n_part; i++) {
        Particle& p = p_arr[i];
        
        if (p.pos < min_x) {
            min_x_idx = i;
            min_x = p.pos;
        } else if (p.pos > max_x) {
            max_x_idx = i;
            max_x = p.pos;
        }
    }

    Particle& left = p_arr[min_x_idx];
    Particle& right = p_arr[max_x_idx];

    // Collect neighbours of left and right
    std::vector<Particle> l_neighbours; 
    std::vector<Particle> r_neighbours;

    #ifdef USE_VARIABLE_H
    // To know what particles are the neighbours now that we have variable smoothing length, we
    // first need to calculate the density of left and right, since that will set their smoothing
    // lengths.
    auto dc = DensityCalculator(config, p_arr);
    dc(left);
    dc(right);
    #endif

    for (int i = 0; i < config.n_part; i++) {
        Particle p = p_arr[i];

        if (p == left || p == right)
            // Don't collect the particles themselves
            continue;

        double l_dist = std::abs(p.pos - left.pos);
        double r_dist = std::abs(p.pos - right.pos);

        // Creates copies of p, since the vector is of non-reference Particle. 2*smoothing length
        // as that is the limit of the kernel, so ensures that edge particles have full neighbours
        if (l_dist < left.h * 2)
            l_neighbours.push_back(p);
        else if (r_dist < right.h * 2)
            r_neighbours.push_back(p);   
    }



    // Set up neighbours (mirror around edges & set velocity and type)
    for (Particle &p : l_neighbours) {
        // Add twice the vector joining the neighbour and origin particle to the neighbour
        // particle's position
        double vec = left.pos - p.pos;
        p.pos += 2*vec;
        p.vel *= -1;
        p.type = Ghost;
    }

    for (Particle &p : r_neighbours) {
        double vec = right.pos - p.pos;
        p.pos += 2*vec;
        p.vel *= -1;
        p.type = Ghost;
    }

    // Before adding to array, dynamically reallocate. The array is already full with n_part, but we
    // can now grow it because we know how many ghost particles are about to be added

    // It's a bit easier to work with raw pointers when copying with data, so we switch back to
    // shared pointers later.
    int n_ghost = l_neighbours.size() + r_neighbours.size();
    config.n_ghost = n_ghost;

    Particle* old_ptr = p_arr.get();
    Particle* new_ptr;

    try {
        new_ptr = new Particle[config.n_part + n_ghost];
    } catch (std::bad_alloc &e) {
        size_t bytes = (config.n_part + n_ghost) * sizeof(Particle);
        std::cerr << "[ERROR] Failed to reallocate array for creation of ghost particles!" << std::endl;
        std::cerr << "[ERROR] Attempted to allocate " << bytes << " bytes for " << config.n_part
                  << " particles and " << n_ghost << " ghost particles" << std::endl;
        exit(1);
    }

    // Copy over old data
    std::copy(old_ptr, old_ptr + config.n_part, new_ptr);

    // Reinitialize shared ptr. This will free old_ptr as the smart pointer detects it is no longer
    // in use.
    p_arr.reset(new_ptr);
    // Add to end of array
    std::copy(l_neighbours.begin(), l_neighbours.end(), p_arr.get() + config.n_part);

    // Imperative to update n_part, both for copy below and so that ghost particles are not
    // ignored in subsequent iteration e.g. calculators.cpp.
    config.n_part += l_neighbours.size();

    std::copy(r_neighbours.begin(), r_neighbours.end(), p_arr.get() + config.n_part);

    config.n_part += r_neighbours.size();
}

#pragma endregion


