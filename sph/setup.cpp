#include <iostream>
#include "setup.hpp"

Config::Config(std::istream &config_stream) {
    ConfigMap config_map = this->parse_config(config_stream);
    set_property(this->n_part, config_map, "n_part");
    set_property(this->d_unit, config_map, "d_unit");
}

ConfigMap Config::parse_config(std::istream &cfg_stream) {
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

std::string Config::read_config_map(ConfigMap &config_map, const std::string &prop_name) {
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

void Config::set_property(int &prop, ConfigMap &config_map, const std::string &prop_name) {
    std::string prop_value = read_config_map(config_map, prop_name);
    try {
        prop = std::stoi(prop_value);
    } catch (const std::invalid_argument &ia) {
        std::cerr << "[ERROR] Failed to parse value '" << ia.what() << "' for property '" <<
        prop_name << "'." << std::endl;
        exit(1);
    } 
}

void Config::set_property(double &prop, ConfigMap &config_map, const std::string &prop_name) {
    std::string prop_value = read_config_map(config_map, prop_name);
    try {
        prop = std::stod(prop_value);
    } catch (const std::invalid_argument &ia) {
        std::cerr << "[ERROR] Failed to parse value '" << ia.what() << "' for property '" <<
        prop_name << "'." << std::endl;
        exit(1);
    }
}