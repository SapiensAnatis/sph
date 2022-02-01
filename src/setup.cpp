#include "setup.hpp"
#include <iostream>

ConfigMap parse_config(std::istream &cfg_stream, int &status_code) {
    ConfigMap result_map;

    int current_line = 0;
    status_code = 0;

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
                std::cout << "[ERROR] Parsing error on line " << current_line << " of config file." << std::endl;
                status_code = 1;
            }

            propname = l.substr(0, space);
            propvalue = l.substr(space + 1, l.length());

            if (propname.length() == 0 || propvalue.length() == 0) {
                // Error if there was a space but nothing on one side of it
                std::cout << "[ERROR] Parsing error on line " << current_line <<  " of config file." << std::endl;
                status_code = 1;
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

void load_properties(ConfigMap &config_map, int &status_code) {
    // This is a bit monolithic, but I can't really think of a better way of doing it.
    set_property(Config::n_part, "n_part", config_map["n_part"]);
    set_property(Config::d_unit, "d_unit", config_map["d_unit"]);
}

void set_property(uint &prop, const std::string &prop_name, const std::string &prop_value) {
    try {
        prop = std::stoi(prop_value);
    } catch (const std::invalid_argument& ia) {
        std::cout << "[ERROR] Failed to parse value '" << ia.what() << "' for property '" <<
        prop_name << "'." << std::endl;
        exit(1);
    }
}

void set_property(double &prop, const std::string &prop_name, const std::string &prop_value) {
    try {
        prop = std::stod(prop_value);
    } catch (const std::invalid_argument& ia) {
        std::cout << "[ERROR] Failed to parse value '" << ia.what() << "' for property '" <<
        prop_name << "'." << std::endl;
        exit(1);
    }
}