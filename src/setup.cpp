#include "setup.hpp"
#include <iostream>

void parse_config(std::istream &cfg_stream, int &status_code) {
    std::string l;
    int current_line = 0;

    while (std::getline(cfg_stream, l)) {    
        current_line++;

        if (l[0] != '#' && l.length() > 0) { // Ignore comments and empty lines
            // Split on space
            size_t space = 0;
            std::string propname;
            std::string propvalue;

            space = l.find(' ');

            if (space == std::string::npos) {
                // Error if no space found
                std::cout << "[ERROR] Parsing error on line " << current_line << 
                " of config file." << std::endl;
                status_code = 1;
                return;
            }

            propname = l.substr(0, space);
            propvalue = l.substr(space + 1, l.length());

            if (propname.length() == 0 || propvalue.length() == 0) {
                // Error if there was a space but nothing on one side of it
                std::cout << "[ERROR] Parsing error on line " << current_line << 
                " of config file." << std::endl;
                status_code = 1;
                return;
            }
            
            std::cout << propname << " is " << propvalue << std::endl;
        }
    }
    
    status_code = 0;
}