#include "setup.hpp"
#include <iostream>

void parse_config(std::ifstream &cfg_filestream, int &status_code) {
    std::string l;
    int current_line = 0;

    while (std::getline(cfg_filestream, l)) {    
        current_line++;

        if (l[0] != '#' && l.length() > 0) { // Ignore comments and empty lines
            // Split on space
            size_t space = 0;
            std::string propname;
            std::string propvalue;

            space = l.find(' ');

            if (space == std::string::npos) {
                // Something went wrong; no space found
                std::cout << "[ERROR] Parsing error on line " << current_line << 
                " of config file." << std::endl;
                status_code = 1;
                return;
            }

            propname = l.substr(0, space);
            propvalue = l.substr(space + 1, l.length());

            if (propname.length() == 0 || propvalue.length() == 0) {
                // Same generic error. I could go to some effort to give more
                // detail, but the config syntax isn't terribly complicated;
                // hopefully a simple reference to the offending line is more
                // than enough.
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