#include "setup.hpp"
#include <iostream>

void parse_config(std::ifstream &cfg_filestream) {
    char c;
    while (cfg_filestream) {    
        cfg_filestream.get(c);
        std::cout << c;
    }
}