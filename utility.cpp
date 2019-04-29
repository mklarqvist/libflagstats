#include <string>
#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
    std::string str;
    while (std::getline(std::cin, str)) {
        uint16_t val = std::atoi(str.c_str());
        std::cout.write((char*)&val, sizeof(uint16_t));
    }
    std::cout.flush();


    return EXIT_SUCCESS;
}