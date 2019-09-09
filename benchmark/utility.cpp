#include <string>
#include <iostream>
#include <fstream>
#include <cstring>

// Utility function accepting data from cin stream and converting
// text-based FLAG values into uint16_t words.
// Intended use:
// samtools view FILE | cut -f 2 | utility > DEST_FILE.bin
int main(int argc, char** argv) {
    std::string str;
    while (std::getline(std::cin, str)) {
        uint16_t val = std::atoi(str.c_str());
        std::cout.write((char*)&val, sizeof(uint16_t));
    }
    std::cout.flush();


    return EXIT_SUCCESS;
}