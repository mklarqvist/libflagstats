#include <iostream>//out streams
#include <random>//random generator (c++11)
#include <chrono>//time (c++11)
#include <cassert>//assert
#include <cstring>//memset

// Reads integers line-by-line from an input stream of ASCII
// characters. This is meant to be used together with
// samtools view INPUT_FILE | cut -f 2 | utility
int main(int argc, char** argv) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator

    std::cerr << atoi(argv[1]) << std::endl;

    std::uniform_int_distribution<uint16_t> distr(0, std::numeric_limits<uint16_t>::max()-1); // right inclusive
    for (int i = 0; i < std::atoi(argv[1]); ++i) {
        uint16_t x = distr(eng);
        std::cout.write((char*)&x, sizeof(uint16_t));
    }

    return EXIT_SUCCESS;
}