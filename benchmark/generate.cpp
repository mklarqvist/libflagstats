#include <iostream>//out streams
#include <random>//random generator (c++11)
#include <chrono>//time (c++11)

// Generates N random numbers in U(0, 4096) given the first provided
// input argument.
int main(int argc, char** argv) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator

    std::uniform_int_distribution<uint16_t> distr(0, 4096-1); // right inclusive
    for (int i = 0; i < strtoull( argv[1], NULL, 10 ); ++i) {
        uint16_t x = distr(eng);
        std::cout.write((char*)&x, sizeof(uint16_t));
    }

    return EXIT_SUCCESS;
}