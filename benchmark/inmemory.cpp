#include <iostream>
#include <chrono>
#include <memory>
#include <cstdint>
#include <random>
#include <iomanip>

#include "libalgebra.h" // pospopcnt

// force flagstats to define all implementations
#define STORM_HAVE_AVX2
#define STORM_HAVE_AVX512
#define STORM_HAVE_SSE42
#include "../libflagstats.h" // flagstats

using Clock = std::chrono::high_resolution_clock;

template <typename UNIT = std::chrono::microseconds>
Clock::time_point::rep elapsed(const Clock::time_point& t1, const Clock::time_point& t2) {
    return std::chrono::duration_cast<UNIT>(t2 - t1).count();
}

class Application {
private:
    const size_t size;
    std::unique_ptr<uint16_t[]> flags;
    std::ostream& out;
public:
    Application(size_t size)
        : size(size)
        , flags(new uint16_t[size])
        , out(std::cout) {

        initialize_input();
    }

    void run() {
        uint32_t scalar[32];
        run("scalar", FLAGSTAT_scalar, scalar);

        uint32_t sse4[32];
        const uint64_t time_sse4 = run("SSE4", FLAGSTAT_sse4, sse4);

        uint32_t sse4_improved[32];
        run("SSE4 improved", FLAGSTAT_sse4_improved, sse4_improved, sse4, time_sse4);

        uint32_t avx512[32];
        const uint64_t time_avx512 = run("AVX512", FLAGSTAT_avx512, avx512, sse4);

        uint32_t avx512_improved[32];
        const uint64_t time_avx512_imrpvoed = run("AVX512 improved", FLAGSTAT_avx512_improved, avx512_improved, sse4, time_avx512);
    }

private:
    void initialize_input() {
        std::random_device rd;
        std::mt19937 eng(rd());
        eng.seed(0); // make the results repeatable

        std::uniform_int_distribution<uint16_t> flag(0, 4096 - 1);
        for (size_t i=0; i < size; i++)
            flags[i] = flag(eng);
    }

    template <typename FUN>
    uint64_t run(const char* name,
                 FUN function,
                 uint32_t* stats,
                 uint32_t* stats_ref = nullptr,
                 uint64_t time_ref = 0)
    {
        out << "Running function " << name << ": ";
        out << std::flush;

        for (int i=0; i < 32; i++) stats[i] = 0;

        const auto t1 = Clock::now();
        function(flags.get(), size, stats);
        const auto t2 = Clock::now();
        
        const uint16_t time_us = elapsed(t1, t2);
        out << "time " << time_us << " us";
        if (time_ref != 0)
            out << " (speedup: " << double(time_ref)/time_us << ")";
        out << '\n';
        dump_stats(stats);

        if (stats_ref != nullptr) {
            const bool has_error = compare(stats_ref, stats);
        }

        return time_us;
    }

    void dump_array(uint32_t* arr, int size) {
        out << '[';
        for (int i=0; i < size; i++) {
            if (i != 0)
                out << ", ";

            out << std::setw(6);
            out << arr[i];
        }
        out << ']';
    }

    void dump_stats(uint32_t* stats) {
        out << "statistics are: ";
        out << '\n';
        for (int i=0; i < 32; i += 8) {
            out << "    ";
            dump_array(stats + i, 8);
            out << '\n';
        }
    }

    bool compare(uint32_t* reference, uint32_t* stats) {
        bool has_error = false;
        for (int i=0; i < 32; i++) {
            const uint32_t expected = reference[i];
            const uint32_t actual = stats[i];
            if (expected != actual) {
                out << "Difference at " << i << ": expected = " << expected << ", actual = " << actual << '\n';
                has_error = true;
            }
        }

        return has_error;
    }
};

int main() {
    Application app(1024 * 100);
    app.run();
}

