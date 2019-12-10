#include <iostream>
#include <iomanip>
#include <chrono>
#include <memory>
#include <cstdint>
#include <random>
#include <vector>

// force flagstats to define all implementations
#define STORM_HAVE_AVX2
#define STORM_HAVE_AVX512
#define STORM_HAVE_SSE42

#include "libalgebra.h" // pospopcnt
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
#if defined(STORM_HAVE_CPUID)
    #if defined(__cplusplus)
    /* C++11 thread-safe singleton */
    static const int cpuid = STORM_get_cpuid();
    #else
    static int cpuid_ = -1;
    int cpuid = cpuid_;
    if (cpuid == -1) {
        cpuid = STORM_get_cpuid();

    #if defined(_MSC_VER)
        _InterlockedCompareExchange(&cpuid_, cpuid, -1);
    #else
        __sync_val_compare_and_swap(&cpuid_, -1, cpuid);
    #endif
    }
    #endif
#endif

        uint32_t scalar[32];
        run("scalar", FLAGSTAT_scalar, scalar);

#if defined(STORM_HAVE_SSE42)
        if (cpuid & STORM_CPUID_runtime_bit_SSE42) {
            uint32_t sse4[32];
            const uint64_t time_sse4 = run("SSE4", FLAGSTAT_sse4, sse4);

            uint32_t sse4_improved[32];
            run("SSE4 improved", FLAGSTAT_sse4_improved, sse4_improved, sse4, time_sse4);

            uint32_t sse4_improved2[32];
            run("SSE4 improved 2", FLAGSTAT_sse4_improved2, sse4_improved2, sse4, time_sse4);
        }
#endif

#if defined(STORM_HAVE_AVX2)
        if (cpuid & STORM_CPUID_runtime_bit_AVX2) {
            uint32_t avx2[32];
            const uint64_t time_avx2 = run("AVX2", FLAGSTAT_avx2, avx2);

            uint32_t avx2_improved[32];
            run("AVX2 improved", FLAGSTAT_avx2_improved, avx2_improved, avx2, time_avx2);
        }
#endif

#if defined(STORM_HAVE_AVX512)
        if (cpuid & STORM_CPUID_runtime_bit_AVX512BW) {
            uint32_t avx512[32];
            const uint64_t time_avx512 = run("AVX512", FLAGSTAT_avx512, avx512, scalar);

            uint32_t avx512_improved[32];
            const uint64_t time_avx512_improved = run("AVX512 improved", FLAGSTAT_avx512_improved, avx512_improved, scalar, time_avx512);

            uint32_t avx512_improved2[32];
            const uint64_t time_avx512_improved2 = run("AVX512 improved 2", FLAGSTAT_avx512_improved2, avx512_improved2, scalar, time_avx512_improved);

            uint32_t avx512_improved3[32];
            const uint64_t time_avx512_improved3 = run("AVX512 improved 3", FLAGSTAT_avx512_improved3, avx512_improved3, scalar);

            uint32_t avx512_improved4[32];
            const uint64_t time_avx512_improved4 = run("AVX512 improved 4", FLAGSTAT_avx512_improved4, avx512_improved4, scalar);
        }
#endif
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
        // test only the counters actually written by FLAGSTAT_scalar_update
        static const std::vector<int> tested_counters{
            FLAGSTAT_FQCFAIL_OFF,
            FLAGSTAT_FSECONDARY_OFF,
            FLAGSTAT_FSUPPLEMENTARY_OFF,
            FLAGSTAT_BIT12_OFF,
            FLAGSTAT_FREAD1_OFF,
            FLAGSTAT_FREAD2_OFF,
            FLAGSTAT_BIT13_OFF,
            FLAGSTAT_BIT14_OFF,
            FLAGSTAT_FUNMAP_OFF,
            FLAGSTAT_FDUP_OFF,
            FLAGSTAT_FQCFAIL_OFF + 16,
            FLAGSTAT_FSECONDARY_OFF + 16,
            FLAGSTAT_FSUPPLEMENTARY_OFF + 16,
            FLAGSTAT_BIT12_OFF + 16,
            FLAGSTAT_FREAD1_OFF + 16,
            FLAGSTAT_FREAD2_OFF + 16,
            FLAGSTAT_BIT13_OFF + 16,
            FLAGSTAT_BIT14_OFF + 16,
            FLAGSTAT_FUNMAP_OFF + 16,
            FLAGSTAT_FDUP_OFF + 16,
        };

        for (const int index: tested_counters) {
            const uint32_t expected = reference[index];
            const uint32_t actual = stats[index];
            if (expected != actual) {
                out << "Difference at " << index << ": expected = " << expected << ", actual = " << actual << '\n';
                has_error = true;
            }
        }

        return has_error;
    }
};

int main(int argc, char* argv[]) {
    const size_t default_size = 1024 * 100;

    size_t size = default_size;
    if (argc > 1)
        size = atoi(argv[1]);

    Application app(size);
    app.run();
}

