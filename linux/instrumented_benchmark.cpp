#ifdef __linux__

/* ****************************
*  Definitions
******************************/
#include <memory>
#include <cassert>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <libgen.h>
#include <random>
#include <string>
#include <vector>
#include <chrono>

#include "libflagstats.h"
#include "linux-perf-events.h"
#include "aligned_alloc.h"

#ifdef ALIGN
#include "memalloc.h"
#   define memory_allocate(size) aligned_alloc(64, (size))
#else
#   define memory_allocate(size) malloc(size)
#endif

// Definition for microsecond timer.
typedef std::chrono::high_resolution_clock::time_point clockdef;

// FLAGSTATS_func methods[] = {FLAGSTAT_sse4, FLAGSTAT_avx2, FLAGSTAT_avx512};

typedef struct {
    long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
    long long n_sgltn[2], n_read1[2], n_read2[2];
    long long n_dup[2];
    long long n_diffchr[2], n_diffhigh[2];
    long long n_secondary[2], n_supp[2];
} bam_flagstat_t;

#define SAMTOOLS_flagstat_loop(s, c) do {                   \
    int w = (c & FLAGSTAT_FQCFAIL)? 1 : 0;                  \
    ++(s)->n_reads[w];                                      \
    if (c & FLAGSTAT_FSECONDARY ) {                         \
        ++(s)->n_secondary[w];                              \
    } else if (c & FLAGSTAT_FSUPPLEMENTARY ) {              \
        ++(s)->n_supp[w];                                   \
    } else if (c & FLAGSTAT_FPAIRED) {                      \
        ++(s)->n_pair_all[w];                               \
        if ( (c & FLAGSTAT_FPROPER_PAIR) && !(c & FLAGSTAT_FUNMAP) ) ++(s)->n_pair_good[w]; \
        if (c & FLAGSTAT_FREAD1) ++(s)->n_read1[w];         \
        if (c & FLAGSTAT_FREAD2) ++(s)->n_read2[w];         \
        if ((c & FLAGSTAT_FMUNMAP) && !(c & FLAGSTAT_FUNMAP)) ++(s)->n_sgltn[w]; \
        if (!(c & FLAGSTAT_FUNMAP) && !(c & FLAGSTAT_FMUNMAP)) {      \
            ++(s)->n_pair_map[w];                           \
        }                                                   \
    }                                                       \
    if (!(c & FLAGSTAT_FUNMAP)) ++(s)->n_mapped[w];         \
    if (c & FLAGSTAT_FDUP) ++(s)->n_dup[w];                 \
} while (0)

void SAMTOOLS_func(bam_flagstat_t* s, const uint16_t c) {
    int w = (c & FLAGSTAT_FQCFAIL) ? 1 : 0; 
    ++(s)->n_reads[w];
    if (c & FLAGSTAT_FSECONDARY ) {
        ++(s)->n_secondary[w];
    } else if (c & FLAGSTAT_FSUPPLEMENTARY ) {
        ++(s)->n_supp[w];
    } else if (c & FLAGSTAT_FPAIRED) {
        ++(s)->n_pair_all[w];
        if ( (c & FLAGSTAT_FPROPER_PAIR) && !(c & FLAGSTAT_FUNMAP) ) ++(s)->n_pair_good[w];
        if (c & FLAGSTAT_FREAD1) ++(s)->n_read1[w];
        if (c & FLAGSTAT_FREAD2) ++(s)->n_read2[w];
        if ((c & FLAGSTAT_FMUNMAP) && !(c & FLAGSTAT_FUNMAP)) ++(s)->n_sgltn[w];
        if (!(c & FLAGSTAT_FUNMAP) && !(c & FLAGSTAT_FMUNMAP)) {
            ++(s)->n_pair_map[w];
        }
    }
    if (!(c & FLAGSTAT_FUNMAP)) ++(s)->n_mapped[w];
    if (c & FLAGSTAT_FDUP) ++(s)->n_dup[w];
}

#ifdef __GNUC__
__attribute__((optimize("no-tree-vectorize")))
#endif
int samtools_flagstats(const uint16_t* array, uint32_t len, uint32_t* flags) {
    bam_flagstat_t* s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));
    for (int i = 0; i < len; ++i) {
        // SAMTOOLS_flagstat_loop(s, array[i]);
        SAMTOOLS_func(s, array[i]);
    }
    flags[0] = s->n_read1[0]; // prevent optimzie away
    flags[1] += flags[0];
    free(s);
    return 0;
}


void print16(uint32_t *flags) {
    for (int k = 0; k < 16; k++)
        printf(" %8u ", flags[k]);
    printf("\n");
}

std::vector<unsigned long long>
compute_mins(std::vector< std::vector<unsigned long long> > allresults) {
    if (allresults.size() == 0)
        return std::vector<unsigned long long>();
    
    std::vector<unsigned long long> answer = allresults[0];
    
    for (size_t k = 1; k < allresults.size(); k++) {
        assert(allresults[k].size() == answer.size());
        for (size_t z = 0; z < answer.size(); z++) {
            if (allresults[k][z] < answer[z])
                answer[z] = allresults[k][z];
        }
    }
    return answer;
}

std::vector<double>
compute_averages(std::vector< std::vector<unsigned long long> > allresults) {
    if (allresults.size() == 0)
        return std::vector<double>();
    
    std::vector<double> answer(allresults[0].size());
    
    for (size_t k = 0; k < allresults.size(); k++) {
        assert(allresults[k].size() == answer.size());
        for (size_t z = 0; z < answer.size(); z++) {
            answer[z] += allresults[k][z];
        }
    }

    for (size_t z = 0; z < answer.size(); z++) {
        answer[z] /= allresults.size();
    }
    return answer;
}

/**
 * @brief 
 * 
 * @param n          Number of integers.
 * @param iterations Number of iterations.
 * @param fn         Target function pointer.
 * @param verbose    Flag enabling verbose output.
 * @return           Returns true if the results are correct. Returns false if the results
 *                   are either incorrect or the target function is not supported.
 */
bool benchmark(uint32_t n, uint32_t iterations, FLAGSTATS_func fn, bool verbose, bool test) {
    std::vector<int> evts;
    uint16_t* vdata = (uint16_t*)memory_allocate(n * sizeof(uint16_t));
    std::unique_ptr<uint16_t, decltype(&free)> dataholder(vdata, free);
    if(verbose) {
      printf("alignment: %d\n", get_alignment(vdata));
    }
    evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
    evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
    evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
    evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
    evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
    evts.push_back(PERF_COUNT_HW_REF_CPU_CYCLES);
    LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
    std::vector<unsigned long long> results; // tmp buffer
    std::vector< std::vector<unsigned long long> > allresults;
    results.resize(evts.size());
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 0xFFFF);

    bool isok = true;
    for (uint32_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < n; k++) {
            vdata[k] = dis(gen); // random init.
        }
        uint32_t correctflags[16] = {0};
        FLAGSTAT_scalar(vdata, n, correctflags); // this is our gold standard
        uint32_t flags[16] = {0};
        
        unified.start();
        fn(vdata, n, flags);
        unified.end(results);

        uint64_t tot_obs = 0;
        for (size_t k = 0; k < 16; ++k) tot_obs += flags[k];
        if (tot_obs == 0) { // when a method is not supported it returns all zero
            return false;
        }

        for (size_t k = 0; k < 16; k++) {
            if (correctflags[k] != flags[k]) {
                if (test) {
                    printf("bug:\n");
                    printf("expected : ");
                    print16(correctflags);
                    printf("got      : ");
                    print16(flags);
                    return false;
                } else {
                    isok = false;
                }
            }
        }
        allresults.push_back(results);
    }

    std::vector<unsigned long long> mins = compute_mins(allresults);
    std::vector<double> avg = compute_averages(allresults);
    
    if (verbose) {
        printf("instructions per cycle %4.2f, cycles per 16-bit word:  %4.3f, "
               "instructions per 16-bit word %4.3f \n",
                double(mins[1]) / mins[0], double(mins[0]) / n, double(mins[1]) / n);
        // first we display mins
        printf("min: %8llu cycles, %8llu instructions, \t%8llu branch mis., %8llu "
               "cache ref., %8llu cache mis.\n",
                mins[0], mins[1], mins[2], mins[3], mins[4]);
        printf("avg: %8.1f cycles, %8.1f instructions, \t%8.1f branch mis., %8.1f "
               "cache ref., %8.1f cache mis.\n",
                avg[0], avg[1], avg[2], avg[3], avg[4]);
    } else {
        printf("cycles per 16-bit word:  %4.3f; ref cycles per 16-bit word: %4.3f \n", double(mins[0]) / n, double(mins[5]) / n);
    }

    return isok;
}

/**
 * @brief 
 * 
 * @param n          Number of integers.
 * @param m          Number of arrays.
 * @param iterations Number of iterations.
 * @param fn         Target function pointer.
 * @param verbose    Flag enabling verbose output.
 * @return           Returns true if the results are correct. Returns false if the results
 *                   are either incorrect or the target function is not supported.
 */
bool benchmarkMany(const std::string& fn_name, uint32_t n, uint32_t m, uint32_t iterations, FLAGSTATS_func fn, bool verbose, bool test, bool tabular) {
    std::vector<int> evts;
#ifdef ALIGN
    std::vector<std::vector<uint16_t,AlignedSTLAllocator<uint16_t,64>>> vdata(m, std::vector<uint16_t,AlignedSTLAllocator<uint16_t,64>>(n));
#else
    std::vector<std::vector<uint16_t>> vdata(m, std::vector<uint16_t>(n));
#endif
#ifdef ALIGN
    for(auto & x : vdata) {
      assert(get_alignment(x.data()) == 64);
    }
#endif
    if(verbose && !tabular) {
      printf("alignments: ");
      for(auto & x : vdata) {
        printf("%d ", get_alignment(x.data()));
      }
      printf("\n");
    }    
    evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
    evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
    evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
    evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
    evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
    evts.push_back(PERF_COUNT_HW_REF_CPU_CYCLES);
    LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
    std::vector<unsigned long long> results; // tmp buffer
    std::vector< std::vector<unsigned long long> > allresults;
    std::vector<uint32_t> times;
    results.resize(evts.size());
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 0xFFFF);

    bool isok = true;
    for (uint32_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < vdata.size(); k++) {
            for(size_t k2 = 0; k2 < vdata[k].size() ; k2++) { 
               vdata[k][k2] = dis(gen); // random init.
            }
        }

        std::vector<std::vector<uint32_t>> flags(m,std::vector<uint32_t>(16*2));

        const clockdef t1 = std::chrono::high_resolution_clock::now();
        unified.start();
        for (size_t k = 0; k < m ; k++) {
            fn(vdata[k].data(), vdata[k].size(), flags[k].data());
        }
        unified.end(results);
        const clockdef t2 = std::chrono::high_resolution_clock::now();
        allresults.push_back(results);

        const auto time_span = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
        times.push_back(time_span.count());
    }

    uint32_t tot_time = std::accumulate(times.begin(), times.end(), 0);
    double mean_time = tot_time / times.size();

    std::vector<unsigned long long> mins = compute_mins(allresults);
    std::vector<double> avg = compute_averages(allresults);

    double throughput = ((2*n) / (1024*1024.0)) / (mean_time / 1000000000.0);
    
    if (tabular) {
        for (int i = 0; i < iterations; ++i) {
            throughput = ((2*n) / (1024*1024.0)) / (times[i] / 1000000000.0);
            printf("%s\t%u\t%d\t", fn_name.c_str(), n, i);
            printf("%4.2f\t%4.3f\t%4.3f\t",
                    double(allresults[i][1]) / allresults[i][0], double(allresults[i][0]) / (n*m), double(allresults[i][1]) / (n*m));
            printf("%llu\t%llu\t%llu\t%llu\t%llu\t",
                    allresults[i][0], allresults[i][1], allresults[i][2], allresults[i][3], allresults[i][4]);
            printf("%u\t%4.2f\n", times[i], throughput);
        }
    } else if (verbose) {
        printf("instructions per cycle %4.2f, cycles per 16-bit word:  %4.3f, "
               "instructions per 16-bit word %4.3f \n",
                double(mins[1]) / mins[0], double(mins[0]) / (n*m), double(mins[1]) / (n*m));
        // first we display mins
        printf("min: %8llu cycles, %8llu instructions, \t%8llu branch mis., %8llu "
               "cache ref., %8llu cache mis.\n",
                mins[0], mins[1], mins[2], mins[3], mins[4]);
        printf("avg: %8.1f cycles, %8.1f instructions, \t%8.1f branch mis., %8.1f "
               "cache ref., %8.1f cache mis.\n",
                avg[0], avg[1], avg[2], avg[3], avg[4]);
        printf("avg time: %f ns, %4.2f mb/s\n", mean_time, throughput);
    } else {
        printf("cycles per 16-bit word:  %4.3f; ref cycles per 16-bit word: %4.3f \n", double(mins[0]) / (n*m), double(mins[5]) / (n*m));
    }

    return isok;
}

bool benchmarkManyMemoryOptimized(const std::string& fn_name, uint32_t n, uint32_t m, uint32_t iterations, FLAGSTATS_func fn, bool verbose, bool test, bool tabular) {
    std::vector<int> evts;
// #ifdef ALIGN
//     std::vector<std::vector<uint16_t,AlignedSTLAllocator<uint16_t,64>>> vdata(m, std::vector<uint16_t,AlignedSTLAllocator<uint16_t,64>>(n));
// #else
//     std::vector<std::vector<uint16_t>> vdata(m, std::vector<uint16_t>(n));
// #endif
// #ifdef ALIGN
//     for(auto & x : vdata) {
//       assert(get_alignment(x.data()) == 64);
//     }
// #endif

    const uint32_t best_alignment = STORM_get_alignment();
    STORM_ALIGN(64) uint16_t** vdata = (uint16_t**)STORM_aligned_malloc(best_alignment, m*sizeof(uint16_t*));
    for (int i = 0; i < m; ++i)
        vdata[i] = (uint16_t*)STORM_aligned_malloc(best_alignment, n*sizeof(uint16_t));

    // if(verbose) {
    //   printf("alignments: ");
    //   for(auto & x : vdata) {
    //     printf("%d ", get_alignment(x.data()));
    //   }
    //   printf("\n");
    // }

    if (!tabular) printf("alignments: %d\n", best_alignment);

    evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
    evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
    evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
    evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
    evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
    evts.push_back(PERF_COUNT_HW_REF_CPU_CYCLES);
    LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
    std::vector<unsigned long long> results; // tmp buffer
    std::vector< std::vector<unsigned long long> > allresults;
    results.resize(evts.size());
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 0xFFFF);

    bool isok = true;
    for (uint32_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < m; k++) {
            for(size_t k2 = 0; k2 < n ; k2++) { 
                vdata[k][k2] = dis(gen); // random init.
            }
        }

        std::vector<std::vector<uint32_t>> flags(m, std::vector<uint32_t>(16*2));
        
        unified.start();
        for (size_t k = 0; k < m ; k++) {
            fn(vdata[k], n, flags[k].data());
        }
        unified.end(results);
        allresults.push_back(results);
    }

    std::vector<unsigned long long> mins = compute_mins(allresults);
    std::vector<double> avg = compute_averages(allresults);
    
    if (verbose) {
        printf("instructions per cycle %4.2f, cycles per 16-bit word:  %4.3f, "
               "instructions per 16-bit word %4.3f \n",
                double(mins[1]) / mins[0], double(mins[0]) / (n*m), double(mins[1]) / (n*m));
        // first we display mins
        printf("min: %8llu cycles, %8llu instructions, \t%8llu branch mis., %8llu "
               "cache ref., %8llu cache mis.\n",
                mins[0], mins[1], mins[2], mins[3], mins[4]);
        printf("avg: %8.1f cycles, %8.1f instructions, \t%8.1f branch mis., %8.1f "
               "cache ref., %8.1f cache mis.\n",
                avg[0], avg[1], avg[2], avg[3], avg[4]);
    } else {
        printf("cycles per 16-bit word:  %4.3f; ref cycles per 16-bit word: %4.3f \n", double(mins[0]) / (n*m), double(mins[5]) / (n*m));
    }

    if (tabular) {
        for (int i = 0; i < iterations; ++i) {
            printf("%s\t%d\t", fn_name.c_str(), i);
            printf("%4.2f\t%4.3f\t%4.3f\t",
                    double(allresults[i][1]) / allresults[i][0], double(allresults[i][0]) / (n*m), double(allresults[i][1]) / (n*m));
            printf("%llu\t%llu\t%llu\t%llu\t%llu\n",
                    allresults[i][0], allresults[i][1], allresults[i][2], allresults[i][3], allresults[i][4]);
        }
    }

    for (int i = 0; i < m; ++i)
        STORM_aligned_free(vdata[i]);
    STORM_aligned_free(vdata);

    return isok;
}

void measureoverhead(uint32_t n, uint32_t iterations, bool verbose) {
    std::vector<int> evts;
    evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
    evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
    evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
    evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
    evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
    evts.push_back(PERF_COUNT_HW_REF_CPU_CYCLES);
    LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
    std::vector<unsigned long long> results; // tmp buffer
    std::vector< std::vector<unsigned long long> > allresults;
    results.resize(evts.size());
    
    for (uint32_t i = 0; i < iterations; i++) {
        unified.start();
        unified.end(results);
        allresults.push_back(results);
    }

    std::vector<unsigned long long> mins = compute_mins(allresults);
    std::vector<double> avg = compute_averages(allresults);
    printf("%-40s\t","nothing");    
    
    if (verbose) {
        printf("instructions per cycle %4.2f, cycles per 16-bit word:  %4.3f, "
               "instructions per 16-bit word %4.3f \n",
                double(mins[1]) / mins[0], double(mins[0]) / n, double(mins[1]) / n);
        // first we display mins
        printf("min: %8llu cycles, %8llu instructions, \t%8llu branch mis., %8llu "
               "cache ref., %8llu cache mis.\n",
                mins[0], mins[1], mins[2], mins[3], mins[4]);
        printf("avg: %8.1f cycles, %8.1f instructions, \t%8.1f branch mis., %8.1f "
               "cache ref., %8.1f cache mis.\n",
                avg[0], avg[1], avg[2], avg[3], avg[4]);
    } else {
        printf("cycles per 16-bit word:  %4.3f; ref cycles per 16-bit word: %4.3f \n", double(mins[0]) / n, double(mins[5]) / n);
    }
}

bool benchmarkMemoryCopy(const std::string& fn_name, uint32_t n, uint32_t m, uint32_t iterations, bool verbose, bool tabular) {
    std::vector<int> evts;
    const uint32_t best_alignment = STORM_get_alignment();
    STORM_ALIGN(64) uint16_t** vdata = (uint16_t**)STORM_aligned_malloc(best_alignment, m*sizeof(uint16_t*));
    for (int i = 0; i < m; ++i)
        vdata[i] = (uint16_t*)STORM_aligned_malloc(best_alignment, n*sizeof(uint16_t));
    STORM_ALIGN(64) uint16_t* dst = (uint16_t*)STORM_aligned_malloc(best_alignment, n*sizeof(uint16_t*));
        
    if (!tabular) printf("alignments: %d\n", best_alignment);

    evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
    evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
    evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
    evts.push_back(PERF_COUNT_HW_CACHE_REFERENCES);
    evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
    evts.push_back(PERF_COUNT_HW_REF_CPU_CYCLES);
    LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
    std::vector<unsigned long long> results; // tmp buffer
    std::vector< std::vector<unsigned long long> > allresults;
    std::vector<uint32_t> times;
    results.resize(evts.size());
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 0xFFFF);

    bool isok = true;
    for (uint32_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < m; k++) {
            for(size_t k2 = 0; k2 < n ; k2++) {
                vdata[k][k2] = dis(gen); // random init.
            }
        }

        // std::vector<std::vector<uint32_t>> flags(m, std::vector<uint32_t>(16*2));
        
        const clockdef t1 = std::chrono::high_resolution_clock::now();
        unified.start();
        for (size_t k = 0; k < m ; k++) {
            memcpy(dst, vdata[k], n*sizeof(uint16_t));
        }
        unified.end(results);
        const clockdef t2 = std::chrono::high_resolution_clock::now();
        allresults.push_back(results);

        const auto time_span = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
        times.push_back(time_span.count());
    }
    
    uint32_t tot_time = std::accumulate(times.begin(), times.end(), 0);
    double mean_time = tot_time / times.size();

    std::vector<unsigned long long> mins = compute_mins(allresults);
    std::vector<double> avg = compute_averages(allresults);

    double throughput = ((2*n) / (1024*1024.0)) / (mean_time / 1000000000.0);
    
    if (tabular) {
        for (int i = 0; i < iterations; ++i) {
            throughput = ((2*n) / (1024*1024.0)) / (times[i] / 1000000000.0);
            printf("%s\t%u\t%d\t", fn_name.c_str(), n, i);
            printf("%4.2f\t%4.3f\t%4.3f\t",
                    double(allresults[i][1]) / allresults[i][0], double(allresults[i][0]) / (n*m), double(allresults[i][1]) / (n*m));
            printf("%llu\t%llu\t%llu\t%llu\t%llu\t",
                    allresults[i][0], allresults[i][1], allresults[i][2], allresults[i][3], allresults[i][4]);
            printf("%u\t%4.2f\n", times[i], throughput);
        }
    } else if (verbose) {
        printf("instructions per cycle %4.2f, cycles per 16-bit word:  %4.3f, "
               "instructions per 16-bit word %4.3f \n",
                double(mins[1]) / mins[0], double(mins[0]) / (n*m), double(mins[1]) / (n*m));
        // first we display mins
        printf("min: %8llu cycles, %8llu instructions, \t%8llu branch mis., %8llu "
               "cache ref., %8llu cache mis.\n",
                mins[0], mins[1], mins[2], mins[3], mins[4]);
        printf("avg: %8.1f cycles, %8.1f instructions, \t%8.1f branch mis., %8.1f "
               "cache ref., %8.1f cache mis.\n",
                avg[0], avg[1], avg[2], avg[3], avg[4]);
        printf("avg time: %f ns, %4.2f mb/s\n", mean_time, throughput);
    } else {
        printf("cycles per 16-bit word:  %4.3f; ref cycles per 16-bit word: %4.3f \n", double(mins[0]) / (n*m), double(mins[5]) / (n*m));
    }

    return isok;

    for (int i = 0; i < m; ++i) STORM_aligned_free(vdata[i]);
    STORM_aligned_free(vdata);
    STORM_aligned_free(dst);

    return isok;
}

static void print_usage(char *command) {
    printf(" Try %s -n 100000 -i 15 -v \n", command);
    printf("-n is the number of 16-bit words \n");
    printf("-i is the number of tests or iterations \n");
    printf("-v makes things verbose\n");
}

int main(int argc, char **argv) {
    size_t n = 100000;
    size_t m = 1;
    size_t iterations = 0; 
    bool verbose = false;
    bool tabular = false;
    int c;

    while ((c = getopt(argc, argv, "vthm:n:i:")) != -1) {
        switch (c) {
        case 'n':
            n = atoll(optarg);
            break;
        case 'm':
            m = atoll(optarg);
            break;
        case 'v':
            verbose = true;
            break;
        case 't':
            tabular = true;
            break;
        case 'h':
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        case 'i':
            iterations = atoi(optarg);
            break;
        default:
            abort();
        }
    }

    if(n > UINT32_MAX) {
       printf("setting n to %u \n", UINT32_MAX);
       n = UINT32_MAX;
    }

    if(iterations > UINT32_MAX) {
       printf("setting iterations to %u \n", UINT32_MAX);
       iterations = UINT32_MAX;
    }

    if(iterations == 0) {
      if(m*n < 1000000) iterations = 100;
      else iterations = 10;
    }
    
    if (!tabular) {
        printf("n = %zu m = %zu \n", n, m);
        printf("iterations = %zu \n", iterations);
    }
    if(n == 0) {
       printf("n cannot be zero.\n");
       return EXIT_FAILURE;
    }

    if (!tabular) {
        size_t array_in_bytes = sizeof(uint16_t) * n * m;
        if(array_in_bytes < 1024) {
        printf("array size: %zu B\n", array_in_bytes);
        } else if (array_in_bytes < 1024 * 1024) {
        printf("array size: %.3f kB\n", array_in_bytes / 1024.);
        } else {
        printf("array size: %.3f MB\n", array_in_bytes / (1024 * 1024.));
        }
    }

    if (!tabular) measureoverhead(n*m, iterations, verbose);

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

    if (!tabular) printf("libflagstats-scalar\t");
    fflush(NULL);
    bool isok = benchmarkMany("libflagstats-scalar", n, m, iterations, FLAGSTAT_scalar, verbose, true, tabular);
    if (isok == false) {
        printf("Problem detected with %u.\n", 0);
    }
    if (verbose && !tabular)
        printf("\n");

    if (!tabular) printf("samtools\t");
    fflush(NULL);
    isok = benchmarkMany("samtools", n, m, iterations, samtools_flagstats, verbose, true, tabular);
    if (isok == false) {
        printf("Problem detected with %u.\n", 0);
    }
    if (verbose && !tabular)
        printf("\n");

    if (!tabular) printf("memcpy\t");
    fflush(NULL);
    isok = benchmarkMemoryCopy("memcpy", n, m, iterations, verbose, tabular);
    if (isok == false) {
        printf("Problem detected with %u.\n", 0);
    }
    if (verbose && !tabular)
        printf("\n");
    

    #if defined(STORM_HAVE_SSE42)
        if ((cpuid & STORM_CPUID_runtime_bit_SSE42)) {
            if (!tabular) printf("libflagstats-sse4.2\t");
            fflush(NULL);
            bool isok = benchmarkMany("libflagstats-sse4.2", n, m, iterations, FLAGSTAT_sse4, verbose, true, tabular);
            if (isok == false) {
                printf("Problem detected with %u.\n", 0);
            }
            if (verbose && !tabular) printf("\n");

            if (!tabular) printf("libflagstats-sse4.2-optimized\t");
            fflush(NULL);
            isok = benchmarkMany("libflagstats-sse4.2-optimized", n, m, iterations, FLAGSTAT_sse4_improved, verbose, true, tabular);
            if (isok == false) {
                printf("Problem detected with %u.\n", 0);
            }
            if (verbose && !tabular) printf("\n");
        }
    // }
    #endif

    #if defined(STORM_HAVE_AVX2)
        if ((cpuid & STORM_CPUID_runtime_bit_AVX2)) {
            if (!tabular) printf("libflagstats-avx2\t");
            fflush(NULL);
            bool isok = benchmarkMany("libflagstats-avx2", n, m, iterations, FLAGSTAT_avx2, verbose, true, tabular);
            if (isok == false) {
                printf("Problem detected with %u.\n", 1);
            }
            if (verbose && !tabular)
                printf("\n");
        }
    // }
    #endif

    #if defined(STORM_HAVE_AVX512)
        if ((cpuid & STORM_CPUID_runtime_bit_AVX512BW)) {
            if (!tabular) printf("libflagstats-avx512bw\t");
            fflush(NULL);
            bool isok = benchmarkMany("libflagstats-avx512bw", n, m, iterations, FLAGSTAT_avx512, verbose, true, tabular);
            if (isok == false) {
                printf("Problem detected with %u.\n", 2);
            }
            if (verbose && !tabular) printf("\n");

            if (!tabular) printf("libflagstats-avx512bw-improved\t");
            fflush(NULL);
            isok = benchmarkMany("libflagstats-avx512bw-improved", n, m, iterations, FLAGSTAT_avx512_improved, verbose, true, tabular);
            if (isok == false) {
                printf("Problem detected with %u.\n", 2);
            }
            if (verbose && !tabular) printf("\n");

            if (!tabular) printf("libflagstats-avx512bw-improved2\t");
            fflush(NULL);
            isok = benchmarkMany("libflagstats-avx512bw-improved2", n, m, iterations, FLAGSTAT_avx512_improved2, verbose, true, tabular);
            if (isok == false) {
                printf("Problem detected with %u.\n", 2);
            }
            if (verbose && !tabular) printf("\n");
        // }
        }
    #endif

    if (!verbose && !tabular)
        printf("Try -v to get more details.\n");

    return EXIT_SUCCESS;
}
#else //  __linux__

#include <stdio.h>
#include <stdlib.h>

int main() {
    printf("This is a linux-specific benchmark\n");
    return EXIT_SUCCESS;
}

#endif
