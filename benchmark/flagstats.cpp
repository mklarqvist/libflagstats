/*
* Copyright (c) 2019
* Author(s): Marcus D. R. Klarqvist
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*   http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing,
* software distributed under the License is distributed on an
* "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, either express or implied.  See the License for the
* specific language governing permissions and limitations
* under the License.
*/

#include <cstring> // memcpy
#include <stdio.h>  // For printf()
#include <string.h> // For memcmp()
#include <stdlib.h> // For exit()
#include <getopt.h> // options

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <chrono>//timer

#include <unistd.h>//sync
#include <bitset>

#include "lz4.h" // lz4
#include "lz4hc.h"
#include "zstd.h" // zstd
// #include "zstd_errors.h"
#include "libalgebra.h" // pospopcnt
#include "libflagstats.h" // flagstats

// samtools count flagstat different:
// https://github.com/samtools/samtools/blob/master/bam_stat.c#L47
typedef struct {
    long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
    long long n_sgltn[2], n_read1[2], n_read2[2];
    long long n_dup[2];
    long long n_diffchr[2], n_diffhigh[2];
    long long n_secondary[2], n_supp[2];
} bam_flagstat_t;

#define flagstat_loop(s, c) do {                                \
        int w = (c & FLAGSTAT_FQCFAIL)? 1 : 0;                       \
        ++(s)->n_reads[w];                                      \
        if (c & FLAGSTAT_FSECONDARY ) {                              \
            ++(s)->n_secondary[w];                              \
        } else if (c & FLAGSTAT_FSUPPLEMENTARY ) {                   \
            ++(s)->n_supp[w];                                   \
        } else if (c & FLAGSTAT_FPAIRED) {                           \
            ++(s)->n_pair_all[w];                               \
            if ((c & FLAGSTAT_FPROPER_PAIR) && !(c & FLAGSTAT_FUNMAP) ) ++(s)->n_pair_good[w]; \
            if (c & FLAGSTAT_FREAD1) ++(s)->n_read1[w];              \
            if (c & FLAGSTAT_FREAD2) ++(s)->n_read2[w];              \
            if ((c & FLAGSTAT_FMUNMAP) && !(c & FLAGSTAT_FUNMAP)) ++(s)->n_sgltn[w]; \
            if (!(c & FLAGSTAT_FUNMAP) && !(c & FLAGSTAT_FMUNMAP)) {      \
                ++(s)->n_pair_map[w];                           \
            }                                                   \
        }                                                       \
        if (!(c & FLAGSTAT_FUNMAP)) ++(s)->n_mapped[w];              \
        if (c & FLAGSTAT_FDUP) ++(s)->n_dup[w];                      \
} while (0)


static const char* percent(char* buffer, long long n, long long total)
{
    if (total != 0) sprintf(buffer, "%.2f%%", (float)n / total * 100.0);
    else strcpy(buffer, "N/A");
    return buffer;
}


// @see: https://stackoverflow.com/questions/6818606/how-to-programmatically-clear-the-filesystem-memory-cache-in-c-on-a-linux-syst
void clear_cache() {
#ifdef __linux__ 
    sync();
    std::ofstream ofs("/proc/sys/vm/drop_caches");
    ofs << "3" << std::endl;
#endif
}
 
int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int ret = ZSTD_decompress(out, out_capacity, in, n_in);
    return(ret);
}

static const std::string SAM_FLAG_NAME[] = {"FPAIRED","FPROPER_PAIR","FUNMAP","FMUNMAP","FREVERSE","FMREVERSE", "FREAD1","FREAD2","FSECONDARY","FQCFAIL","FDUP","FSUPPLEMENTARY","n_pair_good","n_sgltn","n_pair_map"};

/*
 * Easy show-error-and-bail function.
 */
void run_screaming(const char* message, const int code) {
    printf("%s \n", message);
    exit(code);
}

int lz4f(const std::string& file, const std::string& out_prefix, const int acceleration = 2) {
    std::ifstream f(file, std::ios::in | std::ios::binary);
    if (f.good() == false) return 0;

    std::string outfile = out_prefix + "_fast_a" + std::to_string(acceleration) + ".lz4";
    std::cerr << "Opening=" << outfile << std::endl;
    std::ofstream of(outfile, std::ios::out | std::ios::binary);
    if (of.good() == false) return 0;

    uint8_t buffer[1024000]; // 512k 16-bit ints 
    const int max_dst_size = LZ4_compressBound(1024000);
    uint8_t* out_buffer = new uint8_t[max_dst_size];

    while (f.good()) {
        f.read((char*)buffer, 1024000);
        int32_t bytes_read = f.gcount();

        const int32_t compressed_data_size = LZ4_compress_fast((char*)buffer, (char*)out_buffer, bytes_read, max_dst_size, acceleration);
        // Check return_value to determine what happened.
        
        if (compressed_data_size < 0)
            run_screaming("A negative result from LZ4_compress_default indicates a failure trying to compress the data.  See exit code (echo $?) for value returned.", compressed_data_size);
        
        if (compressed_data_size == 0)
            run_screaming("A result of 0 means compression worked, but was stopped because the destination buffer couldn't hold all the information.", 1);

        of.write((char*)&bytes_read, sizeof(int32_t));
        of.write((char*)&compressed_data_size, sizeof(int32_t));
        of.write((char*)out_buffer, compressed_data_size);
        
        // std::cerr << "Compressed " << bytes_read << "->" << compressed_data_size << std::endl;
    }

    delete[] out_buffer;
    return 1;
}

int lz4hc(const std::string& file, const std::string& out_prefix, int clevel = 9) {
    std::ifstream f(file, std::ios::in | std::ios::binary);
    if (f.good() == false) return 0;

    std::string outfile = out_prefix + "_HC_c" + std::to_string(clevel) + ".lz4";
    std::cerr << "Opening=" << outfile << std::endl;
    std::ofstream of(outfile, std::ios::out | std::ios::binary);
    if (of.good() == false) {
        std::cerr << "outfile not good" << std::endl;
        return 0;
    }

    std::cerr << "here1" << std::endl;

    uint8_t buffer[1024000]; // 512k 16-bit ints 
    const int max_dst_size = LZ4_compressBound(1024000);
    uint8_t* out_buffer = new uint8_t[max_dst_size];

    std::cerr << "here" << std::endl;

    while (f.good()) {
        f.read((char*)buffer, 1024000);
        int32_t bytes_read = f.gcount();
        std::cerr << "read=" << bytes_read << std::endl;

        const int32_t compressed_data_size = LZ4_compress_HC((char*)buffer, (char*)out_buffer, bytes_read, max_dst_size, clevel);
        
        // Check return_value to determine what happened.
        if (compressed_data_size < 0)
            run_screaming("A negative result from LZ4_compress_default indicates a failure trying to compress the data. See exit code (echo $?) for value returned.", compressed_data_size);
        
        if (compressed_data_size == 0)
            run_screaming("A result of 0 means compression worked, but was stopped because the destination buffer couldn't hold all the information.", 1);

        of.write((char*)&bytes_read, sizeof(int32_t));
        of.write((char*)&compressed_data_size, sizeof(int32_t));
        of.write((char*)out_buffer, compressed_data_size);
        
        // std::cerr << "Compressed " << bytes_read << "->" << compressed_data_size << std::endl;
    }

    delete[] out_buffer;
    return 1;
}

int zstd(const std::string& file, const std::string& out_prefix, int clevel = 22) {
    std::ifstream f(file, std::ios::in | std::ios::binary);
    if (f.good() == false) return 0;

    std::string outfile = out_prefix + "_c" + std::to_string(clevel) + ".zst";
    std::cerr << "Opening=" << outfile << std::endl;
    std::ofstream of(outfile, std::ios::out | std::ios::binary);
    if (of.good() == false) return 0;

    uint8_t buffer[1024000]; // 512k 16-bit ints 
    uint8_t out_buffer[1024000];

    while (f.good()) {
        f.read((char*)buffer, 1024000);
        int32_t bytes_read = f.gcount();

        const int32_t compressed_data_size = ZstdCompress(buffer, bytes_read, out_buffer, 1024000, clevel);
        // Check return_value to determine what happened.
        
        if (compressed_data_size < 0)
            run_screaming("A negative result from LZ4_compress_default indicates a failure trying to compress the data.  See exit code (echo $?) for value returned.", compressed_data_size);
        
        if (compressed_data_size == 0)
            run_screaming("A result of 0 means compression worked, but was stopped because the destination buffer couldn't hold all the information.", 1);
        
        of.write((char*)&bytes_read, sizeof(int32_t));
        of.write((char*)&compressed_data_size, sizeof(int32_t));
        of.write((char*)out_buffer, compressed_data_size);

        // std::cerr << "Compressed " << bytes_read << "->" << compressed_data_size << std::endl;
    }
    of.close();

    return 1;
}

int lz4_decompress_only(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    // uint8_t buffer[1024000]; // 512k 16-bit ints 
    // uint8_t out_buffer[1024000];
    uint8_t* buffer = new uint8_t[1024000+65536];
    // uint8_t* out_buffer;
    // assert(!posix_memalign((void**)&out_buffer, SIMD_ALIGNMENT, 1024000));
    // out_buffer = new uint8_t[1024000];
    uint8_t* out_buffer = (uint8_t*)STORM_aligned_malloc(STORM_get_alignment(), 1024000+65536);

    int32_t uncompresed_size, compressed_size;

    uint32_t counters[16] = {0}; // flags
    uint64_t tot_flags = 0;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        // std::cerr << "here" << std::endl;
        f.read((char*)&uncompresed_size, sizeof(int32_t));
        f.read((char*)&compressed_size, sizeof(int32_t));
        f.read((char*)buffer, compressed_size);

        const int32_t decompressed_size = LZ4_decompress_safe((char*)buffer, (char*)out_buffer, compressed_size, uncompresed_size);
        // Check return_value to determine what happened.
        if (decompressed_size < 0)
            run_screaming("A negative result from LZ4_decompress_safe indicates a failure trying to decompress the data.  See exit code (echo $?) for value returned.", decompressed_size);
        if (decompressed_size == 0)
            run_screaming("I'm not sure this function can ever return 0.  Documentation in lz4.h doesn't indicate so.", 1);

        // assert(decompressed_size == uncompresed_size);

        const uint32_t N = uncompresed_size >> 1;
        tot_flags += N;
        // pospopcnt_u16((uint16_t*)out_buffer,N,counters);

        // std::cerr << "Decompressed " << compressed_size << "->" << uncompresed_size << std::endl;
        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[LZ4 " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    // std::cerr << "Tot flags=" << tot_flags << std::endl;
    // for (int i = 0; i < 12; ++i) {
    //     std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << std::endl;
    // }

    delete[] buffer;
    // delete[] out_buffer;
    // free(buffer);
    // free(out_buffer);
    STORM_aligned_free(out_buffer);
    return 1;
}

int lz4_decompress(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    // uint8_t buffer[1024000]; // 512k 16-bit ints 
    // uint8_t out_buffer[1024000];
    
    uint8_t* buffer = new uint8_t[1024000+65536];
    // uint8_t* out_buffer;
    // assert(!posix_memalign((void**)&out_buffer, SIMD_ALIGNMENT, 1024000));
    // out_buffer = new uint8_t[1024000];
    uint8_t* out_buffer = (uint8_t*)STORM_aligned_malloc(STORM_get_alignment(), 1024000+65536);

    int32_t uncompresed_size, compressed_size;

    uint32_t counters[16*2] = {0}; // flags
    uint64_t tot_flags = 0;

    FLAGSTATS_func func;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        while (f.good()) {
        f.read((char*)&uncompresed_size, sizeof(int32_t));
        f.read((char*)&compressed_size, sizeof(int32_t));
        f.read((char*)buffer, compressed_size);

        const int32_t decompressed_size = LZ4_decompress_safe((char*)buffer, (char*)out_buffer, compressed_size, uncompresed_size);
        // Check return_value to determine what happened.
        if (decompressed_size < 0)
            run_screaming("A negative result from LZ4_decompress_safe indicates a failure trying to decompress the data.  See exit code (echo $?) for value returned.", decompressed_size);
        if (decompressed_size == 0)
            run_screaming("I'm not sure this function can ever return 0.  Documentation in lz4.h doesn't indicate so.", 1);

        const uint32_t N = uncompresed_size >> 1;
        tot_flags += N;
        // pospopcnt_u16((uint16_t*)out_buffer,N,counters);
        // FLAGSTAT_avx512((uint16_t*)out_buffer,N,counters);
        // STORM_pospopcnt_u16_avx2_harvey_seal((uint16_t*)out_buffer,N,counters);
        func = FLAGSTATS_get_function(N);
        (*func)((uint16_t*)out_buffer,N,counters);
        // std::cerr << "Decompressed " << compressed_size << "->" << uncompresed_size << std::endl;
        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[LZ4 " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    std::cerr << "Tot flags=" << tot_flags << std::endl;
    // std::cerr << "Pass QC" << std::endl;
    for (int i = 0; i < 15; ++i) {
        std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << "\t" << counters[16+i] << std::endl;
    }
    // for (int i = 12; i < 15; ++i) {
    //     std::cerr << "special-" << i << "\t" << counters[i] << "\t" << counters[16+i] << std::endl;
    // }

    // std::cerr << "Fail QC" << std::endl;
    // for (int i = 0; i < 12; ++i) {
    //     std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[16+i] << std::endl;
    // }

    delete[] buffer;
    // delete[] out_buffer;
    // free(buffer);
    // free(out_buffer);
    STORM_aligned_free(out_buffer);
    return 1;
}

int flagstat_raw_read(const std::string& file) {
    std::size_t found = file.find(".bin");
    std::string file2;
    if (found != std::string::npos) {
        // std::cerr << "first 'needle' found at: " << found << '\n';
        std::cerr << "file new=" << file.substr(0, found + 4) << std::endl;
        file2 = file.substr(0, found + 4);
    } else {
        return -1;
    }
    
    std::ifstream f(file2, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) {
        std::cerr << "file not good" << std::endl;
        return 0;
    }
    int64_t filesize = f.tellg();
    f.seekg(0);
    // std::cerr << "filesize=" << filesize << std::endl;
    
    // uint8_t* buffer = new uint8_t[1024000];
    uint8_t* out_buffer = (uint8_t*)STORM_aligned_malloc(STORM_get_alignment(), 1024000);
    uint32_t counters[16*2] = {0}; // flags
    uint64_t tot_flags = 0;

    // FLAGSTATS_func func;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)out_buffer, 1024000);
        size_t read = f.gcount();
        // std::cerr << "Read " << read << std::endl;

        // FLAGSTAT_avx512((uint16_t*)out_buffer, read/2, counters);
        // func = FLAGSTATS_get_function(read >> 1);
        // (*func)((uint16_t*)out_buffer,read >> 1,counters);

        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[RAW READ " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    // std::cerr << "Tot flags=" << tot_flags << std::endl;
    // std::cerr << "Pass QC" << std::endl;
    // for (int i = 0; i < 15; ++i) {
    //     std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << "\t" << counters[16+i] << std::endl;
    // }

    STORM_aligned_free(out_buffer);
    return 1;
}

int flagstat_raw(const std::string& file) {
    std::size_t found = file.find(".bin");
    std::string file2;
    if (found != std::string::npos) {
        // std::cerr << "first 'needle' found at: " << found << '\n';
        std::cerr << "file new=" << file.substr(0, found + 4) << std::endl;
        file2 = file.substr(0, found + 4);
    } else {
        return -1;
    }
    
    std::ifstream f(file2, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) {
        std::cerr << "file not good" << std::endl;
        return 0;
    }
    int64_t filesize = f.tellg();
    f.seekg(0);
    // std::cerr << "filesize=" << filesize << std::endl;
    
    // uint8_t* buffer = new uint8_t[1024000];
    uint8_t* out_buffer = (uint8_t*)STORM_aligned_malloc(STORM_get_alignment(), 1024000);
    uint32_t counters[16*2] = {0}; // flags
    uint64_t tot_flags = 0;

    FLAGSTATS_func func;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)out_buffer, 1024000);
        size_t read = f.gcount();
        // std::cerr << "Read " << read << std::endl;

        // FLAGSTAT_avx512((uint16_t*)out_buffer, read/2, counters);
        func = FLAGSTATS_get_function(read >> 1);
        (*func)((uint16_t*)out_buffer,read >> 1,counters);

        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[RAW " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    std::cerr << "Tot flags=" << tot_flags << std::endl;
    // std::cerr << "Pass QC" << std::endl;
    for (int i = 0; i < 15; ++i) {
        std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << "\t" << counters[16+i] << std::endl;
    }

    STORM_aligned_free(out_buffer);
    return 1;
}

int flagstat_raw_samtools(const std::string& file) {
    std::size_t found = file.find(".bin");
    std::string file2;
    if (found != std::string::npos) {
        // std::cerr << "first 'needle' found at: " << found << '\n';
        std::cerr << "file new=" << file.substr(0, found + 4) << std::endl;
        file2 = file.substr(0, found + 4);
    } else {
        return -1;
    }
    
    std::ifstream f(file2, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) {
        std::cerr << "file not good" << std::endl;
        return 0;
    }
    int64_t filesize = f.tellg();
    f.seekg(0);
    // std::cerr << "filesize=" << filesize << std::endl;
    
    // uint8_t* buffer = new uint8_t[1024000];
    uint8_t* out_buffer = (uint8_t*)STORM_aligned_malloc(STORM_get_alignment(), 1024000);
    uint32_t counters[16*2] = {0}; // flags
    uint64_t tot_flags = 0;

    // FLAGSTATS_func func;
    bam_flagstat_t* s;
    s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)out_buffer, 1024000);
        size_t read = f.gcount();
        // std::cerr << "Read " << read << std::endl;

        // FLAGSTAT_avx512((uint16_t*)out_buffer, read/2, counters);
        // func = FLAGSTATS_get_function(read >> 1);
        // (*func)((uint16_t*)out_buffer,read >> 1,counters);
        uint16_t* inflags = (uint16_t*)out_buffer;
        for (int i = 0; i < ((uint32_t)read >> 1); ++i) {
            flagstat_loop(s, inflags[i]);
        }

        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[RAW SAMTOOLS " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    std::cerr << "Tot flags=" << tot_flags << std::endl;
    // std::cerr << "Pass QC" << std::endl;
    for (int i = 0; i < 15; ++i) {
        std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << "\t" << counters[16+i] << std::endl;
    }

    STORM_aligned_free(out_buffer);
    free(s);
    return 1;
}

int lz4_decompress_samtools(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    uint8_t buffer[1024000+65536]; // 512k 16-bit ints 
    uint8_t out_buffer[1024000+65536];

    int32_t uncompresed_size, compressed_size;

    // uint32_t counters[16] = {0}; // flags
    uint64_t tot_flags = 0;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    
    bam_flagstat_t* s;
    s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));

    while (f.good()) {
        f.read((char*)&uncompresed_size, sizeof(int32_t));
        f.read((char*)&compressed_size, sizeof(int32_t));
        f.read((char*)buffer, compressed_size);

        const int32_t decompressed_size = LZ4_decompress_safe((char*)buffer, (char*)out_buffer, compressed_size, uncompresed_size);
        // Check return_value to determine what happened.
        if (decompressed_size < 0)
            run_screaming("A negative result from LZ4_decompress_safe indicates a failure trying to decompress the data. See exit code (echo $?) for value returned.", decompressed_size);
        if (decompressed_size == 0)
            run_screaming("I'm not sure this function can ever return 0. Documentation in lz4.h doesn't indicate so.", 1);

        const uint32_t N = uncompresed_size >> 1;
        tot_flags += N;

        uint16_t* inflags = (uint16_t*)out_buffer;
        for (int i = 0; i < N; ++i) {
            flagstat_loop(s, inflags[i]);
        }

        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[LZ4 samtools " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    char b0[16], b1[16];
    printf("%lld + %lld in total (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
    printf("%lld + %lld secondary\n", s->n_secondary[0], s->n_secondary[1]);
    printf("%lld + %lld supplementary\n", s->n_supp[0], s->n_supp[1]);
    printf("%lld + %lld duplicates\n", s->n_dup[0], s->n_dup[1]);
    printf("%lld + %lld mapped (%s : %s)\n", s->n_mapped[0], s->n_mapped[1], percent(b0, s->n_mapped[0], s->n_reads[0]), percent(b1, s->n_mapped[1], s->n_reads[1]));
    printf("%lld + %lld paired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
    printf("%lld + %lld read1\n", s->n_read1[0], s->n_read1[1]);
    printf("%lld + %lld read2\n", s->n_read2[0], s->n_read2[1]);
    printf("%lld + %lld properly paired (%s : %s)\n", s->n_pair_good[0], s->n_pair_good[1], percent(b0, s->n_pair_good[0], s->n_pair_all[0]), percent(b1, s->n_pair_good[1], s->n_pair_all[1]));
    printf("%lld + %lld with itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
    printf("%lld + %lld singletons (%s : %s)\n", s->n_sgltn[0], s->n_sgltn[1], percent(b0, s->n_sgltn[0], s->n_pair_all[0]), percent(b1, s->n_sgltn[1], s->n_pair_all[1]));
    // printf("%lld + %lld with mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
    // printf("%lld + %lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
    free(s);

    return 1;
}

int zstd_decompress_only(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    uint8_t buffer[1024000]; // 512k 16-bit ints 
    uint8_t out_buffer[1024000];

    int32_t uncompresed_size, compressed_size;
    uint64_t tot_flags = 0; 
    
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)&uncompresed_size, sizeof(int32_t));
        f.read((char*)&compressed_size, sizeof(int32_t));
        f.read((char*)buffer, compressed_size);

        const int32_t decompressed_size = ZstdDecompress(buffer, 1024000, out_buffer, uncompresed_size);
        // assert(decompressed_size == uncompresed_size);

        const uint32_t N = uncompresed_size >> 1;
        tot_flags += N;

        // std::cerr << "Decompressed " << compressed_size << "->" << uncompresed_size << std::endl;
        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[ZSTD " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    // std::cerr << "Tot flags=" << tot_flags << std::endl;
    // for (int i = 0; i < 12; ++i) {
    //     std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << std::endl;
    // }

    return 1;
}

int zstd_decompress(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    uint8_t buffer[1024000]; // 512k 16-bit ints 
    uint8_t out_buffer[1024000];

    int32_t uncompresed_size, compressed_size;

    uint32_t counters[16*2] = {0}; // flags
    uint64_t tot_flags = 0;


    FLAGSTATS_func func;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)&uncompresed_size, sizeof(int32_t));
        f.read((char*)&compressed_size, sizeof(int32_t));
        f.read((char*)buffer, compressed_size);

        const int32_t decompressed_size = ZstdDecompress(buffer, 1024000, out_buffer, uncompresed_size);
        // assert(decompressed_size == uncompresed_size);

        const uint32_t N = uncompresed_size >> 1;
        tot_flags += N;
        // pospopcnt_u16((uint16_t*)out_buffer,N,counters);
        func = FLAGSTATS_get_function(N);
        (*func)((uint16_t*)out_buffer,N,counters);

        // std::cerr << "Decompressed " << compressed_size << "->" << uncompresed_size << std::endl;
        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[ZSTD " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    // std::cerr << "Tot flags=" << tot_flags << std::endl;
    // for (int i = 0; i < 12; ++i) {
    //     std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << std::endl;
    // }

    return 1;
}

int zstd_decompress_samtools(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    uint8_t buffer[1024000]; // 512k 16-bit ints 
    uint8_t out_buffer[1024000];

    int32_t uncompresed_size, compressed_size;

    uint32_t counters[16] = {0}; // flags
    uint64_t tot_flags = 0;

    bam_flagstat_t *s;
    s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t)); 

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)&uncompresed_size, sizeof(int32_t));
        f.read((char*)&compressed_size, sizeof(int32_t));
        f.read((char*)buffer, compressed_size);

        const int32_t decompressed_size = ZstdDecompress(buffer, 1024000, out_buffer, uncompresed_size);
        // assert(decompressed_size == uncompresed_size);

        const uint32_t N = uncompresed_size >> 1;
        tot_flags += N;
        // pospopcnt_u16((uint16_t*)out_buffer,N,counters);
        uint16_t* inflags = (uint16_t*)out_buffer;
        for (int i = 0; i < N; ++i) {
            // flagstat_loop(s, c);
            flagstat_loop(s, inflags[i]);
            // flagstat_loop_branchless(s, inflags, i);
        }

        // std::cerr << "Decompressed " << compressed_size << "->" << uncompresed_size << std::endl;
        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[ZSTD samtools " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    // std::cerr << "Tot flags=" << tot_flags << std::endl;
    // for (int i = 0; i < 12; ++i) {
    //     std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << std::endl;
    // }

    free(s);

    return 1;
}

int compress(int argc, char** argv) {
    int c;
	if(argc < 3){
		std::cerr << "usage" << std::endl;
        return(EXIT_FAILURE);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",       required_argument, 0,  'i' },
		{"output",      optional_argument, 0,  'o' },
		{"compression-level", required_argument, 0,  'c' },
        {"lz4", optional_argument, 0,  'l' },
        {"zstd", optional_argument, 0,  'z' },
        {"fast", optional_argument, 0,  'f' },
		{0,0,0,0}
	};
	
    std::string input, output;
    int clevel = 1;
    bool lz4_fast = false;
    bool lz4 = false;
    bool mzstd = false;

	while ((c = getopt_long(argc, argv, "i:o:c:lzf?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			input = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'c':
			clevel = atoi(optarg);
			if(clevel < 0){
				std::cerr << "illegal clevel=" << clevel << std::endl;
				return(EXIT_FAILURE);
			}
			break;
		case 'f':
			lz4_fast = true;
			break;
		case 'z':
			mzstd = true;
			break;
		case 'l':
			lz4 = true;
			break;

		default:
			std::cerr << "Unrecognized option: " << (char)c << std::endl;
			return(EXIT_FAILURE);
		}
	}

    if (lz4 == false && mzstd == false) {
        std::cerr << "must pick a compression algorithm (lz4 or zstd)" << std::endl;
        return EXIT_FAILURE;
    }

    if (mzstd && lz4_fast) {
        std::cerr << "fast mode is only used with LZ4. Ignoring..." << std::endl;
    }

    if (input.size() == 0) {
        std::cerr << "No input file given" << std::endl;
        return EXIT_FAILURE; 
    }

    if (mzstd) {
        if (output.size() == 0) {
            output = input;
        }
        zstd(input,output, clevel);
    }
    
    if (lz4) {
        if (output.size() == 0) {
            output = input;
        }
        if (lz4_fast) lz4f(input, output, clevel);
        else lz4hc(input, output, clevel);
    }

    return EXIT_SUCCESS;
}

int check_file_extension(const std::string& filename) {
    size_t pos = filename.rfind('.');
    if (pos == std::string::npos)
        return 0;

    std::string ext = filename.substr(pos + 1);

    if (ext == "zst") return 1;
    if (ext == "lz4") return 2;

    return 0;
}

int decompress(int argc, char** argv) {
    int c;
	if(argc < 3){
		std::cerr << "usage decompress" << std::endl;
        return(EXIT_FAILURE);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",       required_argument, 0,  'i' },

        {"raw-read",       optional_argument, 0,  'R' },
        {"raw-flagstats",       optional_argument, 0,  'D' },
        {"raw-samtools",       optional_argument, 0,  'S' },

        {"read",       optional_argument, 0,  'r' },
        {"flagstats",       optional_argument, 0,  'd' },
        {"samtools",       optional_argument, 0,  's' },

		{0,0,0,0}
	};
	
    std::string input;
    bool lz4 = false;
    bool mzstd = false;
    int step = 0;

	while ((c = getopt_long(argc, argv, "i:RDSrds?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			input = std::string(optarg);
			break;
        case 'R': step = 0; break;
        case 'D': step = 1; break;
        case 'S': step = 2; break;
        case 'r': step = 3; break;
        case 'd': step = 4; break;
        case 's': step = 5; break;
		default:
			std::cerr << "Unrecognized option: " << (char)c << std::endl;
			return(EXIT_FAILURE);
		}
	}

    if (input.size() == 0) {
        std::cerr << "No input file given" << std::endl;
        return EXIT_FAILURE; 
    }

    // determine suffix
    int method = check_file_extension(input);
    if (method == 0) {
        std::cerr << "unknown file extension" << std::endl;
        return EXIT_FAILURE;
    }

    if (method == 1) {
        // std::cerr << "method1" << std::endl;
        switch(step) {
        case 0: flagstat_raw_read(input); break; // -R
        case 1: flagstat_raw(input); break;// -D
        case 2: flagstat_raw_samtools(input); break;// -S
        case 3: zstd_decompress_only(input); break;// warmup, -r
        case 4: zstd_decompress(input); break;// -d
        case 5: zstd_decompress_samtools(input); break;// -s
        }
    }
    
    if (method == 2) {
        switch(step) {
        case 0: lz4_decompress_only(input); break; // warmup
        case 1: flagstat_raw_read(input); break;
        case 2: flagstat_raw(input); break;
        case 3: flagstat_raw_samtools(input); break;
        case 4: lz4_decompress(input); break;
        case 5: lz4_decompress_samtools(input); break;
        }
    }

    return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
    if (argc == 1) {
        std::cerr << "usage" << std::endl;
        return EXIT_SUCCESS;
    }

    if(strcmp(&argv[1][0], "compress") == 0){
		// return(compress(argc, argv));
        std::cerr << "compress" << std::endl;
        return(compress(argc,argv));
	} else if(strcmp(&argv[1][0], "decompress") == 0){
		// return(compress(argc, argv));
        std::cerr << "decompress" << std::endl;
        return(decompress(argc,argv));
	}
    else {
        std::cerr << "unknown" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}