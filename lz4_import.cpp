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
#include "lz4.h" // lz4
#include "lz4hc.h"
#include "zstd.h" // zstd
#include "zstd_errors.h"
#include <cstring> // memcpy
#include "pospopcnt.h" // pospopcnt

#include <stdio.h>  // For printf()
#include <string.h> // For memcmp()
#include <stdlib.h> // For exit()

//
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <chrono>//timer

#include <unistd.h>//sync
#include <bitset>

#include "getopt.h" // options

/****************************
*  SIMD definitions
****************************/
#if defined(_MSC_VER)
     /* Microsoft C/C++-compatible compiler */
     #include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
     /* GCC-compatible compiler, targeting x86/x86-64 */
     #include <x86intrin.h>
#elif defined(__GNUC__) && defined(__ARM_NEON__)
     /* GCC-compatible compiler, targeting ARM with NEON */
     #include <arm_neon.h>
#elif defined(__GNUC__) && defined(__IWMMXT__)
     /* GCC-compatible compiler, targeting ARM with WMMX */
     #include <mmintrin.h>
#elif (defined(__GNUC__) || defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
     /* XLC or GCC-compatible compiler, targeting PowerPC with VMX/VSX */
     #include <altivec.h>
#elif defined(__GNUC__) && defined(__SPE__)
     /* GCC-compatible compiler, targeting PowerPC with SPE */
     #include <spe.h>
#endif

#if defined(__AVX512F__) && __AVX512F__ == 1
#define SIMD_AVAILABLE  1
#define SIMD_VERSION    6
#define SIMD_WIDTH      512
#define SIMD_ALIGNMENT  64
#elif defined(__AVX2__) && __AVX2__ == 1
#define SIMD_AVAILABLE  1
#define SIMD_VERSION    5
#define SIMD_WIDTH      256
#define SIMD_ALIGNMENT  32
#elif defined(__AVX__) && __AVX__ == 1
#define SIMD_AVAILABLE  1
#define SIMD_VERSION    4
#define SIMD_ALIGNMENT  16
#define SIMD_WIDTH      128
#elif defined(__SSE4_1__) && __SSE4_1__ == 1
#define SIMD_AVAILABLE  1
#define SIMD_VERSION    3
#define SIMD_ALIGNMENT  16
#define SIMD_WIDTH      128
#elif defined(__SSE2__) && __SSE2__ == 1
#define SIMD_AVAILABLE  0 // unsupported version
#define SIMD_VERSION    0
#define SIMD_ALIGNMENT  16
#define SIMD_WIDTH      128
#elif defined(__SSE__) && __SSE__ == 1
#define SIMD_AVAILABLE  0 // unsupported version
#define SIMD_VERSION    0
#define SIMD_ALIGNMENT  16
#define SIMD_WIDTH      0
#else
#define SIMD_AVAILABLE  0
#define SIMD_VERSION    0
#define SIMD_ALIGNMENT  16
#define SIMD_WIDTH      0
#endif

// From sam.h
/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048

// Support lookup for SIMD masks.
static const uint16_t lookup_sec[2]  = {65535, 256};
static const uint16_t lookup_sup[2]  = {65535, 2048};
static const uint16_t lookup_pair[2] = {256+2048, 65535};

void samtools_single_update(uint16_t val, uint32_t* flags) {
    const int offset = ( (val & BAM_FQCFAIL) == 0 ) ? 0 : 16;

    // Always.
    flags[offset+3] += (val & BAM_FUNMAP) == 0; // this is implicit
    flags[offset+10] += (val & BAM_FDUP) >> 10;
    // Rule 1.
    flags[offset+8] += (val & BAM_FSECONDARY) >> 8;
    val &= lookup_sec[(val & BAM_FSECONDARY) >> 8];
    // Rule 2.
    flags[offset+11] += (val & BAM_FSUPPLEMENTARY) >> 11;
    val &= lookup_sup[(val & BAM_FSUPPLEMENTARY) >> 11];
    val &= lookup_pair[val & BAM_FPAIRED];
    
    // skip 4,10,8,11
    for (int i = 0; i < 4; ++i) {
        flags[offset+i] += ((val & (1 << i)) >> i);
    }

    for (int i = 5; i < 8; ++i) {
        flags[offset+i] += ((val & (1 << i)) >> i);
    }

    for (int i = 9; i < 11; ++i) {
        flags[offset+i] += ((val & (1 << i)) >> i);
    }

    for (int i = 12; i < 16; ++i) {
        flags[offset+i] += ((val & (1 << i)) >> i);
    }

}

int pospopcnt_u16_avx2_harley_seal_internal(const uint16_t* array, uint32_t len, uint32_t* flags) {
    for (uint32_t i = len - (len % (16 * 16)); i < len; ++i) {
        samtools_single_update(array[i], flags);
    }

    const __m256i* data = (const __m256i*)array;
    size_t size = len / 16;
    __m256i v1  = _mm256_setzero_si256();
    __m256i v2  = _mm256_setzero_si256();
    __m256i v4  = _mm256_setzero_si256();
    __m256i v8  = _mm256_setzero_si256();
    __m256i v16 = _mm256_setzero_si256();
    __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

    __m256i v1U  = _mm256_setzero_si256();
    __m256i v2U  = _mm256_setzero_si256();
    __m256i v4U  = _mm256_setzero_si256();
    __m256i v8U  = _mm256_setzero_si256();
    __m256i v16U = _mm256_setzero_si256();
    __m256i twosAU, twosBU, foursAU, foursBU, eightsAU, eightsBU;

    const uint64_t limit = size - size % 16;
    uint64_t i = 0;
    uint16_t buffer[16];
    __m256i counter[16]; __m256i counterU[16];
    const __m256i one = _mm256_set1_epi16(1);
    const __m256i zero = _mm256_set1_epi16(0);

    const __m256i mask1 = _mm256_set1_epi16(256 + 2048);
    const __m256i mask2 = _mm256_set1_epi16(4 + 256 + 1024 + 2048); // Need to keep 4 + 1024 for Always rule and 256 and 2048 for Rule 1 and Rule 2.
    const __m256i mask3 = _mm256_set1_epi16(512); // QC fail

    while (i < limit) {        
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm256_setzero_si256();
            counterU[i] = _mm256_setzero_si256();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

        // Mask operation: x[i] & ((x[i] & (256 + 2048)) == 0)
        // Need to keep 4 + 1024 for Always rule and 256 and 2048 for Rule 1 and Rule 2.
        // Mask operation: x[i] & (((x[i] & (256 + 2048)) > 0) | (4 + 1024 + 256 + 2048))
        
        //_mm256_loadu_si256(data + i +  0) & (_mm256_cmpeq_epi16( _mm256_loadu_si256(data + i +  0) & mask1, 0 ) | mask2);
        // _mm256_cmpeq_epi16( _mm256_loadu_si256(data + i +  0) & one, one )
#define LOAD(j) __m256i data##j = _mm256_loadu_si256(data + i + j) & (_mm256_cmpeq_epi16( _mm256_loadu_si256(data + i + j) & mask1, zero ) | mask2);
#define L(j) data##j & _mm256_cmpeq_epi16( data##j & mask3, zero )
#define LU(j) data##j & _mm256_cmpeq_epi16( data##j & mask3, mask3 )
        //   _mm256_cmpeq_epi16( data##j & mask3, mask3 )

        // Rule 3.
        // Mask can be written as ((x[i] & 1) == 1)
        // True map to FFFF and False to 0000.
        // Therefore: x[i] & (((x[i] & 1) == 1) is either x[i] or 0.
        
        for (/**/; i < thislimit; i += 16) {
#define U(pos) {                     \
    counter[pos] = _mm256_add_epi16(counter[pos], _mm256_and_si256(v16, one)); \
    v16 = _mm256_srli_epi16(v16, 1); \
}
#define UU(pos) {                     \
    counterU[pos] = _mm256_add_epi16(counterU[pos], _mm256_and_si256(v16U, one)); \
    v16U = _mm256_srli_epi16(v16U, 1); \
}
            LOAD(0) LOAD(1)
            pospopcnt_csa_avx2(&twosA,  &v1, L( 0), L( 1));
            pospopcnt_csa_avx2(&twosAU,  &v1U, LU( 0), LU( 1));
            LOAD(2) LOAD(3)
            pospopcnt_csa_avx2(&twosB,  &v1, L( 2), L( 3));
            pospopcnt_csa_avx2(&twosBU,  &v1U, LU( 2), LU( 3));
            pospopcnt_csa_avx2(&foursA, &v2, twosA, twosB);
            pospopcnt_csa_avx2(&foursAU, &v2U, twosAU, twosBU);
            LOAD(4) LOAD(5)
            pospopcnt_csa_avx2(&twosA,  &v1, L( 4), L( 5));
            pospopcnt_csa_avx2(&twosAU,  &v1U, LU( 4), LU( 5));
            LOAD(6) LOAD(7)
            pospopcnt_csa_avx2(&twosB,  &v1, L( 6), L( 7));
            pospopcnt_csa_avx2(&twosBU,  &v1U, LU( 6), LU( 7));
            pospopcnt_csa_avx2(&foursB, &v2, twosA, twosB);
            pospopcnt_csa_avx2(&foursBU, &v2U, twosAU, twosBU);
            pospopcnt_csa_avx2(&eightsA,&v4, foursA, foursB);
            pospopcnt_csa_avx2(&eightsAU,&v4U, foursAU, foursBU);
            LOAD(8) LOAD(9)
            pospopcnt_csa_avx2(&twosA,  &v1, L( 8),  L( 9));
            pospopcnt_csa_avx2(&twosAU,  &v1U, LU( 8),  LU( 9));
            LOAD(10) LOAD(11)
            pospopcnt_csa_avx2(&twosB,  &v1, L(10),  L(11));
            pospopcnt_csa_avx2(&twosBU,  &v1U, LU(10),  LU(11));
            pospopcnt_csa_avx2(&foursA, &v2, twosA, twosB);
            pospopcnt_csa_avx2(&foursAU, &v2U, twosAU, twosBU);
            LOAD(12) LOAD(13)
            pospopcnt_csa_avx2(&twosA,  &v1, L(12),  L(13));
            pospopcnt_csa_avx2(&twosAU,  &v1U, LU(12),  LU(13));
            LOAD(14) LOAD(15)
            pospopcnt_csa_avx2(&twosB,  &v1, L(14),  L(15));
            pospopcnt_csa_avx2(&twosBU,  &v1U, LU(14),  LU(15));
            pospopcnt_csa_avx2(&foursB, &v2, twosA, twosB);
            pospopcnt_csa_avx2(&foursBU, &v2U, twosAU, twosBU);
            pospopcnt_csa_avx2(&eightsB,&v4, foursA, foursB);
            pospopcnt_csa_avx2(&eightsBU,&v4U, foursAU, foursBU);
            U(0)  U(1)  U(2)  U(3)  U(4)  U(5)  U(6)  U(7)  U(8)  U(9)  U(10)  U(11)  U(12)  U(13)  U(14)  U(15)  // Updates
            UU(0) UU(1) UU(2) UU(3) UU(4) UU(5) UU(6) UU(7) UU(8) UU(9) UU(10) UU(11) UU(12) UU(13) UU(14) UU(15) // Updates
            pospopcnt_csa_avx2(&v16,    &v8, eightsA, eightsB);
            pospopcnt_csa_avx2(&v16U,    &v8U, eightsAU, eightsBU);
#undef U
#undef LOAD
#undef L
        }

        // update the counters after the last iteration
        for (size_t i = 0; i < 16; ++i) {
            counter[i] = _mm256_add_epi16(counter[i], _mm256_and_si256(v16, one));
            v16 = _mm256_srli_epi16(v16, 1);
            counterU[i] = _mm256_add_epi16(counterU[i], _mm256_and_si256(v16U, one));
            v16U = _mm256_srli_epi16(v16U, 1);
        }
        
        for (size_t i = 0; i < 16; ++i) {
            _mm256_storeu_si256((__m256i*)buffer, counter[i]);
            for (size_t z = 0; z < 16; z++) {
                flags[i] += 16 * (uint32_t)buffer[z];
            }

            _mm256_storeu_si256((__m256i*)buffer, counterU[i]);
            for (size_t z = 0; z < 16; z++) {
                flags[16+i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm256_storeu_si256((__m256i*)buffer, v1);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm256_storeu_si256((__m256i*)buffer, v1U);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm256_storeu_si256((__m256i*)buffer, v2);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm256_storeu_si256((__m256i*)buffer, v2U);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm256_storeu_si256((__m256i*)buffer, v4);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm256_storeu_si256((__m256i*)buffer, v4U);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm256_storeu_si256((__m256i*)buffer, v8);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm256_storeu_si256((__m256i*)buffer, v8U);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    return 0;
}

// @see: https://stackoverflow.com/questions/6818606/how-to-programmatically-clear-the-filesystem-memory-cache-in-c-on-a-linux-syst
void clear_cache() {
    sync();
    std::ofstream ofs("/proc/sys/vm/drop_caches");
    ofs << "3" << std::endl;
}
 
int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int ret = ZSTD_decompress(out, out_capacity, in, n_in);
    return(ret);
}

static const std::string SAM_FLAG_NAME[] = {"FPAIRED","FPROPER_PAIR","FUNMAP","FMUNMAP","FREVERSE","FMREVERSE", "FREAD1","FREAD2","FSECONDARY","FQCFAIL","FDUP","FSUPPLEMENTARY"};

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
    if (of.good() == false) return 0;

    uint8_t buffer[1024000]; // 512k 16-bit ints 
    const int max_dst_size = LZ4_compressBound(1024000);
    uint8_t* out_buffer = new uint8_t[max_dst_size];

    while (f.good()) {
        f.read((char*)buffer, 1024000);
        int32_t bytes_read = f.gcount();

        const int32_t compressed_data_size = LZ4_compress_HC((char*)buffer, (char*)out_buffer, bytes_read, max_dst_size, clevel);
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
    
    uint8_t* buffer;
    // assert(!posix_memalign((void**)buffer, SIMD_ALIGNMENT, 1024000));
    buffer = new uint8_t[1024000];
    uint8_t* out_buffer;
    assert(!posix_memalign((void**)&out_buffer, SIMD_ALIGNMENT, 1024000));
    // out_buffer = new uint8_t[1024000];

    int32_t uncompresed_size, compressed_size;

    uint32_t counters[16] = {0}; // flags
    uint64_t tot_flags = 0;

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
    delete[] out_buffer;
    // free(buffer);
    // free(out_buffer);
    return 1;
}

int lz4_decompress(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    // uint8_t buffer[1024000]; // 512k 16-bit ints 
    // uint8_t out_buffer[1024000];
    
    uint8_t* buffer = new uint8_t[1024000];
    uint8_t* out_buffer;
    assert(!posix_memalign((void**)&out_buffer, SIMD_ALIGNMENT, 1024000));
    // out_buffer = new uint8_t[1024000];

    int32_t uncompresed_size, compressed_size;

    uint32_t counters[16*2] = {0}; // flags
    uint64_t tot_flags = 0;

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
        pospopcnt_u16_avx2_harley_seal_internal((uint16_t*)out_buffer,N,counters);

        // std::cerr << "Decompressed " << compressed_size << "->" << uncompresed_size << std::endl;
        if (f.tellg() == filesize) break;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[LZ4 " << file << "] Time elapsed " << time_span.count() << " ms " << tot_flags << std::endl;

    std::cerr << "Tot flags=" << tot_flags << std::endl;
    std::cerr << "Pass QC" << std::endl;
    for (int i = 0; i < 12; ++i) {
        std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[i] << std::endl;
    }

    std::cerr << "Fail QC" << std::endl;
    for (int i = 0; i < 12; ++i) {
        std::cerr << SAM_FLAG_NAME[i] << "\t" << counters[16+i] << std::endl;
    }

    delete[] buffer;
    delete[] out_buffer;
    // free(buffer);
    // free(out_buffer);
    return 1;
}

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
        int w = (c & BAM_FQCFAIL)? 1 : 0;                       \
        ++(s)->n_reads[w];                                      \
        if (c & BAM_FSECONDARY ) {                              \
            ++(s)->n_secondary[w];                              \
        } else if (c & BAM_FSUPPLEMENTARY ) {                   \
            ++(s)->n_supp[w];                                   \
        } else if (c & BAM_FPAIRED) {                           \
            ++(s)->n_pair_all[w];                               \
            if ((c & BAM_FPROPER_PAIR) && !(c & BAM_FUNMAP) ) ++(s)->n_pair_good[w]; \
            if (c & BAM_FREAD1) ++(s)->n_read1[w];              \
            if (c & BAM_FREAD2) ++(s)->n_read2[w];              \
            if ((c & BAM_FMUNMAP) && !(c & BAM_FUNMAP)) ++(s)->n_sgltn[w]; \
            if (!(c & BAM_FUNMAP) && !(c & BAM_FMUNMAP)) {      \
                ++(s)->n_pair_map[w];                           \
            }                                                   \
        }                                                       \
        if (!(c & BAM_FUNMAP)) ++(s)->n_mapped[w];              \
        if (c & BAM_FDUP) ++(s)->n_dup[w];                      \
} while (0)

// Beware data is modified in-place.
static inline
void flagstat_loop_branchless(bam_flagstat_t* s, uint16_t* inflags, const uint32_t i) {
    // QC redirect.
    int w = (inflags[i] & BAM_FQCFAIL) ? 1 : 0;
    
    // Always.
    ++(s)->n_reads[w];
    (s)->n_mapped[w] += (inflags[i] & BAM_FUNMAP) == 0; // this is implicit
    (s)->n_dup[w] += (inflags[i] & BAM_FDUP) >> 10;

    // Rule 1.
    (s)->n_secondary[w] += (inflags[i] & BAM_FSECONDARY) >> 8;
    inflags[i] &= lookup_sec[(inflags[i] & BAM_FSECONDARY) >> 8];
    
    // Rule 2.
    (s)->n_supp[w] += (inflags[i] & BAM_FSUPPLEMENTARY) >> 11;
    inflags[i] &= lookup_sup[(inflags[i] & BAM_FSUPPLEMENTARY) >> 11];
    // Mask operation: x[i] & ((x[i] & (256 + 2048)) > 0)
    // Need to keep 4 + 1024 for Always rule and 256 and 2048 for Rule 1 and Rule 2.
    // Mask operation: x[i] & (((x[i] & (256 + 2048)) > 0) | (4 + 1024 + 256 + 2048))
    
    // Rule 3.
    // Mask can be written as ((x[i] & 1) == 1)
    // True map to FFFF and False to 0000.
    // Therefore: x[i] & (((x[i] & 1) == 1) is either x[i] or 0.
    inflags[i] &= lookup_pair[inflags[i] & BAM_FPAIRED];
    (s)->n_pair_all[w]  += (inflags[i] & BAM_FPAIRED);
    (s)->n_pair_map[w]  += ((inflags[i] & (BAM_FMUNMAP + BAM_FUNMAP + BAM_FSECONDARY + BAM_FSUPPLEMENTARY)) == 0); // 8 + 4 + 256 + 2048
    (s)->n_pair_good[w] += ((inflags[i] & (BAM_FPROPER_PAIR + BAM_FUNMAP)) == BAM_FPROPER_PAIR);
    (s)->n_read1[w] += (inflags[i] & BAM_FREAD1) >> 6;
    (s)->n_read2[w] += (inflags[i] & BAM_FREAD2) >> 7;
    (s)->n_sgltn[w] += ((inflags[i] & (BAM_FMUNMAP + BAM_FUNMAP)) == BAM_FMUNMAP);
}

static const char* percent(char* buffer, long long n, long long total)
{
    if (total != 0) sprintf(buffer, "%.2f%%", (float)n / total * 100.0);
    else strcpy(buffer, "N/A");
    return buffer;
}

int lz4_decompress_samtools(const std::string& file) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) return 0;
    int64_t filesize = f.tellg();
    f.seekg(0);
    uint8_t buffer[1024000]; // 512k 16-bit ints 
    uint8_t out_buffer[1024000];

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

        // assert(decompressed_size == uncompresed_size);

        const uint32_t N = uncompresed_size >> 1;
        tot_flags += N;

        uint16_t* inflags = (uint16_t*)out_buffer;
        // Rules:
        // Always: !4,1024
        // Rule1: 256 only
        // Rule2: 2048 only
        // Rule3: 1
        //        2 + 4 == 2
        //        64
        //        128
        //        8 + 4 == 8
        //        8 + 4 == 0
        for (int i = 0; i < N; ++i) {
                flagstat_loop(s, inflags[i]);
                // flagstat_loop_branchless(s, inflags, i);
            
            // if (1) {
            //     int w = (inflags[i] & BAM_FQCFAIL) ? 1 : 0;
            //     ++(s)->n_reads[w];
            //     if (inflags[i] & BAM_FSECONDARY) { // If secondary then count nothing
            //         ++(s)->n_secondary[w];
            //     } else if (inflags[i] & BAM_FSUPPLEMENTARY) { // If supplementary then count nothing
            //         ++(s)->n_supp[w];
            //     } else if (inflags[i] & BAM_FPAIRED) { // If paired then count the paired information
            //         ++(s)->n_pair_all[w];

            //         // If proper pair and not unmapped.
            //         if ((inflags[i] & (BAM_FPROPER_PAIR + BAM_FUNMAP)) == BAM_FPROPER_PAIR) ++(s)->n_pair_good[w];
            //         // if ((inflags[i] & BAM_FPROPER_PAIR) && !(inflags[i] & BAM_FUNMAP) ) ++(s)->n_pair_good[w];
            //         // Read 1.
            //         if (inflags[i] & BAM_FREAD1) ++(s)->n_read1[w];
            //         // Read 2.
            //         if (inflags[i] & BAM_FREAD2) ++(s)->n_read2[w];
            //         // Mate is unmapped but read is mapped.
            //         if ((inflags[i] & (BAM_FMUNMAP + BAM_FUNMAP)) == BAM_FMUNMAP) ++(s)->n_sgltn[w];
            //         // (s)->n_sgltn[w] += !!((inflags[i] & (BAM_FMUNMAP + BAM_FUNMAP)) == BAM_FMUNMAP);

            //         // if ((inflags[i] & BAM_FMUNMAP) && !(inflags[i] & BAM_FUNMAP))  ++(s)->n_sgltn[w];
            //         // Mate is mapped and read is mapped.
            //         // if (!(inflags[i] & (BAM_FMUNMAP + BAM_FUNMAP)) == 0) ++(s)->n_pair_map[w];
            //         if (!(inflags[i] & BAM_FUNMAP) && !(inflags[i] & BAM_FMUNMAP)) ++(s)->n_pair_map[w];
            //     }
            //     if (!(inflags[i] & BAM_FUNMAP)) ++(s)->n_mapped[w];
            //     if (inflags[i] & BAM_FDUP) ++(s)->n_dup[w];
            // }
            
            /*
            // QC redirect.
            int w = (inflags[i] & BAM_FQCFAIL) ? 1 : 0;
            
            // Always.
            ++(s)->n_reads[w];
            (s)->n_mapped[w] += (inflags[i] & BAM_FUNMAP) == 0; // this is implicit
            (s)->n_dup[w] += (inflags[i] & BAM_FDUP) >> 10;

            // Rule 1.
            (s)->n_secondary[w] += (inflags[i] & BAM_FSECONDARY) >> 8;
            inflags[i] &= lookup_sec[(inflags[i] & BAM_FSECONDARY) >> 8];
            
            // Rule 2.
            (s)->n_supp[w] += (inflags[i] & BAM_FSUPPLEMENTARY) >> 11;
            inflags[i] &= lookup_sup[(inflags[i] & BAM_FSUPPLEMENTARY) >> 11];
            // Mask operation: x[i] & ((x[i] & (256 + 2048)) > 0)
            // Need to keep 4 + 1024 for Always rule and 256 and 2048 for Rule 1 and Rule 2.
            // Mask operation: x[i] & (((x[i] & (256 + 2048)) > 0) | (4 + 1024 + 256 + 2048))
            
            // Rule 3.
            // Mask can be written as ((x[i] & 1) == 1)
            // True map to FFFF and False to 0000.
            // Therefore: x[i] & (((x[i] & 1) == 1) is either x[i] or 0.
            inflags[i] &= lookup_pair[inflags[i] & BAM_FPAIRED];
            (s)->n_pair_all[w]  += (inflags[i] & BAM_FPAIRED);
            (s)->n_pair_map[w]  += ((inflags[i] & (BAM_FMUNMAP + BAM_FUNMAP + BAM_FSECONDARY + BAM_FSUPPLEMENTARY)) == 0); // 8 + 4 + 256 + 2048
            (s)->n_pair_good[w] += ((inflags[i] & (BAM_FPROPER_PAIR + BAM_FUNMAP)) == BAM_FPROPER_PAIR);
            (s)->n_read1[w] += (inflags[i] & BAM_FREAD1) >> 6;
            (s)->n_read2[w] += (inflags[i] & BAM_FREAD2) >> 7;
            (s)->n_sgltn[w] += ((inflags[i] & (BAM_FMUNMAP + BAM_FUNMAP)) == BAM_FMUNMAP);
            */
        }

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

    // s = bam_flagstat_core(fp, header);
    if (1) {
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
    }
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

    uint32_t counters[16] = {0}; // flags
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
        pospopcnt_u16((uint16_t*)out_buffer,N,counters);

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
            // flagstat_loop(s, inflags[i]);
            flagstat_loop_branchless(s, inflags, i);
        }

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
        if (lz4_fast) lz4f(input,output, clevel);
        else lz4hc(input,output,clevel);
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
		{0,0,0,0}
	};
	
    std::string input;
    bool lz4 = false;
    bool mzstd = false;

	while ((c = getopt_long(argc, argv, "i:?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			input = std::string(optarg);
			break;
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
        clear_cache();
        zstd_decompress_only(input); // warmup
        clear_cache();
        zstd_decompress_only(input);
        clear_cache();
        zstd_decompress_samtools(input);
        clear_cache();
        zstd_decompress(input);
    }
    
    if (method == 2) {
        clear_cache();
        lz4_decompress_only(input); // warmup
        clear_cache();
        lz4_decompress_only(input);
        clear_cache();
        lz4_decompress_samtools(input);
        clear_cache();
        lz4_decompress(input);
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