// License for libflagstats.h
/*
* Copyright (c) 2019 Marcus D. R. Klarqvist
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
// License for libalgebra.h
/*
* Copyright (c) 2019 Marcus D. R. Klarqvist
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
// License for pospopcnt.h
/*
* Copyright (c) 2019
* Author(s): Marcus D. R. Klarqvist, Wojciech Muła, and Daniel Lemire
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
#ifndef LIBFLAGSTATS_H_98345934843
#define LIBFLAGSTATS_H_98345934843

/* *************************************
*  Includes
***************************************/
#include "libalgebra.h"

/* *************************************
*  Target SAM fields
***************************************/

// Modified from sam.h
/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define FLAGSTAT_FPAIRED               1
#define FLAGSTAT_FPAIRED_OFF           0
/*! @abstract the read is mapped in a proper pair */
#define FLAGSTAT_FPROPER_PAIR          2
#define FLAGSTAT_FPROPER_PAIR_OFF      1
/*! @abstract the read itself is unmapped; conflictive with FLAGSTAT_FPROPER_PAIR */
#define FLAGSTAT_FUNMAP                4
#define FLAGSTAT_FUNMAP_OFF            2
/*! @abstract the mate is unmapped */
#define FLAGSTAT_FMUNMAP               8
#define FLAGSTAT_FMUNMAP_OFF           3
/*! @abstract the read is mapped to the reverse strand */
#define FLAGSTAT_FREVERSE             16
#define FLAGSTAT_FREVERSE_OFF          4
/*! @abstract the mate is mapped to the reverse strand */
#define FLAGSTAT_FMREVERSE            32
#define FLAGSTAT_FMREVERSE_OFF         5
/*! @abstract this is read1 */
#define FLAGSTAT_FREAD1               64
#define FLAGSTAT_FREAD1_OFF            6
/*! @abstract this is read2 */
#define FLAGSTAT_FREAD2              128
#define FLAGSTAT_FREAD2_OFF            7
/*! @abstract not primary alignment */
#define FLAGSTAT_FSECONDARY          256
#define FLAGSTAT_FSECONDARY_OFF        8
/*! @abstract QC failure */
#define FLAGSTAT_FQCFAIL             512
#define FLAGSTAT_FQCFAIL_OFF           9
/*! @abstract optical or PCR duplicate */
#define FLAGSTAT_FDUP               1024
#define FLAGSTAT_FDUP_OFF             10
/*! @abstract supplementary alignment */
#define FLAGSTAT_FSUPPLEMENTARY     2048
#define FLAGSTAT_FSUPPLEMENTARY_OFF   11

#ifdef __cplusplus
extern "C" {
#endif

void FLAGSTAT_scalar_update(uint16_t val, uint32_t* flags) {
    // If the FLAGSTAT_FQCFAIL is set the data is shift 16 values to
    // the right to distinguish between statistics for data
    // that failed and passed quality control.
    const int offset = ( (val & FLAGSTAT_FQCFAIL) == 0 ) ? 0 : 16;
    uint32_t* f = &flags[offset];
    // Count only reads that with FLAGSTAT_FQCFAIL set. The other
    // reads are implicitly known and computed at the end of
    // FLAGSTAT_* functions.
    if (offset) ++f[FLAGSTAT_FQCFAIL_OFF];

    if (val & FLAGSTAT_FSECONDARY) ++f[FLAGSTAT_FSECONDARY_OFF];
    else if (val & FLAGSTAT_FSUPPLEMENTARY) ++f[FLAGSTAT_FSUPPLEMENTARY_OFF];
    else if (val & FLAGSTAT_FPAIRED) {
        // ++(s)->n_pair_all[w];              
        if ( (val & FLAGSTAT_FPROPER_PAIR) && !(val & FLAGSTAT_FUNMAP) ) ++f[12];
        if (val & FLAGSTAT_FREAD1) ++f[FLAGSTAT_FREAD1_OFF];
        if (val & FLAGSTAT_FREAD2) ++f[FLAGSTAT_FREAD2_OFF];
        if ((val & FLAGSTAT_FMUNMAP) && !(val & FLAGSTAT_FUNMAP))  ++f[13];
        if (!(val & FLAGSTAT_FUNMAP) && !(val & FLAGSTAT_FMUNMAP)) ++f[14];
    }
    // Count as is FUNMAP then use arithmetic to compute N - FUNMAP
    if (val & FLAGSTAT_FUNMAP) ++f[FLAGSTAT_FUNMAP_OFF];
    if (val & FLAGSTAT_FDUP)      ++f[FLAGSTAT_FDUP_OFF];
}

// #define SAMTOOLS_flagstat_loop(s, c) do {                   \
//     int w = (c & FLAGSTAT_FQCFAIL)? 1 : 0;                  \
//     ++(s)->n_reads[w];                                      \
//     if (c & FLAGSTAT_FSECONDARY ) {                         \
//         ++(s)->n_secondary[w];                              \
//     } else if (c & FLAGSTAT_FSUPPLEMENTARY ) {              \
//         ++(s)->n_supp[w];                                   \
//     } else if (c & FLAGSTAT_FPAIRED) {                      \
//         ++(s)->n_pair_all[w];                               \
//         if ( (c & FLAGSTAT_FPROPER_PAIR) && !(c & FLAGSTAT_FUNMAP) ) ++(s)->n_pair_good[w]; \
//         if (c & FLAGSTAT_FREAD1) ++(s)->n_read1[w];         \
//         if (c & FLAGSTAT_FREAD2) ++(s)->n_read2[w];         \
//         if ((c & FLAGSTAT_FMUNMAP) && !(c & FLAGSTAT_FUNMAP)) ++(s)->n_sgltn[w]; \
//         if (!(c & FLAGSTAT_FUNMAP) && !(c & FLAGSTAT_FMUNMAP)) {      \
//             ++(s)->n_pair_map[w];                           \
//         }                                                   \
//     }                                                       \
//     if (!(c & FLAGSTAT_FUNMAP)) ++(s)->n_mapped[w];         \
//     if (c & FLAGSTAT_FDUP) ++(s)->n_dup[w];                 \
// } while (0)

// FLAGSTAT_FPROPER_PAIR & !FLAGSTAT_FUNMAP
// x |= (x & (FLAGSTAT_FPROPER_PAIR + FLAGSTAT_FUNMAP) == FLAGSTAT_FPROPER_PAIR) & 1 << 13
// FLAGSTAT_FMUNMAP & !FLAGSTAT_FUNMAP
// x |= (x & (FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP) == FLAGSTAT_FMUNMAP) & 1 << 14

static
int FLAGSTAT_scalar(const uint16_t* array, uint32_t len, uint32_t* flags) {
    for (uint32_t i = 0; i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], flags);
    }
}

#if defined(STORM_HAVE_SSE42)

#include <immintrin.h>

STORM_TARGET("sse4.2")
static
int FLAGSTAT_sse4(const uint16_t* array, uint32_t len, uint32_t* flags) {
    const uint32_t start_qc = flags[FLAGSTAT_FQCFAIL_OFF + 16];
    
    for (uint32_t i = len - (len % (16 * 8)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], flags);
    }

    const __m128i* data = (const __m128i*)array;
    size_t size = len / 8;
    __m128i v1  = _mm_setzero_si128();
    __m128i v2  = _mm_setzero_si128();
    __m128i v4  = _mm_setzero_si128();
    __m128i v8  = _mm_setzero_si128();
    __m128i v16 = _mm_setzero_si128();
    __m128i twosA, twosB, foursA, foursB, eightsA, eightsB;

    __m128i v1U  = _mm_setzero_si128();
    __m128i v2U  = _mm_setzero_si128();
    __m128i v4U  = _mm_setzero_si128();
    __m128i v8U  = _mm_setzero_si128();
    __m128i v16U = _mm_setzero_si128();
    __m128i twosAU, twosBU, foursAU, foursBU, eightsAU, eightsBU;

    const uint64_t limit = size - size % 16;
    uint64_t i = 0;
    uint16_t buffer[8];
    __m128i counter[16];
    __m128i counterU[16];
    
    // Masks and mask selectors.
    const __m128i m1   = _mm_set1_epi16(FLAGSTAT_FSECONDARY);
    const __m128i m1S  = _mm_set1_epi16(FLAGSTAT_FQCFAIL + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP);
    const __m128i m2   = _mm_set1_epi16(FLAGSTAT_FSUPPLEMENTARY);
    const __m128i m2S  = _mm_set1_epi16(FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP);
    const __m128i m3   = _mm_set1_epi16(FLAGSTAT_FPAIRED);
    const __m128i m4   = _mm_set1_epi16(FLAGSTAT_FQCFAIL);
    const __m128i one  = _mm_set1_epi16(1); // (00...1) vector
    const __m128i zero = _mm_set1_epi16(0); // (00...0) vector

    // Main body.
    while (i < limit) {        
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm_setzero_si128();
            counterU[i] = _mm_setzero_si128();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

        ///////////////////////////////////////////////////////////////////////
        // We load a register of data (data + i + j) and then construct the
        // conditional bits: 
        // 12: FLAGSTAT_FPROPER_PAIR + FLAGSTAT_FUNMAP == FLAGSTAT_FPROPER_PAIR
        // 13: FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP == FLAGSTAT_FMUNMAP
        // 14: FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP == 0
        //
        // These construction of these bits can be described for data x as:
        // x |= (x & LEFT_MASK == RIGHT_MASK) & 1 << TARGET_BIT
        // with the assumption that predicate evaluatons result in the selection
        // masks (00...0) or (11...1) for FALSE and TRUE, respectively. These
        // construction macros are named O1, O2, and O3.
        //
        // The original SAMtools method is also heavily branched with three
        // main branch points:
        // If FLAGSTAT_FSECONDARY then count FLAGSTAT_FSECONDARY
        // If FLAGSTAT_FSUPPLEMENTARY then count FLAGSTAT_FSUPPLEMENTARY
        // Else then count FLAGSTAT_FREAD1, 
        //                 FLAGSTAT_FREAD2,
        //                 Special bit 12, 13, and 14
        // Always count FLAGSTAT_FUNMAP, 
        //              FLAGSTAT_FDUP, 
        //              FLAGSTAT_FQCFAIL
        //
        // These bits can be selected using a mask-select propagate-carry approach:
        // x &= x & ((x == MASK) | CARRY_BITS)
        // with the arguments for MASK and CARRY_BITS as follows:
        //    1. {FLAGSTAT_FSECONDARY, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //    2. {FLAGSTAT_FSUPPLEMENTARY, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY 
        //        + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //    3. {FLAGSTAT_FPAIRED, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY 
        //        + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //
        // FLAGSTATS outputs summary statistics separately for reads that pass
        // QC and those that do not. Therefore we need to partition the data
        // into these two classes. For data that pass QC, the L registers, we
        // first bit-select the target FLAGSTAT_FQCFAIL bit using the mask
        // mask3. The resulting data is used to perform another mask-select
        // using VPCMPEQW against the empty vector (00...0). As above, if the
        // data has the FLAGSTAT_FQCFAIL bit set then this register will be
        // zeroed out. The exact process is performed for reads that fail QC,
        // the LU registers, with the difference that mask-selection is based on
        // the one vector (00...1).

#define W(j) __m128i data##j = _mm_loadu_si128(data + i + j);
#define O1(j) data##j = data##j | _mm_slli_epi16(data##j & _mm_cmpeq_epi16((data##j & _mm_set1_epi16(FLAGSTAT_FPROPER_PAIR + FLAGSTAT_FUNMAP)), _mm_set1_epi16(FLAGSTAT_FPROPER_PAIR)) & one, 12); 
#define O2(j) data##j = data##j | _mm_slli_epi16(data##j & _mm_cmpeq_epi16((data##j & _mm_set1_epi16(FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP)), _mm_set1_epi16(FLAGSTAT_FMUNMAP)) & one, 13); 
#define O3(j) data##j = data##j | _mm_slli_epi16(data##j & _mm_cmpeq_epi16((data##j & _mm_set1_epi16(FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP)), zero) & one, 14); 
#define L1(j) data##j = data##j & (_mm_cmpeq_epi16((data##j & m1), zero) | m1S);
#define L2(j) data##j = data##j & (_mm_cmpeq_epi16((data##j & m2), zero) | m2S);
#define L3(j) data##j = data##j & (_mm_cmpeq_epi16((data##j & m3), m3)   | m2S);
#define LOAD(j) W(j) O1(j) O2(j) O3(j) L1(j) L2(j) L3(j)
#define L(j)  data##j & _mm_cmpeq_epi16( data##j & m4, zero )
#define LU(j) data##j & _mm_cmpeq_epi16( data##j & m4, m4 )
        
        for (/**/; i < thislimit; i += 16) {
#define U(pos) {                     \
    counter[pos] = _mm_add_epi16(counter[pos], _mm_and_si128(v16, one));    \
    v16 = _mm_srli_epi16(v16, 1); \
}
#define UU(pos) {                      \
    counterU[pos] = _mm_add_epi16(counterU[pos], _mm_and_si128(v16U, one)); \
    v16U = _mm_srli_epi16(v16U, 1); \
}
            LOAD(0) LOAD(1)
            STORM_pospopcnt_csa_sse(&twosA,   &v1,  L( 0),  L( 1));
            STORM_pospopcnt_csa_sse(&twosAU,  &v1U, LU( 0), LU( 1));
            LOAD(2) LOAD(3)
            STORM_pospopcnt_csa_sse(&twosB,   &v1,  L( 2),  L( 3));
            STORM_pospopcnt_csa_sse(&twosBU,  &v1U, LU( 2), LU( 3));
            STORM_pospopcnt_csa_sse(&foursA,  &v2,  twosA, twosB);
            STORM_pospopcnt_csa_sse(&foursAU, &v2U, twosAU, twosBU);
            LOAD(4) LOAD(5)
            STORM_pospopcnt_csa_sse(&twosA,   &v1,  L( 4),  L( 5));
            STORM_pospopcnt_csa_sse(&twosAU,  &v1U, LU( 4), LU( 5));
            LOAD(6) LOAD(7)
            STORM_pospopcnt_csa_sse(&twosB,   &v1,  L( 6),  L( 7));
            STORM_pospopcnt_csa_sse(&twosBU,  &v1U, LU( 6), LU( 7));
            STORM_pospopcnt_csa_sse(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_sse(&foursBU, &v2U, twosAU,  twosBU);
            STORM_pospopcnt_csa_sse(&eightsA, &v4,  foursA,  foursB);
            STORM_pospopcnt_csa_sse(&eightsAU,&v4U, foursAU, foursBU);
            LOAD(8) LOAD(9)
            STORM_pospopcnt_csa_sse(&twosA,   &v1,  L( 8),   L( 9));
            STORM_pospopcnt_csa_sse(&twosAU,  &v1U, LU( 8),  LU( 9));
            LOAD(10) LOAD(11)
            STORM_pospopcnt_csa_sse(&twosB,   &v1,  L(10),   L(11));
            STORM_pospopcnt_csa_sse(&twosBU,  &v1U, LU(10),  LU(11));
            STORM_pospopcnt_csa_sse(&foursA,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_sse(&foursAU, &v2U, twosAU,  twosBU);
            LOAD(12) LOAD(13)
            STORM_pospopcnt_csa_sse(&twosA,   &v1,  L(12),   L(13));
            STORM_pospopcnt_csa_sse(&twosAU,  &v1U, LU(12),  LU(13));
            LOAD(14) LOAD(15)
            STORM_pospopcnt_csa_sse(&twosB,   &v1,  L(14),   L(15));
            STORM_pospopcnt_csa_sse(&twosBU,  &v1U, LU(14),  LU(15));
            STORM_pospopcnt_csa_sse(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_sse(&foursBU, &v2U, twosAU,  twosBU);
            STORM_pospopcnt_csa_sse(&eightsB, &v4,  foursA,  foursB);
            STORM_pospopcnt_csa_sse(&eightsBU,&v4U, foursAU, foursBU);
             U(0)  U(1)  U(2)  U(3)  U(4)  U(5)  U(6)  U(7)  U(8)  U(9)  U(10)  U(11)  U(12)  U(13)  U(14)  U(15) // Updates
            UU(0) UU(1) UU(2) UU(3) UU(4) UU(5) UU(6) UU(7) UU(8) UU(9) UU(10) UU(11) UU(12) UU(13) UU(14) UU(15) // Updates
            STORM_pospopcnt_csa_sse(&v16,     &v8,  eightsA,  eightsB);
            STORM_pospopcnt_csa_sse(&v16U,    &v8U, eightsAU, eightsBU);
#undef U
#undef UU
#undef LOAD
#undef L
#undef LU
#undef W
#undef O1
#undef O2
#undef O3
#undef L1
#undef L2
#undef L3
        }

        // Update the counters after the last iteration
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm_add_epi16(counter[i], _mm_and_si128(v16, one));
            v16  = _mm_srli_epi16(v16, 1);
            counterU[i] = _mm_add_epi16(counterU[i], _mm_and_si128(v16U, one));
            v16U = _mm_srli_epi16(v16U, 1);
        }
        
        for (size_t i = 0; i < 16; ++i) {
            _mm_storeu_si128((__m128i*)buffer, counter[i]);
            for (size_t z = 0; z < 8; z++) {
                flags[i] += 16 * (uint32_t)buffer[z];
            }

            _mm_storeu_si128((__m128i*)buffer, counterU[i]);
            for (size_t z = 0; z < 8; z++) {
                flags[16+i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm_storeu_si128((__m128i*)buffer, v1);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm_storeu_si128((__m128i*)buffer, v1U);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm_storeu_si128((__m128i*)buffer, v2);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm_storeu_si128((__m128i*)buffer, v2U);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm_storeu_si128((__m128i*)buffer, v4);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm_storeu_si128((__m128i*)buffer, v4U);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm_storeu_si128((__m128i*)buffer, v8);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm_storeu_si128((__m128i*)buffer, v8U);
    for (size_t i = 0; i < 8; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    // QC
    flags[FLAGSTAT_FQCFAIL_OFF] += len - (flags[FLAGSTAT_FQCFAIL_OFF+16] - start_qc);

    return 0;
}

#endif

#if defined(STORM_HAVE_AVX2)

#include <immintrin.h>

STORM_TARGET("avx2")
static
int FLAGSTAT_avx2(const uint16_t* array, uint32_t len, uint32_t* flags) {
    const uint32_t start_qc = flags[FLAGSTAT_FQCFAIL_OFF + 16];
    
    for (uint32_t i = len - (len % (16 * 16)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], flags);
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
    __m256i counter[16]; 
    __m256i counterU[16];
    
    // Masks and mask selectors.
    const __m256i m1   = _mm256_set1_epi16(FLAGSTAT_FSECONDARY);
    const __m256i m1S  = _mm256_set1_epi16(FLAGSTAT_FQCFAIL + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP);
    const __m256i m2   = _mm256_set1_epi16(FLAGSTAT_FSUPPLEMENTARY);
    const __m256i m2S  = _mm256_set1_epi16(FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP);
    const __m256i m3   = _mm256_set1_epi16(FLAGSTAT_FPAIRED);
    const __m256i m4   = _mm256_set1_epi16(FLAGSTAT_FQCFAIL);
    const __m256i one  = _mm256_set1_epi16(1); // (00...1) vector
    const __m256i zero = _mm256_set1_epi16(0); // (00...0) vector

    // Main body.
    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm256_setzero_si256();
            counterU[i] = _mm256_setzero_si256();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

        ///////////////////////////////////////////////////////////////////////
        // We load a register of data (data + i + j) and then construct the
        // conditional bits: 
        // 12: FLAGSTAT_FPROPER_PAIR + FLAGSTAT_FUNMAP == FLAGSTAT_FPROPER_PAIR
        // 13: FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP == FLAGSTAT_FMUNMAP
        // 14: FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP == 0
        //
        // These construction of these bits can be described for data x as:
        // x |= (x & LEFT_MASK == RIGHT_MASK) & 1 << TARGET_BIT
        // with the assumption that predicate evaluatons result in the selection
        // masks (00...0) or (11...1) for FALSE and TRUE, respectively. These
        // construction macros are named O1, O2, and O3.
        //
        // The original SAMtools method is also heavily branched with three
        // main branch points:
        // If FLAGSTAT_FSECONDARY then count FLAGSTAT_FSECONDARY
        // If FLAGSTAT_FSUPPLEMENTARY then count FLAGSTAT_FSUPPLEMENTARY
        // Else then count FLAGSTAT_FREAD1, 
        //                 FLAGSTAT_FREAD2,
        //                 Special bit 12, 13, and 14
        // Always count FLAGSTAT_FUNMAP, 
        //              FLAGSTAT_FDUP, 
        //              FLAGSTAT_FQCFAIL
        //
        // These bits can be selected using a mask-select propagate-carry approach:
        // x &= x & ((x == MASK) | CARRY_BITS)
        // with the arguments for MASK and CARRY_BITS as follows:
        //    1. {FLAGSTAT_FSECONDARY, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //    2. {FLAGSTAT_FSUPPLEMENTARY, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY 
        //        + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //    3. {FLAGSTAT_FPAIRED, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY 
        //        + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //
        // FLAGSTATS outputs summary statistics separately for reads that pass
        // QC and those that do not. Therefore we need to partition the data
        // into these two classes. For data that pass QC, the L registers, we
        // first bit-select the target FLAGSTAT_FQCFAIL bit using the mask
        // mask3. The resulting data is used to perform another mask-select
        // using VPCMPEQW against the empty vector (00...0). As above, if the
        // data has the FLAGSTAT_FQCFAIL bit set then this register will be
        // zeroed out. The exact process is performed for reads that fail QC,
        // the LU registers, with the difference that mask-selection is based on
        // the one vector (00...1).

#define W(j) __m256i data##j = _mm256_loadu_si256(data + i + j);
#define O1(j) data##j = data##j | _mm256_slli_epi16(data##j & _mm256_cmpeq_epi16((data##j & _mm256_set1_epi16(FLAGSTAT_FPROPER_PAIR + FLAGSTAT_FUNMAP)), _mm256_set1_epi16(FLAGSTAT_FPROPER_PAIR)) & one, 12); 
#define O2(j) data##j = data##j | _mm256_slli_epi16(data##j & _mm256_cmpeq_epi16((data##j & _mm256_set1_epi16(FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP)), _mm256_set1_epi16(FLAGSTAT_FMUNMAP)) & one, 13); 
#define O3(j) data##j = data##j | _mm256_slli_epi16(data##j & _mm256_cmpeq_epi16((data##j & _mm256_set1_epi16(FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP)), zero) & one, 14); 
#define L1(j) data##j = data##j & (_mm256_cmpeq_epi16((data##j & m1), zero) | m1S);
#define L2(j) data##j = data##j & (_mm256_cmpeq_epi16((data##j & m2), zero) | m2S);
#define L3(j) data##j = data##j & (_mm256_cmpeq_epi16((data##j & m3), m3)   | m2S);
#define LOAD(j) W(j) O1(j) O2(j) O3(j) L1(j) L2(j) L3(j)
#define L(j)  data##j & _mm256_cmpeq_epi16( data##j & m4, zero )
#define LU(j) data##j & _mm256_cmpeq_epi16( data##j & m4, m4 )
        
        for (/**/; i < thislimit; i += 16) {
#define U(pos) {                     \
    counter[pos] = _mm256_add_epi16(counter[pos], _mm256_and_si256(v16, one));    \
    v16 = _mm256_srli_epi16(v16, 1); \
}
#define UU(pos) {                      \
    counterU[pos] = _mm256_add_epi16(counterU[pos], _mm256_and_si256(v16U, one)); \
    v16U = _mm256_srli_epi16(v16U, 1); \
}
            LOAD(0) LOAD(1)
            STORM_pospopcnt_csa_avx2(&twosA,   &v1,  L( 0),  L( 1));
            STORM_pospopcnt_csa_avx2(&twosAU,  &v1U, LU( 0), LU( 1));
            LOAD(2) LOAD(3)
            STORM_pospopcnt_csa_avx2(&twosB,   &v1,  L( 2),  L( 3));
            STORM_pospopcnt_csa_avx2(&twosBU,  &v1U, LU( 2), LU( 3));
            STORM_pospopcnt_csa_avx2(&foursA,  &v2,  twosA, twosB);
            STORM_pospopcnt_csa_avx2(&foursAU, &v2U, twosAU, twosBU);
            LOAD(4) LOAD(5)
            STORM_pospopcnt_csa_avx2(&twosA,   &v1,  L( 4),  L( 5));
            STORM_pospopcnt_csa_avx2(&twosAU,  &v1U, LU( 4), LU( 5));
            LOAD(6) LOAD(7)
            STORM_pospopcnt_csa_avx2(&twosB,   &v1,  L( 6),  L( 7));
            STORM_pospopcnt_csa_avx2(&twosBU,  &v1U, LU( 6), LU( 7));
            STORM_pospopcnt_csa_avx2(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx2(&foursBU, &v2U, twosAU,  twosBU);
            STORM_pospopcnt_csa_avx2(&eightsA, &v4,  foursA,  foursB);
            STORM_pospopcnt_csa_avx2(&eightsAU,&v4U, foursAU, foursBU);
            LOAD(8) LOAD(9)
            STORM_pospopcnt_csa_avx2(&twosA,   &v1,  L( 8),   L( 9));
            STORM_pospopcnt_csa_avx2(&twosAU,  &v1U, LU( 8),  LU( 9));
            LOAD(10) LOAD(11)
            STORM_pospopcnt_csa_avx2(&twosB,   &v1,  L(10),   L(11));
            STORM_pospopcnt_csa_avx2(&twosBU,  &v1U, LU(10),  LU(11));
            STORM_pospopcnt_csa_avx2(&foursA,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx2(&foursAU, &v2U, twosAU,  twosBU);
            LOAD(12) LOAD(13)
            STORM_pospopcnt_csa_avx2(&twosA,   &v1,  L(12),   L(13));
            STORM_pospopcnt_csa_avx2(&twosAU,  &v1U, LU(12),  LU(13));
            LOAD(14) LOAD(15)
            STORM_pospopcnt_csa_avx2(&twosB,   &v1,  L(14),   L(15));
            STORM_pospopcnt_csa_avx2(&twosBU,  &v1U, LU(14),  LU(15));
            STORM_pospopcnt_csa_avx2(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx2(&foursBU, &v2U, twosAU,  twosBU);
            STORM_pospopcnt_csa_avx2(&eightsB, &v4,  foursA,  foursB);
            STORM_pospopcnt_csa_avx2(&eightsBU,&v4U, foursAU, foursBU);
             U(0)  U(1)  U(2)  U(3)  U(4)  U(5)  U(6)  U(7)  U(8)  U(9)  U(10)  U(11)  U(12)  U(13)  U(14)  U(15) // Updates
            UU(0) UU(1) UU(2) UU(3) UU(4) UU(5) UU(6) UU(7) UU(8) UU(9) UU(10) UU(11) UU(12) UU(13) UU(14) UU(15) // Updates
            STORM_pospopcnt_csa_avx2(&v16,     &v8,  eightsA,  eightsB);
            STORM_pospopcnt_csa_avx2(&v16U,    &v8U, eightsAU, eightsBU);
#undef U
#undef UU
#undef LOAD
#undef L
#undef LU
#undef W
#undef O1
#undef O2
#undef O3
#undef L1
#undef L2
#undef L3
        }

        // Update the counters after the last iteration
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm256_add_epi16(counter[i], _mm256_and_si256(v16, one));
            v16  = _mm256_srli_epi16(v16, 1);
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

    // QC
    flags[FLAGSTAT_FQCFAIL_OFF] += len - (flags[FLAGSTAT_FQCFAIL_OFF+16] - start_qc);

    return 0;
}
#endif // end avx2

#if defined(STORM_HAVE_AVX512)

#include <immintrin.h>

STORM_TARGET("avx512bw")
static
int FLAGSTAT_avx512(const uint16_t* array, size_t len, uint32_t* out) {
    for (uint32_t i = len - (len % (32 * 16)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], out);
    }

    const __m512i* data = (const __m512i*)array;
    size_t size = len / 32;
    __m512i v1  = _mm512_setzero_si512();
    __m512i v2  = _mm512_setzero_si512();
    __m512i v4  = _mm512_setzero_si512();
    __m512i v8  = _mm512_setzero_si512();
    __m512i v16 = _mm512_setzero_si512();
    __m512i twosA, twosB, foursA, foursB, eightsA, eightsB;

    __m512i v1U  = _mm512_setzero_si512();
    __m512i v2U  = _mm512_setzero_si512();
    __m512i v4U  = _mm512_setzero_si512();
    __m512i v8U  = _mm512_setzero_si512();
    __m512i v16U = _mm512_setzero_si512();
    __m512i twosAU, twosBU, foursAU, foursBU, eightsAU, eightsBU;

    const uint64_t limit = size - size % 16;
    uint64_t i = 0;
    uint16_t buffer[32];
    __m512i counter[16]; __m512i counterU[16];
    const __m512i one  = _mm512_set1_epi16(1); // (11...1) vector
    const __m512i zero = _mm512_set1_epi16(0); // (00...0) vector

    // Mask selectors.
    const __m512i mask1 = _mm512_set1_epi16(256 + 2048); // FLAGSTAT_FSECONDARY (256) + FLAGSTAT_FSUPPLEMENTARY (2048)
    const __m512i mask2 = _mm512_set1_epi16(4 + 256 + 1024 + 2048); // Need to keep FLAGSTAT_FUNMAP (4) + FLAGSTAT_FDUP (1024) 
        // for Always rule and FLAGSTAT_FSECONDARY (256) and FLAGSTAT_FSUPPLEMENTARY (2048) for Rule 1 and Rule 2.
    const __m512i mask3 = _mm512_set1_epi16(512); // FLAGSTAT_FQCFAIL (512)

    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_setzero_si512();
            counterU[i] = _mm512_setzero_si512();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

        ///////////////////////////////////////////////////////////////////////
        // We load a register of data (data + i + j) and then construct the
        // conditional bits: 
        // 12: FLAGSTAT_FPROPER_PAIR + FLAGSTAT_FUNMAP == FLAGSTAT_FPROPER_PAIR
        // 13: FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP == FLAGSTAT_FMUNMAP
        // 14: FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP == 0
        //
        // These construction of these bits can be described for data x as:
        // x |= (x & LEFT_MASK == RIGHT_MASK) & 1 << TARGET_BIT
        // with the assumption that predicate evaluatons result in the selection
        // masks (00...0) or (11...1) for FALSE and TRUE, respectively. These
        // construction macros are named O1, O2, and O3.
        //
        // The original SAMtools method is also heavily branched with three
        // main branch points:
        // If FLAGSTAT_FSECONDARY then count FLAGSTAT_FSECONDARY
        // If FLAGSTAT_FSUPPLEMENTARY then count FLAGSTAT_FSUPPLEMENTARY
        // Else then count FLAGSTAT_FREAD1, 
        //                 FLAGSTAT_FREAD2,
        //                 Special bit 12, 13, and 14
        // Always count FLAGSTAT_FUNMAP, 
        //              FLAGSTAT_FDUP, 
        //              FLAGSTAT_FQCFAIL
        //
        // These bits can be selected using a mask-select propagate-carry approach:
        // x &= x & ((x == MASK) | CARRY_BITS)
        // with the arguments for MASK and CARRY_BITS as follows:
        //    1. {FLAGSTAT_FSECONDARY, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //    2. {FLAGSTAT_FSUPPLEMENTARY, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY 
        //        + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //    3. {FLAGSTAT_FPAIRED, 
        //        FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY 
        //        + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP}
        //
        // FLAGSTATS outputs summary statistics separately for reads that pass
        // QC and those that do not. Therefore we need to partition the data
        // into these two classes. For data that pass QC, the L registers, we
        // first bit-select the target FLAGSTAT_FQCFAIL bit using the mask
        // mask3. The resulting data is used to perform another mask-select
        // using VPCMPEQW against the empty vector (00...0). As above, if the
        // data has the FLAGSTAT_FQCFAIL bit set then this register will be
        // zeroed out. The exact process is performed for reads that fail QC,
        // the LU registers, with the difference that mask-selection is based on
        // the one vector (00...1).

#define LOAD(j) __m512i data##j = _mm512_loadu_si512(data + i + j) & (_mm512_cmpeq_epi16_mask( _mm512_loadu_si512(data + i + j) & mask1, zero ) | mask2);
#define L(j)  data##j & _mm512_cmpeq_epi16_mask( data##j & mask3, zero )
#define LU(j) data##j & _mm512_cmpeq_epi16_mask( data##j & mask3, mask3 )



        for (/**/; i < thislimit; i += 16) {
#define U(pos) {                     \
    counter[pos] = _mm512_add_epi16(counter[pos], _mm512_and_si512(v16, one)); \
    v16 = _mm512_srli_epi16(v16, 1); \
}
#define UU(pos) {                      \
    counterU[pos] = _mm512_add_epi16(counterU[pos], _mm512_and_si512(v16U, one)); \
    v16U = _mm512_srli_epi16(v16U, 1); \
}
            LOAD(0) LOAD(1)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 0),  L( 1));
            STORM_pospopcnt_csa_avx512(&twosAU,  &v1U, LU( 0), LU( 1));
            LOAD(2) LOAD(3)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L( 2),  L( 3));
            STORM_pospopcnt_csa_avx512(&twosBU,  &v1U, LU( 2), LU( 3));
            STORM_pospopcnt_csa_avx512(&foursA,  &v2,  twosA, twosB);
            STORM_pospopcnt_csa_avx512(&foursAU, &v2U, twosAU, twosBU);
            LOAD(4) LOAD(5)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 4),  L( 5));
            STORM_pospopcnt_csa_avx512(&twosAU,  &v1U, LU( 4), LU( 5));
            LOAD(6) LOAD(7)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L( 6),  L( 7));
            STORM_pospopcnt_csa_avx512(&twosBU,  &v1U, LU( 6), LU( 7));
            STORM_pospopcnt_csa_avx512(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx512(&foursBU, &v2U, twosAU,  twosBU);
            STORM_pospopcnt_csa_avx512(&eightsA, &v4,  foursA,  foursB);
            STORM_pospopcnt_csa_avx512(&eightsAU,&v4U, foursAU, foursBU);
            LOAD(8) LOAD(9)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 8),   L( 9));
            STORM_pospopcnt_csa_avx512(&twosAU,  &v1U, LU( 8),  LU( 9));
            LOAD(10) LOAD(11)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L(10),   L(11));
            STORM_pospopcnt_csa_avx512(&twosBU,  &v1U, LU(10),  LU(11));
            STORM_pospopcnt_csa_avx512(&foursA,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx512(&foursAU, &v2U, twosAU,  twosBU);
            LOAD(12) LOAD(13)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L(12),   L(13));
            STORM_pospopcnt_csa_avx512(&twosAU,  &v1U, LU(12),  LU(13));
            LOAD(14) LOAD(15)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L(14),   L(15));
            STORM_pospopcnt_csa_avx512(&twosBU,  &v1U, LU(14),  LU(15));
            STORM_pospopcnt_csa_avx512(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx512(&foursBU, &v2U, twosAU,  twosBU);
            STORM_pospopcnt_csa_avx512(&eightsB, &v4,  foursA,  foursB);
            STORM_pospopcnt_csa_avx512(&eightsBU,&v4U, foursAU, foursBU);
             U(0)  U(1)  U(2)  U(3)  U(4)  U(5)  U(6)  U(7)  U(8)  U(9)  U(10)  U(11)  U(12)  U(13)  U(14)  U(15) // Updates
            UU(0) UU(1) UU(2) UU(3) UU(4) UU(5) UU(6) UU(7) UU(8) UU(9) UU(10) UU(11) UU(12) UU(13) UU(14) UU(15) // Updates
            STORM_pospopcnt_csa_avx512(&v16,     &v8,  eightsA,  eightsB);
            STORM_pospopcnt_csa_avx512(&v16U,    &v8U, eightsAU, eightsBU);
#undef U
#undef UU
#undef LOAD
#undef L
#undef LU
        }
        
        // Update the counters after the last iteration
        for (size_t i = 0; i < 32; ++i) {
            counter[i]  = _mm512_add_epi16(counter[i], _mm512_and_si512(v16, one));
            v16  = _mm512_srli_epi16(v16, 1);
            counterU[i] = _mm512_add_epi16(counterU[i], _mm512_and_si512(v16U, one));
            v16U = _mm512_srli_epi16(v16U, 1);
        }
        
        for (size_t i = 0; i < 32; ++i) {
            _mm512_storeu_si512((__m512i*)buffer, counter[i]);
            for (size_t z = 0; z < 16; z++) {
                out[i] += 16 * (uint32_t)buffer[z];
            }

            _mm512_storeu_si512((__m512i*)buffer, counterU[i]);
            for (size_t z = 0; z < 16; z++) {
                out[16+i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v1);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v1U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[16+j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v2);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v2U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[16+j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v4);
    for (size_t i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v4U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[16+j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v8);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v8U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            out[16+j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    return 0;
}
#endif // end AVX512

static
uint64_t FLAGSTATS_u16(const uint16_t* array, uint32_t n_len, uint32_t* flags)
{

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

#if defined(STORM_HAVE_AVX512)
    if ((cpuid & STORM_CPUID_runtime_bit_AVX512BW) && n_len >= 1024) { // 16*512
        return FLAGSTAT_avx512(array, n_len, flags);
    }
#endif

#if defined(STORM_HAVE_AVX2)
    if ((cpuid & STORM_CPUID_runtime_bit_AVX2) && n_len >= 512) { // 16*256
        return FLAGSTAT_avx2(array, n_len, flags);
    }
    
    if ((cpuid & STORM_CPUID_runtime_bit_SSE42) && n_len >= 256) { // 16*128
        return FLAGSTAT_sse4(array, n_len, flags);
    }
#endif

#if defined(STORM_HAVE_SSE42)
    if ((cpuid & STORM_CPUID_runtime_bit_SSE42) && n_len >= 256) { // 16*128
        return FLAGSTAT_sse4(array, n_len, flags);
    }
#endif

    return FLAGSTAT_scalar(array, n_len, flags);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif