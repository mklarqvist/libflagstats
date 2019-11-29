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
/*! @abstract auxilary bit #12 set by SIMD procedures */
#define FLAGSTAT_BIT12              (1 << 12)
#define FLAGSTAT_BIT12_OFF          12
/*! @abstract auxilary bit #13 set by SIMD procedures */
#define FLAGSTAT_BIT13              (1 << 13)
#define FLAGSTAT_BIT13_OFF          13
/*! @abstract auxilary bit #14 set by SIMD procedures */
#define FLAGSTAT_BIT14              (1 << 14)
#define FLAGSTAT_BIT14_OFF          14

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
        if ( (val & FLAGSTAT_FPROPER_PAIR) && !(val & FLAGSTAT_FUNMAP) ) ++f[FLAGSTAT_BIT12_OFF];
        if (val & FLAGSTAT_FREAD1) ++f[FLAGSTAT_FREAD1_OFF];
        if (val & FLAGSTAT_FREAD2) ++f[FLAGSTAT_FREAD2_OFF];
        if ((val & FLAGSTAT_FMUNMAP) && !(val & FLAGSTAT_FUNMAP))  ++f[FLAGSTAT_BIT13_OFF];
        if (!(val & FLAGSTAT_FUNMAP) && !(val & FLAGSTAT_FMUNMAP)) ++f[FLAGSTAT_BIT14_OFF];
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
    return 0;
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

STORM_TARGET("sse4.2")
static
int FLAGSTAT_sse4_improved(const uint16_t* array, uint32_t len, uint32_t* flags) {
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

    const __m128i complete_bits_lookup = _mm_setr_epi8( // generated by expand_data.py
        0x00, 0x40, 0x00, 0x50, 0x00, 0x00, 0x00, 0x00, 0x00, 0x20, 0x00, 0x30, 0x00, 0x00, 0x00, 0x00);

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
        /*
            The bits span the lower 4 bits of flag set and this makes that
            suitable for use pshufb to quickly lookup desired bits. The lookup
            must be done only for the higer byte of 16-bit word, the lower
            byte must be zeroed. Such a transformed input must be or-ed with
            the input word.
                               FLAGSTAT_FMUNMAP
                               |FLAGSTAT_FUNMAP
                               ||FLAGSTAT_FPROPER_PAIR
                               |||FLAGSTAT_FPAIRED
                               ||||
            input word:  [0000|xxxx|1000|0000]
                                    |
                                    forces zeoring
                           bit14
                           |bit13
                           ||bit12
                           |||
            output word: [0abc|0000|0000|0000]
        */
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
#define O1(j) const __m128i complete_index##j = \
                (_mm_slli_epi16(data##j, 8) & _mm_set1_epi16(0x0f00)) | _mm_set1_epi16(0x0080);
#define O2(j) data##j = data##j | _mm_shuffle_epi8(complete_bits_lookup, complete_index##j);
#define O3(j) 
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
int FLAGSTAT_avx512(const uint16_t* array, uint32_t len, uint32_t* flags) {
    const uint32_t start_qc = flags[FLAGSTAT_FQCFAIL_OFF + 16];
    
    for (uint32_t i = len - (len % (32 * 16)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], flags);
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
    __m512i counter[16]; 
    __m512i counterU[16];
    
    // Masks and mask selectors.
    const __m512i m1   = _mm512_set1_epi16(FLAGSTAT_FSECONDARY);
    const __m512i m1S  = _mm512_set1_epi16(FLAGSTAT_FQCFAIL + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP);
    const __m512i m2   = _mm512_set1_epi16(FLAGSTAT_FSUPPLEMENTARY);
    const __m512i m2S  = _mm512_set1_epi16(FLAGSTAT_FQCFAIL + FLAGSTAT_FSUPPLEMENTARY + FLAGSTAT_FSECONDARY + FLAGSTAT_FUNMAP + FLAGSTAT_FDUP);
    const __m512i m3   = _mm512_set1_epi16(FLAGSTAT_FPAIRED);
    const __m512i m4   = _mm512_set1_epi16(FLAGSTAT_FQCFAIL);
    const __m512i one  = _mm512_set1_epi16(1); // (00...1) vector
    const __m512i zero = _mm512_set1_epi16(0); // (00...0) vector

    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_setzero_si512();
            counterU[i] = _mm512_setzero_si512();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

#define W(j) __m512i data##j = _mm512_loadu_si512(data + i + j);
#define O1(j) data##j = data##j |   _mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask(data##j & _mm512_set1_epi16(FLAGSTAT_FPROPER_PAIR + FLAGSTAT_FUNMAP), _mm512_set1_epi16(FLAGSTAT_FPROPER_PAIR)), (uint16_t)1 << 12); 
#define O2(j) data##j = data##j |   _mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask(data##j & _mm512_set1_epi16(FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP), _mm512_set1_epi16(FLAGSTAT_FMUNMAP)), (uint16_t)1 << 13);
#define O3(j) data##j = data##j |   _mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask(data##j & _mm512_set1_epi16(FLAGSTAT_FMUNMAP + FLAGSTAT_FUNMAP), zero), (uint16_t)1 << 14);
#define L1(j) data##j = data##j & (_mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask((data##j & m1), zero),65535) | m1S);
#define L2(j) data##j = data##j & (_mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask((data##j & m2), zero),65535) | m2S);
#define L3(j) data##j = data##j & (_mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask((data##j & m3), m3),  65535) | m2S);
#define LOAD(j) W(j) O1(j) O2(j) O3(j) L1(j) L2(j) L3(j)
#define L(j)  data##j & _mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask( data##j & m4, zero ), 65535)
#define LU(j) data##j & _mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask( data##j & m4, m4 ),   65535)

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
            counter[i]  = _mm512_add_epi16(counter[i], _mm512_and_si512(v16, one));
            v16  = _mm512_srli_epi16(v16, 1);
            counterU[i] = _mm512_add_epi16(counterU[i], _mm512_and_si512(v16U, one));
            v16U = _mm512_srli_epi16(v16U, 1);
        }
        
        for (size_t i = 0; i < 16; ++i) {
            _mm512_storeu_si512((__m512i*)buffer, counter[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[i] += 16 * (uint32_t)buffer[z];
            }

            _mm512_storeu_si512((__m512i*)buffer, counterU[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[16+i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v1);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v1U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v2);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v2U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v4);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v4U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v8);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v8U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    // QC
    flags[FLAGSTAT_FQCFAIL_OFF] += len - (flags[FLAGSTAT_FQCFAIL_OFF+16] - start_qc);

    return 0;
}

STORM_TARGET("avx512bw")
static
int FLAGSTAT_avx512_improved(const uint16_t* array, uint32_t len, uint32_t* flags) {
    const uint32_t start_qc = flags[FLAGSTAT_FQCFAIL_OFF + 16];
    
    for (uint32_t i = len - (len % (32 * 16)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], flags);
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
    __m512i counter[16]; 
    __m512i counterU[16];
    
    // Masks and mask selectors.
    const __m512i one  = _mm512_set1_epi16(1); // (00...1) vector

    // generated by expand_data.py
    const __m512i complete_bits_lookup = _mm512_setr_epi32(
        0x40000000, 0x50000000, 0x00000000, 0x00000000, 0x20000000, 0x30000000, 0x00000000, 0x00000000,
        0x40000000, 0x50000000, 0x00000000, 0x00000000, 0x20000000, 0x30000000, 0x00000000, 0x00000000);

    // generated by mask_data.py
    const __m512i qcfail_1_lookup = _mm512_setr_epi32(
        0x00000000, 0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04,
        0x00000000, 0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04);

    // generated by mask_data.py
    const __m512i qcfail_0_lookup = _mm512_setr_epi32(
        0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04, 0x00000000,
        0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04, 0x00000000);

    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_setzero_si512();
            counterU[i] = _mm512_setzero_si512();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;


#define W(j) __m512i data##j = _mm512_loadu_si512(data + i + j);
#define O1(j) data##j = data##j | _mm512_permutexvar_epi16(data##j, complete_bits_lookup);
        /*
            We're gathering bits that decide about the major control flow.
            The resulting value will be issued to lookup instructions, that
            provide masks for data.

                  FLAGSTAT_FSUPPLEMENTARY
                  | FLAGSTAT_FQCFAIL
                  | |FLAGSTAT_FSECONDARY
                  | ||         FLAGSTAT_FPAIRED
                  | ||         |
            [....|x.yw|....|...z]  (. = garbage)

            This is the layout of words obtained in L1

            [0000|0000|0000|xzyw]
        */

// We're using 32-bit shifts, as they're faster than 16-bit shifts.
#define L1(j) const __m512i mask_index##j = \
                (_mm512_srli_epi32(data##j, 8) & _mm512_set1_epi16(0x0b)) | \
                (_mm512_slli_epi32(data##j, 2) & _mm512_set1_epi16(0x04));
#define L2(j) const __m512i mask_qcfail0_##j = _mm512_permutexvar_epi16(mask_index##j, qcfail_0_lookup);
#define L3(j) const __m512i mask_qcfail1_##j = _mm512_permutexvar_epi16(mask_index##j, qcfail_1_lookup);
#define LOAD(j) W(j) O1(j) L1(j) L2(j) L3(j)
#define L(j)  _mm512_and_si512(data##j, mask_qcfail0_##j)
#define LU(j) _mm512_and_si512(data##j, mask_qcfail1_##j)

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
#undef W
#undef O1
#undef L1
#undef L2
#undef L3
        }

        // Update the counters after the last iteration
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_add_epi16(counter[i], _mm512_and_si512(v16, one));
            v16  = _mm512_srli_epi16(v16, 1);
            counterU[i] = _mm512_add_epi16(counterU[i], _mm512_and_si512(v16U, one));
            v16U = _mm512_srli_epi16(v16U, 1);
        }
        
        for (size_t i = 0; i < 16; ++i) {
            _mm512_storeu_si512((__m512i*)buffer, counter[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[i] += 16 * (uint32_t)buffer[z];
            }

            _mm512_storeu_si512((__m512i*)buffer, counterU[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[16+i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v1);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v1U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v2);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v2U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v4);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v4U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v8);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v8U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    // QC
    flags[FLAGSTAT_FQCFAIL_OFF] += len - (flags[FLAGSTAT_FQCFAIL_OFF+16] - start_qc);

    return 0;
}

STORM_TARGET("avx512bw")
static
int FLAGSTAT_avx512_improved2(const uint16_t* array, uint32_t len, uint32_t* flags) {
    const uint32_t start_qc = flags[FLAGSTAT_FQCFAIL_OFF + 16];
    
    for (uint32_t i = len - (len % (32 * 16)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], flags);
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
    __m512i counter[16]; 
    __m512i counterU[16];
    
    // Masks and mask selectors.
    const __m512i one  = _mm512_set1_epi16(1); // (00...1) vector

    /*
                     FLAGSTAT_FSUPPLEMENTARY
                     | FLAGSTAT_FQCFAIL
                     | |FLAGSTAT_FSECONDARY
                     | ||      FLAGSTAT_FMUNMAP
                     | ||      |FLAGSTAT_FUNMAP
                     | ||      ||FLAGSTAT_FPROPER_PAIR
                     | ||      |||FLAGSTAT_FPAIRED
                     | ||      ||||
               [....|e.fg|....|abcd]  (. = garbage)

            From the marked bits we need two extra masks:

            F1 [0ABC|0000|0000|0000] (A, B, C are calculated from a, b, c, d)
            F2 [0000|0000|0000|edfg]

            The first mask is produced by VPSHUFW. But at this step we can
            also place bit d at proper position in the 3rd nibble:

            F3 [0ABC|0d00|0000|0000]

            The mask F1 must be or'ed with the input data, but we can use
            F3 and ternary logic operator to perform merge:

            // 0xca defines function (m & F3) | (~m & data)
            data = _mm512_ternarylogic_epi32(m=[0111|0000|0000|0000], F3=[0ABC|0d00|0000|0000], data, 0xca)

            The mask F2 can be also obtained from data vector and F3

            // F2 = [....|edfg|....|....]
            F2 = _mm512_ternarylogic_epi32(m=[0000|0100|0000|0000], F3=[0ABC|0d00|0000|0000], data, 0xca)
            // F2 = [0000|0000|....|edfg] 
            F2 = _mm512_srli_epi32(tmp, 8);
    */

    // generated by expand_data.py
    const __m512i complete_bits_lookup = _mm512_setr_epi32(
        0x44000000, 0x54000000, 0x04000000, 0x04000000, 0x24000000, 0x34000000, 0x04000000, 0x04000000,
        0x44000000, 0x54000000, 0x04000000, 0x04000000, 0x24000000, 0x34000000, 0x04000000, 0x04000000);

    // generated by mask_data.py
    const __m512i qcfail_1_lookup = _mm512_setr_epi32(
        0x00000000, 0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04,
        0x00000000, 0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04);

    // generated by mask_data.py
    const __m512i qcfail_0_lookup = _mm512_setr_epi32(
        0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04, 0x00000000,
        0x07040604, 0x00000000, 0x070476c4, 0x00000000, 0x07040e04, 0x00000000, 0x07040e04, 0x00000000);

    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_setzero_si512();
            counterU[i] = _mm512_setzero_si512();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

#define W(j) __m512i data##j = _mm512_loadu_si512(data + i + j);
#define O1(j) const __m512i mask##j = _mm512_permutexvar_epi16(data##j, complete_bits_lookup);
#define O2(j) data##j = _mm512_ternarylogic_epi32(_mm512_set1_epi16(0x7000), mask##j, data##j, 0xca);
#define L1(j) const __m512i mask_index##j = _mm512_srli_epi16(_mm512_ternarylogic_epi32(_mm512_set1_epi16(0x0400), mask##j, data##j, 0xca), 8);
#define L2(j) const __m512i mask_qcfail0_##j = _mm512_permutexvar_epi16(mask_index##j, qcfail_0_lookup);
#define L3(j) const __m512i mask_qcfail1_##j = _mm512_permutexvar_epi16(mask_index##j, qcfail_1_lookup);
#define LOAD(j) W(j) O1(j) O2(j) L1(j) L2(j) L3(j)
#define L(j)  _mm512_and_si512(data##j, mask_qcfail0_##j)
#define LU(j) _mm512_and_si512(data##j, mask_qcfail1_##j)

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
#undef W
#undef O1
#undef L1
#undef L2
#undef L3
        }

        // Update the counters after the last iteration
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_add_epi16(counter[i], _mm512_and_si512(v16, one));
            v16  = _mm512_srli_epi16(v16, 1);
            counterU[i] = _mm512_add_epi16(counterU[i], _mm512_and_si512(v16U, one));
            v16U = _mm512_srli_epi16(v16U, 1);
        }
        
        for (size_t i = 0; i < 16; ++i) {
            _mm512_storeu_si512((__m512i*)buffer, counter[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[i] += 16 * (uint32_t)buffer[z];
            }

            _mm512_storeu_si512((__m512i*)buffer, counterU[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[16+i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v1);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v1U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v2);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v2U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v4);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v4U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v8);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }
    _mm512_storeu_si512((__m512i*)buffer, v8U);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[16+j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    // QC
    flags[FLAGSTAT_FQCFAIL_OFF] += len - (flags[FLAGSTAT_FQCFAIL_OFF+16] - start_qc);

    return 0;
}

/*
    Variant 3
    --------------------------------------------------------------

    This variant is based on the observation that we're doing two
    full 16-bit positional popcounts (32 counters), while in fact
    we have to count only 19 counters. The FLAGSTAT_FQCFAIL is counted
    unconditionally and then following 8 bits are counted twice for
    FLAGSTAT_FQCFAIL = 0 or 1: FLAGSTAT_FUNMAP, FLAGSTAT_FDUP,
    FLAGSTAT_FSUPPLEMENTARY, FLAGSTAT_FSECONDARY, FLAGSTAT_FREAD1,
    FLAGSTAT_FREAD2, FLAGSTAT_BIT12, FLAGSTAT_BIT13, FLAGSTAT_BIT14.

    The idea is to duplicate most of these bits in a 16-bit word
    and then mask them depending on FQCFAIL. The 16-bit word is
    then issued to pospopcnt procedure. The remaining three bits
    are counted separately.

    The bits we want to re-shuffle (it's layout from the previous version,
    so bits #12, 13 and 14 are not actually occupied in the input word):

             BIT#14
             |   BIT#13
             |   |   BIT#12
             |   |   |   FSUPPLEMENTARY
             |   |   |   |   FDUP
             |   |   |   |   |   FQCFAIL
             |   |   |   |   |   |   FSECONDARY
             |   |   |   |   |   |   |   FREAD2
             |   |   |   |   |   |   |   |   FREAD1          FUNMAP  FPAIRED
             |   |   |   |   |   |   |   |   |               |       |
       [ . | a | b | c | d | e | f | g | h | i | . | . | . | j | . | p ] = input
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    1. In the first step we need to calculate bits #12, #13 and #14 using
       lookup as in the previous approaches (bits 0..3 are indices for lookup).
       (notation: a0 - bit for case when FQCFAIL=0, a1 - FQCFAIL=1):

       [ . | . | . | p | . | . | j | . | . | . | a1| b1| c1| a0| b0| c0] = t0
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    2. In the second step we need to duplicate bits FREAD1 (i) and FREAD2 (h):

       [ . | . | . | . | . | . | h1| i1| h0| i0| . | . | . | . | . | . ] = t1
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    3. We marge vector 'input' with 't0' and shift it right to form an
       index for following lookups:

       [ . | . | . | p | d | e | j | g | . | . | . | . | . | . | . | . ] --- merged
       [ . | . | . | . | . | . | . | . | . | . | . | p | d | e | j | g ] = idx
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    4. The index 'idx' is now used to build two vectors with:
       a. duplicated FSUPPLEMENTARY (d), FDUP (e) and FSECONDARY (g), and also
          a copy of FUNMAP (j)

       [ d1| d0| e1| e0| g1| g0| . | . | . | . | . | . | . | . | . | j ] = t2
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

       b. build masks for condition depending on flags: FSUPPLEMENTARY, FSECONDARY
          and FPAIRED: let's call this 'cond_mask'

    5. Build vector with flags depending on 'cond_mask' from intermediate vectors
       t0, t1 and t2.

       [ d1| d0| e1| e0| g1| g0| h1| i1| h0| i0| . | . | . | . | . | j ] = t1 | t2
       [ d1| d0| e1| e0| g1| g0| h1| i1| h0| i0| a1| a0| b1| b0| c1| c0] = t3 = bitcond(0x003f, t0, (t1 | t2))
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    6. Populate FQCFAIL (f) [note: as decimal 0 or -1]:

       [ f | f | f | f | f | f | f | f | f | f | f | f | f | f | f | f ] = qcfail

    7. Mask out bits from t3 depending on FQCFAIL -- we get mask 'qcfail_mask',
       equals either:

       [ 0 | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 1 | 1 | 0 | 1 | 0 | 1 | 0 | 1 ] --- qcfail=0
       [ 1 | 0 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 ] --- qcfail=1
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

       qcfail_mask = (qcfail & ~0x54d5) | (~qcfail & 0x54d5) = qcfail ^ 0x54d5

    8. Call pospopcnt(t3 & qcfail_mask)

    7. Separately count FUNMAP (j), which is now present in t2:
       - increment local 16-bit counter: t2 & qcfail & vector(0x0001)   # ternary func: 0x80
       - increment local 16-bit counter: t2 & ~qcfail & vector(0x0001)  # ternary func: 0x20

    8. add to local 16-bit counter: qcfail

    9. After pospopcnt, update the global counters with these local ones updated in #7 and #8.
*/
#define AVX512_BIT12_FQCFAIL_0          0
#define AVX512_BIT12_FQCFAIL_1          1
#define AVX512_BIT13_FQCFAIL_0          2
#define AVX512_BIT13_FQCFAIL_1          3
#define AVX512_BIT14_FQCFAIL_0          4
#define AVX512_BIT14_FQCFAIL_1          5
#define AVX512_FREAD1_FQCFAIL_0         6
#define AVX512_FREAD2_FQCFAIL_0         7
#define AVX512_FREAD1_FQCFAIL_1         8
#define AVX512_FREAD2_FQCFAIL_1         9
#define AVX512_FSECONDARY_FQCFAIL_0     10
#define AVX512_FSECONDARY_FQCFAIL_1     11
#define AVX512_FDUP_FQCFAIL_0           12
#define AVX512_FDUP_FQCFAIL_1           13
#define AVX512_FSUPPLEMENTARY_FQCFAIL_0 14
#define AVX512_FSUPPLEMENTARY_FQCFAIL_1 15

#define avx512_bitcond(condition, true_val, false_val) _mm512_ternarylogic_epi32((condition), (true_val), (false_val), 0xca)

STORM_TARGET("avx512bw")
static
uint32_t avx512_sum_epu32(__m512i v) {
    uint32_t tmp[16];
    uint32_t sum = 0;

    _mm512_storeu_si512(tmp, v);
    for (int i=0; i < 16; i++)
        sum += tmp[i];

    return sum;
}

STORM_TARGET("avx512bw")
static
int FLAGSTAT_avx512_improved3(const uint16_t* array, uint32_t len, uint32_t* user_flags) {
    for (uint32_t i = len - (len % (32 * 16)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], user_flags);
    }

    uint32_t flags[16];
    for (int i=0; i < 16; i++)
        flags[i] = 0;

    const __m512i* data = (const __m512i*)array;
    size_t size = len / 32;
    __m512i v1  = _mm512_setzero_si512();
    __m512i v2  = _mm512_setzero_si512();
    __m512i v4  = _mm512_setzero_si512();
    __m512i v8  = _mm512_setzero_si512();
    __m512i v16 = _mm512_setzero_si512();
    __m512i twosA, twosB, foursA, foursB, eightsA, eightsB;

    const uint64_t limit = size - size % 16;
    uint64_t i = 0;
    uint16_t buffer[32];
    __m512i counter[16];

    // Masks and mask selectors.
    const __m512i one  = _mm512_set1_epi16(1); // (00...1) vector

    // generated by scripts/expand_data3.py
    const __m512i complete_bits_lookup = _mm512_setr_epi32(
        0x10300000, 0x10330000, 0x12000200, 0x12000200, 0x100c0000, 0x100f0000, 0x12000200, 0x12000200,
        0x10300000, 0x10330000, 0x12000200, 0x12000200, 0x100c0000, 0x100f0000, 0x12000200, 0x12000200
    );

    // generated by scripts/mask_data3.py
    const __m512i duplicate_bits_lookup = _mm512_setr_epi32(
        0x0c000000, 0x0c010001, 0x3c003000, 0x3c013001, 0xcc00c000, 0xcc01c001, 0xfc00f000, 0xfc01f001,
        0x0c000000, 0x0c010001, 0x3c003000, 0x3c013001, 0xcc00c000, 0xcc01c001, 0xfc00f000, 0xfc01f001
    );

    // generated by scripts/mask_data3.py
    const __m512i condition_mask_lookup = _mm512_setr_epi32(
        0x3c003000, 0x3c003000, 0x3c003000, 0x3c003000, 0x3c00f000, 0x3c00f000, 0x3c00f000, 0x3c00f000,
        0x3c0033ff, 0x3c0033ff, 0x3c0033ff, 0x3c0033ff, 0x3c00f000, 0x3c00f000, 0x3c00f000, 0x3c00f000
    );

    // 32-bit counters
    __m512i qcfail_global_counter = _mm512_setzero_si512();
    __m512i funmap_global_counter_qcfail_0 = _mm512_setzero_si512();
    __m512i funmap_global_counter_qcfail_1 = _mm512_setzero_si512();

    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_setzero_si512();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

#define LOAD(j) \
        __m512i data##j              = _mm512_loadu_si512(data + i + j); \
        const __m512i t0_##j         = _mm512_permutexvar_epi16(data##j, complete_bits_lookup); \
        const __m512i t1a_##j        = _mm512_and_si512(data##j, _mm512_set1_epi16(0x00c0)); \
        const __m512i t1_##j         = _mm512_or_si512(t1a_##j, _mm512_slli_epi32(t1a_##j, 2)); \
        const __m512i idx##j         = _mm512_srli_epi32(avx512_bitcond(_mm512_set1_epi16(0x1200), t0_##j, data##j), 8); \
        const __m512i t2_##j         = _mm512_permutexvar_epi16(idx##j, duplicate_bits_lookup); \
        const __m512i cond_mask##j   = _mm512_permutexvar_epi16(idx##j, condition_mask_lookup); \
        const __m512i t3_##j         = avx512_bitcond(_mm512_set1_epi16(0x003f), t0_##j, (t1_##j | t2_##j)); \
        const __m512i qcfail##j      = _mm512_srai_epi16(_mm512_slli_epi32(data##j, 6), 16); \
        const __m512i qcfail_mask##j = _mm512_xor_si512(qcfail##j, _mm512_set1_epi16(0x54d5)); \
        qcfail_counter = _mm512_sub_epi16(qcfail_counter, qcfail##j); \
        funmap_counter_qcfail_0 = _mm512_add_epi16(funmap_counter_qcfail_0, \
                                          _mm512_ternarylogic_epi32(t2_##j, qcfail##j, one, 0x20)); \
        funmap_counter_qcfail_1 = _mm512_add_epi16(funmap_counter_qcfail_1, \
                                          _mm512_ternarylogic_epi32(t2_##j, qcfail##j, one, 0x80));

#define L(j) (t3_##j & cond_mask##j & qcfail_mask##j)

        for (/**/; i < thislimit; i += 16) {

            // 16-bit counters
            __m512i qcfail_counter = _mm512_setzero_si512();
            __m512i funmap_counter_qcfail_0 = _mm512_setzero_si512();
            __m512i funmap_counter_qcfail_1 = _mm512_setzero_si512();

#define U(pos) {                     \
    counter[pos] = _mm512_add_epi16(counter[pos], _mm512_and_si512(v16, one)); \
    v16 = _mm512_srli_epi16(v16, 1); \
}
            LOAD(0) LOAD(1)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 0),  L( 1));
            LOAD(2) LOAD(3)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L( 2),  L( 3));
            STORM_pospopcnt_csa_avx512(&foursA,  &v2,  twosA, twosB);
            LOAD(4) LOAD(5)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 4),  L( 5));
            LOAD(6) LOAD(7)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L( 6),  L( 7));
            STORM_pospopcnt_csa_avx512(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx512(&eightsA, &v4,  foursA,  foursB);
            LOAD(8) LOAD(9)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 8),   L( 9));
            LOAD(10) LOAD(11)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L(10),   L(11));
            STORM_pospopcnt_csa_avx512(&foursA,  &v2,  twosA,   twosB);
            LOAD(12) LOAD(13)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L(12),   L(13));
            LOAD(14) LOAD(15)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L(14),   L(15));
            STORM_pospopcnt_csa_avx512(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx512(&eightsB, &v4,  foursA,  foursB);
             U(0)  U(1)  U(2)  U(3)  U(4)  U(5)  U(6)  U(7)  U(8)  U(9)  U(10)  U(11)  U(12)  U(13)  U(14)  U(15) // Updates
            STORM_pospopcnt_csa_avx512(&v16,     &v8,  eightsA,  eightsB);

            qcfail_global_counter = _mm512_add_epi32(qcfail_global_counter,
                                                     _mm512_madd_epi16(qcfail_counter, one));
            funmap_global_counter_qcfail_0 =
                _mm512_add_epi32(funmap_global_counter_qcfail_0,
                _mm512_madd_epi16(funmap_counter_qcfail_0, one));
            funmap_global_counter_qcfail_1 =
                _mm512_add_epi32(funmap_global_counter_qcfail_1,
                _mm512_madd_epi16(funmap_counter_qcfail_1, one));
#undef U
#undef UU
#undef LOAD
#undef L
#undef LU
#undef W
#undef O1
#undef L1
#undef L2
#undef L3
        }

        // Update the counters after the last iteration
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_add_epi16(counter[i], _mm512_and_si512(v16, one));
            v16  = _mm512_srli_epi16(v16, 1);
        }

        for (size_t i = 0; i < 16; ++i) {
            _mm512_storeu_si512((__m512i*)buffer, counter[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v1);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v2);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v4);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v8);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    // update flags from our custom flags table
    user_flags[FLAGSTAT_FQCFAIL_OFF + 16]        += avx512_sum_epu32(qcfail_global_counter);
    user_flags[FLAGSTAT_FUNMAP_OFF + 0]          += avx512_sum_epu32(funmap_global_counter_qcfail_0);
    user_flags[FLAGSTAT_FUNMAP_OFF + 16]         += avx512_sum_epu32(funmap_global_counter_qcfail_1);
    user_flags[FLAGSTAT_FDUP_OFF + 0]            += flags[AVX512_FDUP_FQCFAIL_0];
    user_flags[FLAGSTAT_FDUP_OFF + 16]           += flags[AVX512_FDUP_FQCFAIL_1];
    user_flags[FLAGSTAT_FSECONDARY_OFF + 0]      += flags[AVX512_FSECONDARY_FQCFAIL_0];
    user_flags[FLAGSTAT_FSECONDARY_OFF + 16]     += flags[AVX512_FSECONDARY_FQCFAIL_1];
    user_flags[FLAGSTAT_FSUPPLEMENTARY_OFF + 0]  += flags[AVX512_FSUPPLEMENTARY_FQCFAIL_0];
    user_flags[FLAGSTAT_FSUPPLEMENTARY_OFF + 16] += flags[AVX512_FSUPPLEMENTARY_FQCFAIL_1];
    user_flags[FLAGSTAT_BIT12_OFF + 0]           += flags[AVX512_BIT12_FQCFAIL_0];
    user_flags[FLAGSTAT_BIT12_OFF + 16]          += flags[AVX512_BIT12_FQCFAIL_1];
    user_flags[FLAGSTAT_BIT13_OFF + 0]           += flags[AVX512_BIT13_FQCFAIL_0];
    user_flags[FLAGSTAT_BIT13_OFF + 16]          += flags[AVX512_BIT13_FQCFAIL_1];
    user_flags[FLAGSTAT_BIT14_OFF + 0]           += flags[AVX512_BIT14_FQCFAIL_0];
    user_flags[FLAGSTAT_BIT14_OFF + 16]          += flags[AVX512_BIT14_FQCFAIL_1];
    user_flags[FLAGSTAT_FREAD1_OFF + 0]          += flags[AVX512_FREAD1_FQCFAIL_0];
    user_flags[FLAGSTAT_FREAD1_OFF + 16]         += flags[AVX512_FREAD1_FQCFAIL_1];
    user_flags[FLAGSTAT_FREAD2_OFF + 0]          += flags[AVX512_FREAD2_FQCFAIL_0];
    user_flags[FLAGSTAT_FREAD2_OFF + 16]         += flags[AVX512_FREAD2_FQCFAIL_1];

    return 0;
}
#undef AVX512_BIT12_FQCFAIL_0
#undef AVX512_BIT12_FQCFAIL_1
#undef AVX512_BIT13_FQCFAIL_0
#undef AVX512_BIT13_FQCFAIL_1
#undef AVX512_BIT14_FQCFAIL_0
#undef AVX512_BIT14_FQCFAIL_1
#undef AVX512_FREAD1_FQCFAIL_0
#undef AVX512_FREAD2_FQCFAIL_0
#undef AVX512_FREAD1_FQCFAIL_1
#undef AVX512_FREAD2_FQCFAIL_1
#undef AVX512_FSECONDARY_FQCFAIL_0
#undef AVX512_FSECONDARY_FQCFAIL_1
#undef AVX512_FDUP_FQCFAIL_0
#undef AVX512_FDUP_FQCFAIL_1
#undef AVX512_FSUPPLEMENTARY_FQCFAIL_0
#undef AVX512_FSUPPLEMENTARY_FQCFAIL_1

/*
    Variant 4
    --------------------------------------------------------------

    The idea is similar to #3, but we don't scatter the bits of interest
    over the whole 16-bit word, but put them in lower byte and then
    duplicate that byte.

    The bits we want to re-shuffle (it's layout from the previous versions,
    so bits #12, 13 and 14 are not actually occupied in the input word):

             BIT#14
             |   BIT#13
             |   |   BIT#12
             |   |   |   FSUPPLEMENTARY
             |   |   |   |   FDUP
             |   |   |   |   |   FQCFAIL
             |   |   |   |   |   |   FSECONDARY
             |   |   |   |   |   |   |   FREAD2
             |   |   |   |   |   |   |   |   FREAD1          FUNMAP  FPAIRED
             |   |   |   |   |   |   |   |   |               |       |
       [ . | a | b | c | d | e | f | g | h | i | . | . | . | j | . | p ] = input
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    1. In the first step we need to calculate bits #12, #13 and #14 using
       lookup as in the previous approaches (bits 0..3 are indices for lookup).

       [ . | . | . | p | . | . | . | . | . | . | a | b | c | . | . | . ] = t0
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    2. We marge vector 'input' with 't0' and shift it right to form an
       index for following lookups:

       [ . | . | . | p | d | e | f | g | . | . | . | . | . | . | . | . ] --- merged
       [ . | . | . | . | . | . | . | . | . | . | . | p | d | e | f | g ] = idx
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    3. The index 'idx' is now used to build two vectors with:
       a. shuffled FSUPPLEMENTARY (d), FDUP (e), FSECONDARY (g), and FQCFAIL (d).

       [ f | . | . | . | . | . | . | . | . | . | . | . | . | d | e | g ] = t1
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

       b. masks for condition depending on flags: FSUPPLEMENTARY, FSECONDARY
          and FPAIRED: let's call it 'cond_mask'

    5. Populate FQCFAIL (f) which is now MSB of t1 [note: as decimal 0 or -1]:

       [ f | f | f | f | f | f | f | f | f | f | f | f | f | f | f | f ] = qcfail
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    4. Merge t0 and t1. The lower byte has got almost all flags: FREAD1(h),
       FREAD2(i), BIT12(a), BIT13(b), BIT14(c), FSUPPLEMENTARY(d), FDUP(e),
       FSECONDARY(g). Only FUNMAP(j) and FQCFAIL(f) have to be counted separately.

       [ 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | h | i | a | b | c | d | e | g ] = t2
         15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   0

    5. If FQCFAIL=1 then shift t2 by 8.

    6. Call pospopcnt(t3 & qcfail)

    7. Separately count FUNMAP (j), which is bit #2 of input
       - increment local 16-bit counter: t2 & ~qcfail & vector(0x0004)
       - increment local 16-bit counter: t2 & qcfail & vector(0x0004)

    8. add to local 16-bit counter: qcfail

    9. After pospopcnt, update the global counters with these local ones updated in #7 and #8.
       FUNMAP local counter should be divided by 4, and FQCFAIL's higer byte must be negated.
*/
#define AVX512_BIT12_FQCFAIL_0          5
#define AVX512_BIT12_FQCFAIL_1          13
#define AVX512_BIT13_FQCFAIL_0          4
#define AVX512_BIT13_FQCFAIL_1          12
#define AVX512_BIT14_FQCFAIL_0          3
#define AVX512_BIT14_FQCFAIL_1          11
#define AVX512_FREAD1_FQCFAIL_0         7
#define AVX512_FREAD2_FQCFAIL_0         15
#define AVX512_FREAD1_FQCFAIL_1         6
#define AVX512_FREAD2_FQCFAIL_1         14
#define AVX512_FSECONDARY_FQCFAIL_0     0
#define AVX512_FSECONDARY_FQCFAIL_1     8
#define AVX512_FDUP_FQCFAIL_0           1
#define AVX512_FDUP_FQCFAIL_1           9
#define AVX512_FSUPPLEMENTARY_FQCFAIL_0 2
#define AVX512_FSUPPLEMENTARY_FQCFAIL_1 10

#define avx512_bitcond(condition, true_val, false_val) _mm512_ternarylogic_epi32((condition), (true_val), (false_val), 0xca)
#define epi16 _mm512_set1_epi16

STORM_TARGET("avx512bw")
static
int FLAGSTAT_avx512_improved4(const uint16_t* array, uint32_t len, uint32_t* user_flags) {
    for (uint32_t i = len - (len % (32 * 16)); i < len; ++i) {
        FLAGSTAT_scalar_update(array[i], user_flags);
    }

    uint32_t flags[16];
    for (int i=0; i < 16; i++)
        flags[i] = 0;

    const __m512i* data = (const __m512i*)array;
    size_t size = len / 32;
    __m512i v1  = _mm512_setzero_si512();
    __m512i v2  = _mm512_setzero_si512();
    __m512i v4  = _mm512_setzero_si512();
    __m512i v8  = _mm512_setzero_si512();
    __m512i v16 = _mm512_setzero_si512();
    __m512i twosA, twosB, foursA, foursB, eightsA, eightsB;

    const uint64_t limit = size - size % 16;
    uint64_t i = 0;
    uint16_t buffer[32];
    __m512i counter[16];

    // Masks and mask selectors.
    const __m512i one  = _mm512_set1_epi16(1); // (00...1) vector

    // generated by scripts/version5.py
    const __m512i complete_bits_lookup = _mm512_setr_epi32(
        0x18080000, 0x38280000, 0x10000000, 0x10000000, 0x10100000, 0x30300000, 0x10000000, 0x10000000,
        0x18080000, 0x38280000, 0x10000000, 0x10000000, 0x10100000, 0x30300000, 0x10000000, 0x10000000
    );

    const __m512i reshuffle_bits_lookup = _mm512_setr_epi32(
        0x00010000, 0x80018000, 0x00030002, 0x80038002, 0x00050004, 0x80058004, 0x00070006, 0x80078006,
        0x00010000, 0x80018000, 0x00030002, 0x80038002, 0x00050004, 0x80058004, 0x00070006, 0x80078006
    );

    const __m512i condition_mask_lookup = _mm512_setr_epi32(
        0x03030202, 0x03030202, 0x03030202, 0x03030202, 0x03030606, 0x03030606, 0x03030606, 0x03030606,
        0x0303fafa, 0x0303fafa, 0x0303fafa, 0x0303fafa, 0x03030606, 0x03030606, 0x03030606, 0x03030606
    );
    // end of autogenered content

    // 32-bit counters
    __m512i qcfail_global_counter = _mm512_setzero_si512();
    __m512i funmap_global_counter_qcfail_0 = _mm512_setzero_si512();
    __m512i funmap_global_counter_qcfail_1 = _mm512_setzero_si512();

    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_setzero_si512();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

#define LOAD(j) \
        __m512i data##j              = _mm512_loadu_si512(data + i + j); \
        const __m512i t0_##j         = _mm512_permutexvar_epi16(data##j, complete_bits_lookup); \
        const __m512i idx##j         = _mm512_srli_epi32(avx512_bitcond(epi16(0x1000), t0_##j, data##j), 8); \
        const __m512i t1_##j         = _mm512_permutexvar_epi16(idx##j, reshuffle_bits_lookup); \
        const __m512i cond_mask##j   = _mm512_permutexvar_epi16(idx##j, condition_mask_lookup); \
        const __m512i qcfail##j      = _mm512_srai_epi16(t1_##j, 16); \
        const __m512i t2_##j         = avx512_bitcond(epi16(0x00c0), data##j, (t0_##j | t1_##j)) & epi16(0x00ff); \
        const __m512i t3_##j         = _mm512_sllv_epi16(t2_##j, qcfail##j & _mm512_set1_epi16(8)); \
        qcfail_counter = _mm512_sub_epi16(qcfail_counter, qcfail##j); \
        funmap_counter_qcfail_0 = _mm512_add_epi16(funmap_counter_qcfail_0, \
                                          _mm512_ternarylogic_epi32(data##j, qcfail##j, _mm512_set1_epi16(FLAGSTAT_FUNMAP), 0x20)); \
        funmap_counter_qcfail_1 = _mm512_add_epi16(funmap_counter_qcfail_1, \
                                          _mm512_ternarylogic_epi32(data##j, qcfail##j, _mm512_set1_epi16(FLAGSTAT_FUNMAP), 0x80));

#define L(j) (t3_##j & cond_mask##j)

        for (/**/; i < thislimit; i += 16) {
            // 16-bit counters
            __m512i qcfail_counter = _mm512_setzero_si512();
            __m512i funmap_counter_qcfail_0 = _mm512_setzero_si512();
            __m512i funmap_counter_qcfail_1 = _mm512_setzero_si512();

#define U(pos) {                     \
    counter[pos] = _mm512_add_epi16(counter[pos], _mm512_and_si512(v16, one)); \
    v16 = _mm512_srli_epi16(v16, 1); \
}
            LOAD(0) LOAD(1)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 0),  L( 1));
            LOAD(2) LOAD(3)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L( 2),  L( 3));
            STORM_pospopcnt_csa_avx512(&foursA,  &v2,  twosA, twosB);
            LOAD(4) LOAD(5)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 4),  L( 5));
            LOAD(6) LOAD(7)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L( 6),  L( 7));
            STORM_pospopcnt_csa_avx512(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx512(&eightsA, &v4,  foursA,  foursB);
            LOAD(8) LOAD(9)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L( 8),   L( 9));
            LOAD(10) LOAD(11)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L(10),   L(11));
            STORM_pospopcnt_csa_avx512(&foursA,  &v2,  twosA,   twosB);
            LOAD(12) LOAD(13)
            STORM_pospopcnt_csa_avx512(&twosA,   &v1,  L(12),   L(13));
            LOAD(14) LOAD(15)
            STORM_pospopcnt_csa_avx512(&twosB,   &v1,  L(14),   L(15));
            STORM_pospopcnt_csa_avx512(&foursB,  &v2,  twosA,   twosB);
            STORM_pospopcnt_csa_avx512(&eightsB, &v4,  foursA,  foursB);
             U(0)  U(1)  U(2)  U(3)  U(4)  U(5)  U(6)  U(7)  U(8)  U(9)  U(10)  U(11)  U(12)  U(13)  U(14)  U(15) // Updates
            STORM_pospopcnt_csa_avx512(&v16,     &v8,  eightsA,  eightsB);

            qcfail_global_counter = _mm512_add_epi32(qcfail_global_counter,
                                                     _mm512_madd_epi16(qcfail_counter, one));
            funmap_global_counter_qcfail_0 =
                _mm512_add_epi32(funmap_global_counter_qcfail_0,
                _mm512_srli_epi16(_mm512_madd_epi16(funmap_counter_qcfail_0, one), 2));
            funmap_global_counter_qcfail_1 =
                _mm512_add_epi32(funmap_global_counter_qcfail_1,
                _mm512_srli_epi16(_mm512_madd_epi16(funmap_counter_qcfail_1, one), 2));
#undef U
#undef UU
#undef LOAD
#undef L
#undef LU
#undef W
#undef O1
#undef L1
#undef L2
#undef L3
        }

        // Update the counters after the last iteration
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_add_epi16(counter[i], _mm512_and_si512(v16, one));
            v16  = _mm512_srli_epi16(v16, 1);
        }

        for (size_t i = 0; i < 16; ++i) {
            _mm512_storeu_si512((__m512i*)buffer, counter[i]);
            for (size_t z = 0; z < 32; z++) {
                flags[i] += 16 * (uint32_t)buffer[z];
            }
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v1);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v2);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 2 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v4);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 4 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    _mm512_storeu_si512((__m512i*)buffer, v8);
    for (size_t i = 0; i < 32; ++i) {
        for (int j = 0; j < 16; ++j) {
            flags[j] += 8 * ((buffer[i] & (1 << j)) >> j);
        }
    }

    // update flags from our custom flags table
    user_flags[FLAGSTAT_FQCFAIL_OFF + 16]        += avx512_sum_epu32(qcfail_global_counter);
    user_flags[FLAGSTAT_FUNMAP_OFF + 0]          += avx512_sum_epu32(funmap_global_counter_qcfail_0);
    user_flags[FLAGSTAT_FUNMAP_OFF + 16]         += avx512_sum_epu32(funmap_global_counter_qcfail_1);
    user_flags[FLAGSTAT_FDUP_OFF + 0]            += flags[AVX512_FDUP_FQCFAIL_0];
    user_flags[FLAGSTAT_FDUP_OFF + 16]           += flags[AVX512_FDUP_FQCFAIL_1];
    user_flags[FLAGSTAT_FSECONDARY_OFF + 0]      += flags[AVX512_FSECONDARY_FQCFAIL_0];
    user_flags[FLAGSTAT_FSECONDARY_OFF + 16]     += flags[AVX512_FSECONDARY_FQCFAIL_1];
    user_flags[FLAGSTAT_FSUPPLEMENTARY_OFF + 0]  += flags[AVX512_FSUPPLEMENTARY_FQCFAIL_0];
    user_flags[FLAGSTAT_FSUPPLEMENTARY_OFF + 16] += flags[AVX512_FSUPPLEMENTARY_FQCFAIL_1];
    user_flags[FLAGSTAT_BIT12_OFF + 0]           += flags[AVX512_BIT12_FQCFAIL_0];
    user_flags[FLAGSTAT_BIT12_OFF + 16]          += flags[AVX512_BIT12_FQCFAIL_1];
    user_flags[FLAGSTAT_BIT13_OFF + 0]           += flags[AVX512_BIT13_FQCFAIL_0];
    user_flags[FLAGSTAT_BIT13_OFF + 16]          += flags[AVX512_BIT13_FQCFAIL_1];
    user_flags[FLAGSTAT_BIT14_OFF + 0]           += flags[AVX512_BIT14_FQCFAIL_0];
    user_flags[FLAGSTAT_BIT14_OFF + 16]          += flags[AVX512_BIT14_FQCFAIL_1];
    user_flags[FLAGSTAT_FREAD1_OFF + 0]          += flags[AVX512_FREAD1_FQCFAIL_0];
    user_flags[FLAGSTAT_FREAD1_OFF + 16]         += flags[AVX512_FREAD1_FQCFAIL_1];
    user_flags[FLAGSTAT_FREAD2_OFF + 0]          += flags[AVX512_FREAD2_FQCFAIL_0];
    user_flags[FLAGSTAT_FREAD2_OFF + 16]         += flags[AVX512_FREAD2_FQCFAIL_1];

    return 0;
}
#undef AVX512_BIT12_FQCFAIL_0
#undef AVX512_BIT12_FQCFAIL_1
#undef AVX512_BIT13_FQCFAIL_0
#undef AVX512_BIT13_FQCFAIL_1
#undef AVX512_BIT14_FQCFAIL_0
#undef AVX512_BIT14_FQCFAIL_1
#undef AVX512_FREAD1_FQCFAIL_0
#undef AVX512_FREAD2_FQCFAIL_0
#undef AVX512_FREAD1_FQCFAIL_1
#undef AVX512_FREAD2_FQCFAIL_1
#undef AVX512_FSECONDARY_FQCFAIL_0
#undef AVX512_FSECONDARY_FQCFAIL_1
#undef AVX512_FDUP_FQCFAIL_0
#undef AVX512_FDUP_FQCFAIL_1
#undef AVX512_FSUPPLEMENTARY_FQCFAIL_0
#undef AVX512_FSUPPLEMENTARY_FQCFAIL_1
#endif // end AVX512

/* *************************************
*  Function pointer definitions.
***************************************/
typedef int (*FLAGSTATS_func)(const uint16_t*, uint32_t, uint32_t*);

/* *************************************
*  Wrapper functions
***************************************/

static
FLAGSTATS_func FLAGSTATS_get_function(uint32_t n_len)
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
        return &FLAGSTAT_avx512;
    }
#endif

#if defined(STORM_HAVE_AVX2)
    if ((cpuid & STORM_CPUID_runtime_bit_AVX2) && n_len >= 512) { // 16*256
        return &FLAGSTAT_avx2;
    }
    
    if ((cpuid & STORM_CPUID_runtime_bit_SSE42) && n_len >= 256) { // 16*128
        return &FLAGSTAT_sse4;
    }
#endif

#if defined(STORM_HAVE_SSE42)
    if ((cpuid & STORM_CPUID_runtime_bit_SSE42) && n_len >= 256) { // 16*128
        return &FLAGSTAT_sse4;
    }
#endif

    return &FLAGSTAT_scalar;
}

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
