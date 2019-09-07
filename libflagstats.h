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

#ifdef __cplusplus
extern "C" {
#endif

static const uint16_t FLAGSTAT_lookup_sec[2]  = {65535, 256}; // Complete selector, BAM_FSECONDARY
static const uint16_t FLAGSTAT_lookup_sup[2]  = {65535, 2048}; // Complete selector, BAM_FSUPPLEMENTARY
static const uint16_t FLAGSTAT_lookup_pair[2] = {256+2048, 65535}; // BAM_FSECONDARY + BAM_FSUPPLEMENTARY, complete selector

void FLAGSTAT_samtools_single_update(uint16_t val, uint32_t* flags) {
    const int offset = ( (val & BAM_FQCFAIL) == 0 ) ? 0 : 16;

    // Always.
    flags[offset+3]  += (val & BAM_FUNMAP) == 0; // this is implicit
    flags[offset+10] += (val & BAM_FDUP) >> 10;
    // Rule 1.
    flags[offset+8]  += (val & BAM_FSECONDARY) >> 8;
    val &= FLAGSTAT_lookup_sec[(val & BAM_FSECONDARY) >> 8];
    // Rule 2.
    flags[offset+11] += (val & BAM_FSUPPLEMENTARY) >> 11;
    val &= FLAGSTAT_lookup_sup[(val & BAM_FSUPPLEMENTARY) >> 11];
    val &= FLAGSTAT_lookup_pair[val & BAM_FPAIRED];
    
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

#if defined(STORM_HAVE_AVX2)

#include <immintrin.h>

STORM_TARGET("avx2")
static
int FLAGSTAT_avx2(const uint16_t* array, uint32_t len, uint32_t* flags) {
    for (uint32_t i = len - (len % (16 * 16)); i < len; ++i) {
        FLAGSTAT_samtools_single_update(array[i], flags);
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
    const __m256i one  = _mm256_set1_epi16(1); // (11...1) vector
    const __m256i zero = _mm256_set1_epi16(0); // (00...0) vector

    // Mask selectors.
    const __m256i mask1 = _mm256_set1_epi16(256 + 2048); // BAM_FSECONDARY (256) + BAM_FSUPPLEMENTARY (2048)
    const __m256i mask2 = _mm256_set1_epi16(4 + 256 + 1024 + 2048); // Need to keep BAM_FUNMAP (4) + BAM_FDUP (1024) 
        // for Always rule and BAM_FSECONDARY (256) and BAM_FSUPPLEMENTARY (2048) for Rule 1 and Rule 2.
    const __m256i mask3 = _mm256_set1_epi16(512); // BAM_FQCFAIL (512)

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
        // We load a register of data (data + i + j) and then using a the resulting mask from
        // a VPCMPEQW instruction comparing equality with the mask mask1 
        // (BAM_FSECONDARY + BAM_FSUPPLEMENTARY). The resulting data is either the original data
        // or empty as DATA & (00...0) is a zero register and DATA & (11...1) is the data itself. 
        // The resulting data is combined (bitwise or) with the mask mask2 as this information
        // is required.
        //
        // FLAGSTATS outputs summary statistics separately for reads that pass QC and those
        // that do not. Therefore we need to partition the data into these two classes.
        // For data that pass QC, the L registers, we first bit-select the target BAM_FQCFAIL bit
        // using the mask mask3. The resulting data is used to perform another mask-select using VPCMPEQW 
        // against the empty vector (00...0). As above, if the data has the BAM_FQCFAIL bit set then
        // this register will be zeroed out.
        // The exact process is performed for reads that fail QC, the LU registers, with the difference
        // that mask-selection is based on the one vector (11...1).

#define LOAD(j) __m256i data##j = _mm256_loadu_si256(data + i + j) & (_mm256_cmpeq_epi16( _mm256_loadu_si256(data + i + j) & mask1, zero ) | mask2);
#define L(j)  data##j & _mm256_cmpeq_epi16( data##j & mask3, zero )
#define LU(j) data##j & _mm256_cmpeq_epi16( data##j & mask3, mask3 )
        
        for (/**/; i < thislimit; i += 16) {
#define U(pos) {                     \
    counter[pos] = _mm256_add_epi16(counter[pos], _mm256_and_si256(v16, one)); \
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

    return 0;
}
#endif // end avx2

#if defined(STORM_HAVE_AVX512)

#include <immintrin.h>

STORM_TARGET("avx512bw")
static
int FLAGSTAT_avx512(const uint16_t* array, size_t len, uint32_t* out) {
    for (uint32_t i = len - (len % (32 * 16)); i < len; ++i) {
        FLAGSTAT_samtools_single_update(array[i], out);
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
    const __m512i mask1 = _mm512_set1_epi16(256 + 2048); // BAM_FSECONDARY (256) + BAM_FSUPPLEMENTARY (2048)
    const __m512i mask2 = _mm512_set1_epi16(4 + 256 + 1024 + 2048); // Need to keep BAM_FUNMAP (4) + BAM_FDUP (1024) 
        // for Always rule and BAM_FSECONDARY (256) and BAM_FSUPPLEMENTARY (2048) for Rule 1 and Rule 2.
    const __m512i mask3 = _mm512_set1_epi16(512); // BAM_FQCFAIL (512)

    while (i < limit) {
        for (size_t i = 0; i < 16; ++i) {
            counter[i]  = _mm512_setzero_si512();
            counterU[i] = _mm512_setzero_si512();
        }

        size_t thislimit = limit;
        if (thislimit - i >= (1 << 16))
            thislimit = i + (1 << 16) - 1;

        ///////////////////////////////////////////////////////////////////////
        // We load a register of data (data + i + j) and then using a the resulting mask from
        // a VPCMPEQW instruction comparing equality with the mask mask1 
        // (BAM_FSECONDARY + BAM_FSUPPLEMENTARY). The resulting data is either the original data
        // or empty as DATA & (00...0) is a zero register and DATA & (11...1) is the data itself. 
        // The resulting data is combined (bitwise or) with the mask mask2 as this information
        // is required.
        //
        // FLAGSTATS outputs summary statistics separately for reads that pass QC and those
        // that do not. Therefore we need to partition the data into these two classes.
        // For data that pass QC, the L registers, we first bit-select the target BAM_FQCFAIL bit
        // using the mask mask3. The resulting data is used to perform another mask-select using VPCMPEQW 
        // against the empty vector (00...0). As above, if the data has the BAM_FQCFAIL bit set then
        // this register will be zeroed out.
        // The exact process is performed for reads that fail QC, the LU registers, with the difference
        // that mask-selection is based on the one vector (11...1).

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


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif