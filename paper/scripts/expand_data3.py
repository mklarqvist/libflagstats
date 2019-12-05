#!/usr/bin/env python3
from avx512 import *

AVX512_BIT12_FQCFAIL_0          = 0
AVX512_BIT12_FQCFAIL_1          = 1
AVX512_BIT13_FQCFAIL_0          = 2
AVX512_BIT13_FQCFAIL_1          = 3
AVX512_BIT14_FQCFAIL_0          = 4
AVX512_BIT14_FQCFAIL_1          = 5
AVX512_FREAD1_FQCFAIL_0         = 6
AVX512_FREAD1_FQCFAIL_1         = 8
AVX512_FREAD2_FQCFAIL_0         = 7
AVX512_FREAD2_FQCFAIL_1         = 9
AVX512_FSECONDARY_FQCFAIL_0     = 10
AVX512_FSECONDARY_FQCFAIL_1     = 11
AVX512_FSUPPLEMENTARY_FQCFAIL_0 = 14
AVX512_FSUPPLEMENTARY_FQCFAIL_1 = 15
AVX512_FDUP_FQCFAIL_0           = 12
AVX512_FDUP_FQCFAIL_1           = 13

def bit12(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP):
    return FPAIRED & (~FUNMAP) & FPROPER_PAIR

def bit13(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP):
    return FPAIRED & (~FUNMAP) & FMUNMAP

def bit14(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP):
    return FPAIRED & (~FUNMAP) & (~FMUNMAP)

print(" # | MNUNMAP | FUNMAP | FPROPER_PAIR | FPAIRED || bit #12 | bit #13 | bit #14 | pshufw word")
print("---+---------+--------+--------------+---------++---------+---------+---------+-------------")

vpshufb_values = []
for k in range(16):
    FPAIRED      = int(k & 0x01 != 0)
    FPROPER_PAIR = int(k & 0x02 != 0)
    FUNMAP       = int(k & 0x04 != 0)
    FMUNMAP      = int(k & 0x08 != 0)

    b12 = bit12(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP) & 0x01
    b13 = bit13(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP) & 0x01
    b14 = bit14(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP) & 0x01

    word = (b12 << AVX512_BIT12_FQCFAIL_0) \
         | (b12 << AVX512_BIT12_FQCFAIL_1) \
         | (b13 << AVX512_BIT13_FQCFAIL_0) \
         | (b13 << AVX512_BIT13_FQCFAIL_1) \
         | (b14 << AVX512_BIT14_FQCFAIL_0) \
         | (b14 << AVX512_BIT14_FQCFAIL_1) \
         | (FPAIRED << 12) \
         | (FUNMAP << 9)

    vpshufb_values.append(word & 0xff)
    vpshufb_values.append(word >> 8)

    print(f"{k:^3x}|{FMUNMAP:^9}|{FUNMAP:^8}|{FPROPER_PAIR:^14}|{FPAIRED:^9}||{b12:^9}|{b13:^9}|{b14:^9}| 0x{word:04x}")

vpshufb_values.extend(vpshufb_values[:]) # the upper half must be the same as lower one
print(avx512_const(vpshufb_values))
