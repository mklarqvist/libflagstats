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
AVX512_FUNMAP = 0

def bit(pos):
    return 1 << pos


def get_condition_mask(FSUPPLEMENTARY, FSECONDARY, FPAIRED):
    mask = 0
    mask |= bit(AVX512_FDUP_FQCFAIL_0) | bit(AVX512_FDUP_FQCFAIL_1)

    if FSECONDARY:
        mask |= bit(AVX512_FSECONDARY_FQCFAIL_0) | bit(AVX512_FSECONDARY_FQCFAIL_1)
    elif FSUPPLEMENTARY:
        mask |= bit(AVX512_FSUPPLEMENTARY_FQCFAIL_0) | bit(AVX512_FSUPPLEMENTARY_FQCFAIL_1)
    elif FPAIRED:
        mask |= bit(AVX512_BIT12_FQCFAIL_0) | bit(AVX512_BIT12_FQCFAIL_1)
        mask |= bit(AVX512_BIT13_FQCFAIL_0) | bit(AVX512_BIT13_FQCFAIL_1)
        mask |= bit(AVX512_BIT14_FQCFAIL_0) | bit(AVX512_BIT14_FQCFAIL_1)
        mask |= bit(AVX512_FREAD1_FQCFAIL_0) | bit(AVX512_FREAD1_FQCFAIL_1)
        mask |= bit(AVX512_FREAD2_FQCFAIL_0) | bit(AVX512_FREAD2_FQCFAIL_1)

    print("{:d} {:d} {:d} {:016b}".format(FSECONDARY, FSUPPLEMENTARY, FPAIRED, mask))

    return mask


def get_duplication_word(FSUPPLEMENTARY, FDUP, FSECONDARY, FUNMAP):
    mask = 0

    if FSUPPLEMENTARY:
        mask |= bit(AVX512_FSUPPLEMENTARY_FQCFAIL_0)
        mask |= bit(AVX512_FSUPPLEMENTARY_FQCFAIL_1)

    if FDUP:
        mask |= bit(AVX512_FDUP_FQCFAIL_0)
        mask |= bit(AVX512_FDUP_FQCFAIL_1)

    if FSECONDARY:
        mask |= bit(AVX512_FSECONDARY_FQCFAIL_0)
        mask |= bit(AVX512_FSECONDARY_FQCFAIL_1)

    if FUNMAP:
        mask |= bit(AVX512_FUNMAP)

    return mask


condition = []
duplication = []
for k in range(2**5):
    FSECONDARY     = int(k & 0x01 != 0)
    FUNMAP         = int(k & 0x02 != 0)
    FDUP           = int(k & 0x04 != 0)
    FSUPPLEMENTARY = int(k & 0x08 != 0)
    FPAIRED        = int(k & 0x10 != 0)

    condition.append(get_condition_mask(FSUPPLEMENTARY, FSECONDARY, FPAIRED))
    duplication.append(get_duplication_word(FSUPPLEMENTARY, FDUP, FSECONDARY, FUNMAP))

print("Duplication lookup")
print(avx512_const(word2byte_array(duplication)))

print()

print("Condition mask")
print(avx512_const(word2byte_array(condition)))
