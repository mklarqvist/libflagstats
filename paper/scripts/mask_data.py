#!/usr/bin/env python3

FLAGSTAT_FPAIRED        =    1 # bit 0
FLAGSTAT_FUNMAP         =    4 # bit 2
FLAGSTAT_FMUNMAP        =    8 # bit 3
FLAGSTAT_FREVERSE       =   16 # bit 4
FLAGSTAT_FMREVERSE      =   32 # bit 5
FLAGSTAT_FREAD1         =   64 # bit 6
FLAGSTAT_FREAD2         =  128 # bit 7
FLAGSTAT_FSECONDARY     =  256 # bit 8
FLAGSTAT_FQCFAIL        =  512 # bit 9
FLAGSTAT_FDUP           = 1024 # bit 10
FLAGSTAT_FSUPPLEMENTARY = 2048 # bit 11

FLAGSTAT_BIT12 = 1 << 12
FLAGSTAT_BIT13 = 1 << 13
FLAGSTAT_BIT14 = 1 << 14

"""
      FLAGSTAT_FSUPPLEMENTARY
      | FLAGSTAT_FQCFAIL
      | |FLAGSTAT_FSECONDARY
      | ||         FLAGSTAT_FPAIRED
      | ||         |
[0000|x0yw|0000|000z]
[0000|0000|0000|xzyw]
"""

def get_mask(FSUPPLEMENTARY, FQCFAIL, FSECONDARY, FPAIRED):
    mask = 0
    mask |= FLAGSTAT_FUNMAP
    mask |= FLAGSTAT_FDUP
    mask |= FLAGSTAT_FQCFAIL

    if FSECONDARY:
        mask |= FLAGSTAT_FSECONDARY
    elif FSUPPLEMENTARY:
        mask |= FLAGSTAT_FSUPPLEMENTARY
    elif FPAIRED:
        mask |= FLAGSTAT_BIT12
        mask |= FLAGSTAT_BIT13
        mask |= FLAGSTAT_BIT14
        mask |= FLAGSTAT_FREAD1
        mask |= FLAGSTAT_FREAD2

    return mask

def qcfail(case):
    print(" # | FSUPPLEMENTARY | FPAIRED | FQCFAIL | FSECONDARY || mask")

    vpshufb_values = []
    for k in range(32):
        FSECONDARY     = int(k & 0x01 != 0)
        FQCFAIL        = int(k & 0x02 != 0)
        FPAIRED        = int(k & 0x04 != 0)
        FSUPPLEMENTARY = int(k & 0x08 != 0)

        m = get_mask(FSUPPLEMENTARY, FQCFAIL, FSECONDARY, FPAIRED)
        if FQCFAIL != case:
            m = 0

        vpshufb_values.append(m)

        print(f"{k:^3x}|{FSUPPLEMENTARY:^16}|{FPAIRED:^9}|{FQCFAIL:^9}|{FSECONDARY:^12}|| 0x{m:04x}")

    array = []
    for word in vpshufb_values:
        array.append(word & 0xff)
        array.append(word >> 8)

    print(avx512_const(array))

def avx512_const(array):
    assert len(array) == 64
    dwords = []
    for i in range(0, 64, 4):
        b0 = array[i + 0]
        b1 = array[i + 1]
        b2 = array[i + 2]
        b3 = array[i + 3]

        dword = (b3 << 24) | (b2 << 16) | (b1 << 8) | b0
        dwords.append(dword)

    return "_mm512_setr_epi32(%s)" % ', '.join('0x%08x' % v for v in dwords)

qcfail(1)
qcfail(0)
