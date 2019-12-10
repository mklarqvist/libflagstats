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

def get_mask16(FSECONDARY, FSUPPLEMENTARY, FPAIRED):
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


def get_mask(FSUPPLEMENTARY, HIGH_BYTE, FPAIRED, FSECONDARY):
    mask = get_mask16(FSECONDARY, FSUPPLEMENTARY, FPAIRED)

    if HIGH_BYTE:
        return mask >> 8
    else:
        return mask & 0xff


def main():
    values = []
    for k in range(2**4):
        FSECONDARY      = int(k & 0x01 != 0)
        FPAIRED         = int(k & 0x02 != 0)
        HIGH_BYTE       = int(k & 0x04 != 0)
        FSUPPLEMENTARY  = int(k & 0x08 != 0)

        values.append(get_mask(FSUPPLEMENTARY, HIGH_BYTE, FPAIRED, FSECONDARY))

    hexstr = ', '.join(f"0x{byte:02x}" for byte in values)
    print(f"const __m128i mask_lookup = _mm_setr_epi8({hexstr});")
    print(f"const __m256i mask_lookup = _mm256_setr_epi8({hexstr}, {hexstr});")


if __name__ == '__main__':
    main()
