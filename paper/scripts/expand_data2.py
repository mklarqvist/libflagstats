#!/usr/bin/env python3

def bit12(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP):
    return FPAIRED & (~FUNMAP) & FPROPER_PAIR

def bit13(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP):
    return FPAIRED & (~FUNMAP) & FMUNMAP

def bit14(FPAIRED, FPROPER_PAIR, FUNMAP, FMUNMAP):
    return FPAIRED & (~FUNMAP) & (~FMUNMAP)

print(" # | MNUNMAP | FUNMAP | FPROPER_PAIR | FPAIRED || bit #12 | bit #13 | bit #14 | pshufb word")
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
    
    pshufb = (b12 << 4) | (b13 << 5) | (b14 << 6) | (FPAIRED << 2)
    vpshufb_values.append(0)
    vpshufb_values.append(pshufb)

    print(f"{k:^3x}|{FMUNMAP:^9}|{FUNMAP:^8}|{FPROPER_PAIR:^14}|{FPAIRED:^9}||{b12:^9}|{b13:^9}|{b14:^9}| 0x{pshufb:02x}")

vpshufb_values.extend(vpshufb_values[:]) # the upper half must be the same as lower one

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

print(avx512_const(vpshufb_values))
