from libc.stdint cimport uint64_t, uint32_t, uint16_t
import numpy as np
cimport numpy as np

cdef extern from "libflagstats.h":
    uint64_t FLAGSTATS_u16(const uint16_t* array, uint32_t n_len, uint32_t* flags)

def flagstats(values):
    if type(values) != np.ndarray:
        raise ValueError("Values must be an numpy.ndarray")

    if values.dtype != "uint16":
        raise ValueError("Values must have the dtype \"uint16\"")

    if not values.flags['C_CONTIGUOUS']:
        print("Input array is not contiguous. Fixing...")
        values = np.ascontiguousarray(values, dtype=np.uint16)
    
    cdef np.ndarray flags  = np.zeros(32, dtype="uint32")
    cdef uint32_t[::1] out = flags
    cdef uint16_t[::1] v   = values
    cdef ret = FLAGSTATS_u16(&v[0], len(values), &out[0])

    SAM_FLAG_NAMES = ["FPAIRED","FPROPER_PAIR","FUNMAP","FMUNMAP","FREVERSE","FMREVERSE", "FREAD1","FREAD2","FSECONDARY","FQCFAIL","FDUP","FSUPPLEMENTARY","n_pair_good","n_sgltn","n_pair_map"]

    ret = {
    "n_values": len(values),
    "passed": 
        dict(zip(SAM_FLAG_NAMES, flags[0:15,])), 
    "failed": 
        dict(zip(SAM_FLAG_NAMES, flags[16:31,])), 
    }

    ret["passed"]["mapped"] = len(values) - ret["passed"]["FUNMAP"] - ret["failed"]["FUNMAP"]
    ret["passed"]["paired_in_seq"] = ret["passed"]["FREAD1"] + ret["passed"]["FREAD2"]

    return ret

