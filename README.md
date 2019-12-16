# libflagstats

[![PyPI version](https://badge.fury.io/py/pyflagstats.svg)](https://badge.fury.io/py/pyflagstats)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mklarqvist/libflagstats/blob/master/python/pyflagstats.ipynb)


These functions compute the summary statistics for the SAM FLAG field (flagstat)
using fast SIMD instructions. They are based on the
[positional-popcount](https://github.com/mklarqvist/positional-popcount)
(`pospopcnt`) subroutines implemented in
[libalgebra](https://github.com/mklarqvist/libalgebra) and assumes that the data
is available as contiguous streams (e.g. [column projection](https://en.wikipedia.org/wiki/Projection_(relational_algebra))). The heavily
branched and bit-dependent [SAMtools](https://github.com/samtools/samtools/) code can be rewritten using a mask-select
propagate-carry approach that feeds into bit-transposed Harley-Seal-based carry-save adder
networks.

The positional population count is described in the paper:
* [Efficient Computation of Positional Population Counts Using SIMD Instructions](https://arxiv.org/abs/1911.02696) by Marcus D. R. Klarqvist and Wojciech Muła and Daniel Lemire

## Speedup

This benchmark shows the speedup of the pospopcnt algorithms compared to
[samtools](https://github.com/samtools/samtools) using a human HiSeqX readset
with 824,541,892 reads. See [Results](#results) for additional information. On
this readset, the pospopcnt-based functions are >3802-fold faster compared to
BAM and 402.6-fold faster compared to CRAM. Note that the BAM format does not
have column projection capabilities. Around 80% of the CPU time is spent
retrieving data from disk for the `pospopcnt`-functions. In this example, we
compress data blocks of 512k records (1 MB) using LZ4. 

| Approach        | Time       | Speedup |
|-----------------|------------|---------|
| Samtools – BAM  | 30m 50.06s | 1       |
| Samtools – CRAM | 4m 50.68s  | 6.36    |
| libflagstats-LZ9       | 0.72s      | 2569.53 |
| libflagstats-raw       | 0.48s      | 3802.90 |

We also directly compared the potential speed of the naive flagstat subroutine
in samtools againt these functions if samtools would be rewritten with efficient
column projection and the compression method changed to LZ4. In this context,
these functions are still 6.58-fold faster.

| Approach             | Time   | Speedup |
|----------------------|--------|---------|
| samtools-rewrite+LZ4 | 4.74 s | 1       |
| libflagstats            | 0.72s  | 6.58    |

### Usage

For Linux/Mac: Compile the test suite with: `make`. LZ4 and Zstd needs to be
installed on the target machine or pass the `LZ4_PATH` and `ZSTD_PATH` flags to
make for non-standard locations.

First, import 2-byte FLAG words from a SAM file:
```bash
samtools view NA12878D_HiSeqX_R12_GRCh37.bam | cut -f 2 | ./utility > NA12878D_HiSeqX_R12_GRCh37_flags.bin
```

Compress the binary file using 512k blocks:
```bash
for i in {1..10}; do ./bench compress -i NA12878D_HiSeqX_R12_GRCh37_flags.bin -l -c $i; done
for i in {1..10}; do ./bench compress -i NA12878D_HiSeqX_R12_GRCh37_flags.bin -l -f -c $i; done
for i in {1..20}; do ./bench compress -i NA12878D_HiSeqX_R12_GRCh37_flags.bin -z -c $i; done
```

Evaluate the flagstat subroutines by first decompressing the file twice while
clearing cache and then computing samtools and pospopcnt.
```bash
for i in {1..10}; do ./bench decompress -i NA12878D_HiSeqX_R12_GRCh37_flags.bin_HC_c${i}.lz4; done
for i in {1..10}; do ./bench decompress -i NA12878D_HiSeqX_R12_GRCh37_flags.bin_fast_a${i}.lz4; done
for i in {1..20}; do ./bench decompress -i NA12878D_HiSeqX_R12_GRCh37_flags.bin_c${i}.zst; done
```

### History

These functions were developed for [pil](https://github.com/mklarqvist/pil).

## Problem statement

The FLAG field in the [SAM interchange
format](https://github.com/samtools/hts-specs) is defined as the union of
[1-hot](https://en.wikipedia.org/wiki/One-hot) encoded states for a given read.
For example, the following three states evaluating to true

```
00000001: read paired
01000000: first in pair
00001000: mate unmapped
--------
01001001: Decimal (73)
```

are stored in a packed 16-bit value (only the LSB is shown here). There are 12
states described in the SAM format:

| One-hot           | Description                               |
|-------------------|-------------------------------------------|
| 00000000 00000001 | Read paired                               |
| 00000000 00000010 | Read mapped in proper pair                |
| 00000000 00000100 | Read unmapped                             |
| 00000000 00001000 | Mate unmapped                             |
| 00000000 00010000 | Read reverse strand                       |
| 00000000 00100000 | Mate reverse strand                       |
| 00000000 01000000 | First in pair                             |
| 00000000 10000000 | Second in pair                            |
| 00000001 00000000 | Not primary alignment                     |
| 00000010 00000000 | Read fails platform/vendor quality checks |
| 00000100 00000000 | Read is PCR or optical duplicate          |
| 00001000 00000000 | Supplementary alignment                   |

Computing FLAG statistics from readsets involves iteratively incrementing up to
16 counters. The native properties of a column-oriented storage, specifically
column projection, already deliver good performance because of data locality
(memory contiguity) and value typing. We want to maximize compute on large
arrays of values by exploiting vectorized instructions, if available.

## Goals

* Achieve high-performance on large arrays of values.
* Support machines without SIMD (scalar).
* Specialized algorithms for SSE4 up to AVX512.

## Technical approach

### Results
### Datasets
Aligned data:
* https://dnanexus-rnd.s3.amazonaws.com/NA12878-xten/mappings/NA12878D_HiSeqX_R1.bam

Unaligned data:
* https://dnanexus-rnd.s3.amazonaws.com/NA12878-xten/reads/NA12878D_HiSeqX_R1.fastq.gz
* https://dnanexus-rnd.s3.amazonaws.com/NA12878-xten/reads/NA12878D_HiSeqX_R2.fastq.gz

### Speed

| Comp. Method | Decomp. | Sam-branchless | sam  | flagstat |
|--------------|---------|----------------|------|----------|
| LZ4-HC-c1    | 988     | 8924           | 4991 | 1107     |
| LZ4-HC-c2    | 993     | 8848           | 5076 | 1132     |
| LZ4-HC-c3    | 938     | 8686           | 4930 | 1049     |
| LZ4-HC-c4    | 846     | 8803           | 4876 | 933      |
| LZ4-HC-c5    | 824     | 8525           | 5117 | 974      |
| LZ4-HC-c6    | 770     | 8536           | 4774 | 851      |
| LZ4-HC-c7    | 680     | 8404           | 4748 | 837      |
| LZ4-HC-c8    | 644     | 8453           | 4662 | 755      |
| LZ4-HC-c9    | **580**     | 8434           | 4740 | **722**      |
| LZ4-fast-c2  | 814     | 8658           | 4886 | 990      |
| LZ4-fast-c3  | 837     | 8576           | 4840 | 941      |
| LZ4-fast-c4  | 889     | 8627           | 4861 | 1026     |
| LZ4-fast-c5  | 826     | 8590           | 4885 | 1037     |
| LZ4-fast-c6  | 823     | 8629           | 5034 | 951      |
| LZ4-fast-c7  | 837     | 8606           | 4999 | 985      |
| LZ4-fast-c8  | 834     | 8604           | 4920 | 962      |
| LZ4-fast-c9  | 853     | 8615           | 4944 | 951      |
| LZ4-fast-c10 | 853     | 8648           | 5024 | 950      |
| Zstd-c1      | 3435    | 11438          | 7798 | 3630     |
| Zstd-c2      | 3577    | 11231          | 8110 | 3767     |
| Zstd-c3      | 3403    | 11250          | 7922 | 3553     |
| Zstd-c4      | 3562    | 11223          | 7949 | 3649     |
| Zstd-c5      | 2919    | 10584          | 7263 | 2986     |
| Zstd-c6      | 2964    | 10680          | 7545 | 3015     |
| Zstd-c7      | 2681    | 10591          | 7067 | 2715     |
| Zstd-c8      | 2641    | 10523          | 7103 | 2850     |
| Zstd-c9      | 2352    | 10453          | 6669 | 2463     |
| Zstd-c10     | 2309    | 10094          | 6756 | 2509     |
| Zstd-c11     | 2344    | 10018          | 6430 | 2467     |
| Zstd-c12     | 2116    | 9916           | 6242 | 2252     |
| Zstd-c13     | 2107    | 9844           | 6183 | 2236     |
| Zstd-c14     | 1955    | 9616           | 5969 | 2044     |
| Zstd-c15     | 1716    | 9562           | 5807 | 1808     |
| Zstd-c16     | 1286    | 9208           | 5592 | 1448     |
| Zstd-c17     | 1278    | 8996           | 5592 | 1396     |
| Zstd-c18     | 1192    | 8907           | 5294 | 1306     |
| Zstd-c19     | 1181    | 8931           | 5362 | 1293     |
| Zstd-c20     | 1175    | 8982           | 5369 | 1303     |

```bash
$ time samtools flagstat NA12878D_HiSeqX_R12_GRCh37.bam
824541892 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
5393628 + 0 supplementary
0 + 0 duplicates
805383403 + 0 mapped (97.68% : N/A)
819148264 + 0 paired in sequencing
409574132 + 0 read1
409574132 + 0 read2
781085884 + 0 properly paired (95.35% : N/A)
797950890 + 0 with itself and mate mapped
2038885 + 0 singletons (0.25% : N/A)
9537902 + 0 with mate mapped to a different chr
4425946 + 0 with mate mapped to a different chr (mapQ>=5)

real    30m50.059s
user    30m10.638s
sys 0m38.440s
```

```bash
$ time samtools flagstat NA12878D_HiSeqX_R12_GRCh37.cram
824541892 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
5393628 + 0 supplementary
0 + 0 duplicates
805383403 + 0 mapped (97.68% : N/A)
819148264 + 0 paired in sequencing
409574132 + 0 read1
409574132 + 0 read2
781085884 + 0 properly paired (95.35% : N/A)
797950890 + 0 with itself and mate mapped
2038885 + 0 singletons (0.25% : N/A)
9537902 + 0 with mate mapped to a different chr
4425946 + 0 with mate mapped to a different chr (mapQ>=5)

real    4m50.684s
user    3m37.394s
sys 1m12.396s
```

```c
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
```