# FlagStats

| Comp. Method | Decomp. | samtools* | flagstat |
|--------------|---------|-----------|----------|
| LZ4-HC-c1    | 988     | 9335      | 1107     |
| LZ4-HC-c2    | 993     | 8932      | 1132     |
| LZ4-HC-c3    | 938     | 8757      | 1049     |
| LZ4-HC-c4    | 846     | 9110      | 933      |
| LZ4-HC-c5    | 824     | 9044      | 974      |
| LZ4-HC-c6    | 770     | 9280      | 851      |
| LZ4-HC-c7    | 680     | 8828      | 837      |
| LZ4-HC-c8    | 644     | 9083      | 755      |
| LZ4-HC-c9    | 580     | 8455      | 722      |
| LZ4-fast-c2  | 814     | 8765      | 990      |
| LZ4-fast-c3  | 837     | 8758      | 941      |
| LZ4-fast-c4  | 889     | 8608      | 1026     |
| LZ4-fast-c5  | 826     | 9366      | 1037     |
| LZ4-fast-c6  | 823     | 8843      | 951      |
| LZ4-fast-c7  | 837     | 8781      | 985      |
| LZ4-fast-c8  | 834     | 8738      | 962      |
| LZ4-fast-c9  | 853     | 8979      | 951      |
| LZ4-fast-c10 | 853     | 8982      | 950      |
| Zstd-c1      | 3435    | 7661      | 3630     |
| Zstd-c2      | 3577    | 7763      | 3767     |
| Zstd-c3      | 3403    | 7943      | 3553     |
| Zstd-c4      | 3562    | 7906      | 3649     |
| Zstd-c5      | 2919    | 6970      | 2986     |
| Zstd-c6      | 2964    | 7030      | 3015     |
| Zstd-c7      | 2681    | 7128      | 2715     |
| Zstd-c8      | 2641    | 7027      | 2850     |
| Zstd-c9      | 2352    | 6762      | 2463     |
| Zstd-c10     | 2309    | 6497      | 2509     |
| Zstd-c11     | 2344    | 6665      | 2467     |
| Zstd-c12     | 2116    | 6377      | 2252     |
| Zstd-c13     | 2107    | 6494      | 2236     |
| Zstd-c14     | 1955    | 6342      | 2044     |
| Zstd-c15     | 1716    | 5901      | 1808     |
| Zstd-c16     | 1286    | 5906      | 1448     |
| Zstd-c17     | 1278    | 5547      | 1396     |
| Zstd-c18     | 1192    | 5402      | 1306     |
| Zstd-c19     | 1181    | 5457      | 1293     |
| Zstd-c20     | 1175    | 5433      | 1303     |

| Comp. Method | flagstat | samtools* |
|--------------|----------|-----------|
| LZ4-HC-c1    | 744.84   | 88.33     |
| LZ4-HC-c2    | 728.39   | 92.31     |
| LZ4-HC-c3    | 786.03   | 94.16     |
| LZ4-HC-c4    | 883.75   | 90.51     |
| LZ4-HC-c5    | 846.55   | 91.17     |
| LZ4-HC-c6    | 968.91   | 88.85     |
| LZ4-HC-c7    | 985.12   | 93.4      |
| LZ4-HC-c8    | 1092.11  | 90.78     |
| LZ4-HC-c9    | 1026.83  | 94.33     |
| LZ4-fast-c2  | 832.87   | 94.07     |
| LZ4-fast-c3  | 883.75   | 95.64     |
| LZ4-fast-c4  | 803.65   | 95.79     |
| LZ4-fast-c5  | 795.12   | 88.04     |
| LZ4-fast-c6  | 827.02   | 92.58     |
| LZ4-fast-c7  | 837.1    | 93.9      |
| LZ4-fast-c8  | 857.11   | 94.36     |
| LZ4-fast-c9  | 867.03   | 91.83     |
| LZ4-fast-c10 | 867.94   | 91.8      |
| Zstd-c1      | 219.23   | 106.2     |
| Zstd-c2      | 218.89   | 106.21    |
| Zstd-c3      | 232.07   | 103.81    |
| Zstd-c4      | 225.96   | 104.29    |
| Zstd-c5      | 276.14   | 118.3     |
| Zstd-c6      | 273.48   | 117.29    |
| Zstd-c7      | 303.7    | 115.68    |
| Zstd-c8      | 289.31   | 117.34    |
| Zstd-c9      | 334.77   | 121.94    |
| Zstd-c10     | 328.63   | 126.91    |
| Zstd-c11     | 334.23   | 123.71    |
| Zstd-c12     | 366.14   | 129.3     |
| Zstd-c13     | 368.76   | 126.97    |
| Zstd-c14     | 403.4    | 130.01    |
| Zstd-c15     | 456.05   | 139.73    |
| Zstd-c16     | 569.44   | 139.61    |
| Zstd-c17     | 590.65   | 148.65    |
| Zstd-c18     | 631.35   | 152.64    |
| Zstd-c19     | 637.7    | 151.1     |
| Zstd-c20     | 632.8    | 151.77    |

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

real	30m50.059s
user	30m10.638s
sys	0m38.440s
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

real	4m50.684s
user	3m37.394s
sys	1m12.396s
```