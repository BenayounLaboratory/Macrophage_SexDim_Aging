findPeaks ../BAM_clean/01_0LRG_02ZNUSC_YF-Pool1_IRF2_mm_i1-13.mm10_CLEAN_TAGs  -o 01_0LRG_02ZNUSC_YF-Pool1_IRF2.pos    -style factor
findPeaks ../BAM_clean/02_0LRH_02ZNUSC_YF-Pool2_IRF2_mm_i2-13.mm10_CLEAN_TAGs  -o 02_0LRH_02ZNUSC_YF-Pool2_IRF2.pos    -style factor
findPeaks ../BAM_clean/03_0LRI_02ZNUSC_OF-Pool1_IRF2_mm_i3-13.mm10_CLEAN_TAGs  -o 03_0LRI_02ZNUSC_OF-Pool1_IRF2.pos    -style factor
findPeaks ../BAM_clean/04_0LRJ_02ZNUSC_OF-Pool2_IRF2_mm_i4-13.mm10_CLEAN_TAGs  -o 04_0LRJ_02ZNUSC_OF-Pool2_IRF2.pos    -style factor

parse_pos_for_MSPC.pl 01_0LRG_02ZNUSC_YF-Pool1_IRF2.pos
parse_pos_for_MSPC.pl 02_0LRH_02ZNUSC_YF-Pool2_IRF2.pos
parse_pos_for_MSPC.pl 03_0LRI_02ZNUSC_OF-Pool1_IRF2.pos
parse_pos_for_MSPC.pl 04_0LRJ_02ZNUSC_OF-Pool2_IRF2.pos


#    4552 01_0LRG_02ZNUSC_YF-Pool1_IRF2.cl.bed
#    2344 02_0LRH_02ZNUSC_YF-Pool2_IRF2.cl.bed
#    2287 03_0LRI_02ZNUSC_OF-Pool1_IRF2.cl.bed
#    2932 04_0LRJ_02ZNUSC_OF-Pool2_IRF2.cl.bed
#   12115 total