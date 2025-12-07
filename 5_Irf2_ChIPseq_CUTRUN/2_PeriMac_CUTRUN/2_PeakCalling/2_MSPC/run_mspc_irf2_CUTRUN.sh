/Users/berenice/Softwares/mspc_5.5/mspc -f ./Irf2_CUTandRUN/   -c 4 -r bio -w 1e-5 -s 1e-8 -m Lowest --excludeHeader -o MSPC_Irf2_PeriMac_4samples

cp ./MSPC_Irf2_PeriMac_4samples/ConsensusPeaks.bed     MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.bed  

wc -l *_ConsensusPeaks.bed
#    1565 MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.bed

annotatePeaks.pl MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.bed   mm10 > HOMER_MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.xls  


########
intersectBed -a MSPC_Irf2_PeriMac_4samples_ConsensusPeaks.bed -b /Users/berenice/Desktop/Irf2_CUTRUN/01_ChIPseq/MSPC/MSPC_Irf2_BMDM_ConsensusPeaks.bed > MSPC_Irf2_PeriMac_BMDM_intersect.bed
#     316 MSPC_Irf2_PeriMac_BMDM_intersect.bed

annotatePeaks.pl MSPC_Irf2_PeriMac_BMDM_intersect.bed   mm10 > HOMER_MSPC_Irf2_PeriMac_BMDM_intersect.xls  

intersectBed -a MSPC_Irf2_PeriMac_BMDM_intersect.bed -b /Users/berenice/Desktop/Irf2_CUTRUN/DiffBind/MSPC_4samples/2025-11-28_PeriMac_Irf2_stringent_Aging_diffbind_peaks_FDR5.bed > MSPC_Irf2_PeriMac_BMDM_AgingFDR5_intersect.bed
#     105 MSPC_Irf2_PeriMac_BMDM_AgingFDR5_intersect.bed

annotatePeaks.pl MSPC_Irf2_PeriMac_BMDM_AgingFDR5_intersect.bed   mm10 > HOMER_MSPC_Irf2_PeriMac_BMDM_AgingFDR5_intersect.xls  
