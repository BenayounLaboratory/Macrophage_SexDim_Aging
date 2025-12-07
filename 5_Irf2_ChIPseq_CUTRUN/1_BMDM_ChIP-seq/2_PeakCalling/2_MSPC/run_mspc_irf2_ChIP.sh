/Users/berenice/Softwares/mspc_5.5/mspc -f ./Irf2_ChIP/   -c 2 -r bio -w 1e-5 -s 1e-8 -m Lowest --excludeHeader -o MSPC_Irf2_BMDM

cp ./MSPC_Irf2_BMDM/ConsensusPeaks.bed     MSPC_Irf2_BMDM_ConsensusPeaks.bed  

wc -l *_ConsensusPeaks.bed
#         460 /Users/berenice/Desktop/Irf2_CUTRUN/01_ChIPseq/MSPC/MSPC_Irf2_BMDM_ConsensusPeaks.bed

annotatePeaks.pl MSPC_Irf2_BMDM_ConsensusPeaks.bed   mm10 > HOMERMSPC_Irf2_BMDM_ConsensusPeaks.xls  
