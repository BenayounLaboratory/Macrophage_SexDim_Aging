mspc -f ./Irf2_BMDM/   -c 2 -r bio -w 1e-5 -s 1e-8 -m Lowest --excludeHeader -o MSPC_Irf2_BMDM

cp ./MSPC_Irf2_BMDM/ConsensusPeaks.bed     MSPC_Irf2_BMDM_ConsensusPeaks.bed  

wc -l *_ConsensusPeaks.bed
#     339 MSPC_Irf2_BMDM_ConsensusPeaks.bed

annotatePeaks.pl MSPC_Irf2_BMDM_ConsensusPeaks.bed   mm10 > HOMER_MSPC_Irf2_BMDM_ConsensusPeaks.xls  
