#!/bin/bash

if [[ "$#" -lt 2 ]]
then
    echo "$(basename $0) [Ass] [BamDir] [oDir] "  1>&2
    echo "   [Ass]: assembly" 1>&2 
    echo "   [BamDir]: folder with bam files" 1>&2
    echo "   [oDir]: output directory" 1>&2
    exit 1
fi

Ass=$(echo $1 | sed 's:/$::g')
BamDir=$(echo $2 | sed 's:/$::g')
oDir=$(echo $3 | sed 's:/$::g')

# make output directory if it doesnt exist
[[ ! -d "${oDir}" ]] && mkdir "${oDir}"
[[ ! -d "${oDir}/BAM_clean" ]] && mkdir "${oDir}/BAM_clean"

##########################################
# 1. convert SAM into BAM, clean up alignemnts
for f in $(find "${BamDir}" -name '*.bam')
do
	filePath="${oDir}/BAM_clean"
	fileName=$(basename "${f}" | sed 's/\.bam/\.mp\.bam/g');
    oFname="${filePath}/${fileName}"
    samtools view -b -F 4 -q 15 $f > $oFname # mapped only
    
	outBamPre=$(basename "${f}" | sed 's/\.bam/\.mp\.srt/g');
	oFname2="${filePath}/${outBamPre}"
	samtools sort $oFname $oFname2
	
	outBam2=$(basename "${f}" | sed 's/\.bam/\.mp\.srt\.bam/g');
	oFname2="${filePath}/${outBam2}"
	
	outBam3=$(basename "${f}" | sed 's/\.bam/\.mp\.srt\.rmdup\.bam/g');
	oFname3="${filePath}/${outBam3}"
	samtools rmdup -s $oFname2 $oFname3
	
	# determine duplication rates and store in report file
	getDupRate_PE.pl $oFname $oFname3 >> $filePath/Libraries_PCR_duplication_rates_report.txt


	# make HOMER files
	outHOM=$(basename "${f}" | sed 's/\.bam/_CLEAN_TAGs/g');
	oDir4="${filePath}/${outHOM}"
	makeTagDirectory $oDir4 $oFname3 -genome $Ass -format sam -keepOne -normGC default -single
    
done
echo "Finished bam clean up"
