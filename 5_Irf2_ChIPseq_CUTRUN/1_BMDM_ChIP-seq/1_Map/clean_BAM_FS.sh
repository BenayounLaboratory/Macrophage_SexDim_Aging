#!/bin/bash

if [[ "$#" -lt 3 ]]
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
	
	# filter to retain mapped reads
	fileName=$(basename "${f}" | sed 's/\.bam/\.mp\.bam/g');
    oFname="${filePath}/${fileName}"
    samtools view -b -F 4 -q 15 $f > $oFname # mapped only
    
    # sort mapped reads (space gain)
	outBamPre=$(basename "${f}" | sed 's/\.bam/\.mp\.srt/g');
	oFname2b="${filePath}/${outBamPre}"
	outBam2=$(basename "${f}" | sed 's/\.bam/\.mp\.srt\.bam/g');
	oFname2="${filePath}/${outBam2}"

	#echo $oFname2b
	#echo $oFname2
	
	samtools sort $oFname $oFname2b
	
	# clean up pre-sort file
	rm $oFname
	
	# index sorted bam file
	samtools index $oFname2

	# Create the read count histogram for Fixseq input. 
	# These flags exclude unmapped reads, secondary alignments,
	# and QC failures.
	outFile3=$(basename "${f}" | sed 's/\.bam/\.INPUT\.csv/g');
	oFname3="${filePath}/${outFile3}"
	samtools view -F 772 $oFname2 | cut -f 3,4 | uniq -c | awk '{print $1}' | sort | uniq -c | sort -n | awk '{print $1","$2}' > $oFname3
	
	# Add the number of bases with no reads. 
	# To use a 'mappable fraction of base' estimate change the print int(x*1 - y)
	# to print(x*fraction - y)
	samtools view -H $oFname2 | ggrep -oP "LN:\d+" | awk -v y=`awk -F , '{x+=$2} END {print x}' $oFname3` -F : '{x+=$2} END {print int(x*1 - y)",0"}' >> $oFname3
	
	# Next, run modified methods.r and obtain corresponding output.csv file.
	FIXSEQ_cmdLine.r $oFname3
	
	# From FIXSEQ help: Choose a rounding scheme you'd like to employ, and create a
	# remapping count file. Here, use the default
	cut -f 2,6 -d "," output.csv | tr "," "\t" | tail -n +2 | sort -n > $filePath/counts_map.txt
	
	outFile4=$(basename "${f}" | sed 's/\.bam/\.FIXSEQ_OUTPUT\.csv/g');
	oFname4="${filePath}/${outFile4}"
	
	mv output.csv $oFname4
	
	# Now use BEDTools and FixSEQ Python script to remap read counts
	# by expanding/contracting a bed file.  IMPORTANT: The original bam
	# must be coordinate-sorted, since the remapper expects duplicate
	# reads to be adjacent to each other in the bam (then bed).
	outFile5=$(basename "${f}" | sed 's/\.bam/\.FIXSEQ_CLEANED_READS\.bed/g');
	oFname5="${filePath}/${outFile5}"

	bamToBed -i $oFname2  | python /Users/benayoun/Softwares/fixseq-e2cc14f19764/bed_count_remapper.py $filePath/counts_map.txt > $oFname5
	
	#rm $oFname3
	
	
	# determine duplication rates and store in report file
	getDupRate.pl $oFname2 $oFname5 >> $filePath/Libraries_PCR_duplication_rates_report.txt

	
	outHOM=$(basename "${f}" | sed 's/\.bam/_FIXSEQ_TAGs/g');
	oDir4="${filePath}/${outHOM}"
	makeTagDirectory $oDir4 $oFname5 -genome $Ass -format bed
	
    
done
echo "Finished bam clean up"


# Clean up temporary file
rm $filePath/counts_map.txt
