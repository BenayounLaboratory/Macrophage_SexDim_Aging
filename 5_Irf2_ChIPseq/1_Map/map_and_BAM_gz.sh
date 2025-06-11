#!/bin/bash

if [[ "$#" -lt 3 ]]
then
    echo "$(basename $0) [sraDir] [oDir] "  1>&2
    echo "   [Ass]: assembly" 1>&2 
    echo "   [sraDir]: directory containing SRA files (the SRA files can be in subdirectories)" 1>&2
    echo "   [oDir]: output directory" 1>&2
    exit 1
fi

Ass=$(echo $1 | sed 's:/$::g')
sraDir=$(echo $2 | sed 's:/$::g')
oDir=$(echo $3 | sed 's:/$::g')

# make output directory if it doesnt exist
[[ ! -d "${oDir}" ]] && mkdir "${oDir}"
[[ ! -d "${oDir}/SAM" ]] && mkdir "${oDir}/SAM"
[[ ! -d "${oDir}/BAM" ]] && mkdir "${oDir}/BAM"


# run bowtie to map reads into SAM format
filePath="${sraDir}"
for f in $(find "$filePath" -name '*.fastq.gz')
do
	fileName=$(basename "${f}" | sed 's/\.fastq\.gz/\.sam/g');
 	filePath="${oDir}/SAM"
    oFname="${filePath}/${fileName}"
    gzcat $f | bowtie -p 4 -S -q -v 2 -m 3 --best --strata $BOWTIE_INDEXES/$Ass - > $oFname
done

for f in $(find "$filePath" -name '*.fq.gz')
do
	fileName=$(basename "${f}" | sed 's/\.fq\.gz/\.sam/g');
 	filePath="${oDir}/SAM"
    oFname="${filePath}/${fileName}"
    gzcat $f | bowtie -p 4 -S -q -v 2 -m 3 --best --strata $BOWTIE_INDEXES/$Ass - > $oFname
done
echo "Finished fastq -> sam conversion\n"

# convert SAM into BAM
# Macs do not parse SAM files properly when there is a space in the seq line
filePath="${oDir}/SAM"
for f in $(find "$filePath" -name '*.sam')
do
	fileName=$(basename "${f}" | sed 's/\.sam/\.bam/g');
 	filePath="${oDir}/BAM"
    oFname="${filePath}/${fileName}"
    samtools view -b -S $f > $oFname
done
echo "Finished sam -> bam conversion\n"

