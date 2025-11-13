#Bowtie 2 version 2.4.4 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)

FqDir="/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Mph_CUT_and_RUN/H3K4me3_Aging/NGMerge"
oDir="/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Mph_CUT_and_RUN/H3K4me3_Aging/Bowtie2"

cd $oDir;

for f in $(find $FqDir -name '*_NGMerge_1.fastq.gz')
do
    f2=$(basename "${f}" | sed 's/_NGMerge_1\.fastq\.gz/_NGMerge_2\.fastq\.gz/g');
    inf2="${FqDir}/${f2}"

    of=$(basename "${f}" | sed 's/_NGMerge_1\.fastq\.gz/\.bam/g');
    oFname="${oDir}/${of}"

    bowtie2 --no-discordant --no-mixed --sensitive -p 2 -X 1000 -x $BOWTIE2_INDEXES/mm10 -1 $f -2 $inf2 | samtools view -b -S - > $oFname

done

