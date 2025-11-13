FqDir="/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq/NGMerge"
BowtieDir="/Volumes/BB_Home_HQ/Immune_sex_dimorphism_Aging/Macrophages/Clean_Macrophage_Data/ATACseq/Bowtie2"


for f in $(find $FqDir -name '*_NGmerge_1.fastq.gz')
do
    f2=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/_NGmerge_2\.fastq\.gz/g');
    inf2="${FqDir}/${f2}"

    of=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/\.bam/g');
    oFname="${BowtieDir}/${of}"

    bowtie2 --sensitive --no-discordant --no-mixed  -k 2 -p 2 -X 2500 -x $BOWTIE2_INDEXES/mm10 -1 $f -2 $inf2 | samtools view -b -S - > $oFname

done

