FqDir="/Users/berenice/Desktop/Irf2_CUTRUN/FASTQ"
oDir="/Users/berenice/Desktop/Irf2_CUTRUN/Bowtie2"

cd $oDir;

for f in $(find $FqDir -name '*r1.fastq.gz')
do
    f2=$(basename "${f}" | sed 's/_r1\.fastq\.gz/_r2\.fastq\.gz/g');
    inf2="${FqDir}/${f2}"

	# Map to mouse genome (CUT&RUN core)
    of1=$(basename "${f}" | sed 's/_r1\.fastq\.gz/\.mm10\.bam/g');
    oFname1="${oDir}/${of1}"

    bowtie2 --no-discordant --no-mixed --sensitive -p 2 -X 1000 -x $BOWTIE2_INDEXES/mm10 -1 $f -2 $inf2 | samtools view -b -S - > $oFname1

	# Map to drosophila genome (CUT&RUN spike ins)
    of2=$(basename "${f}" | sed 's/_r1\.fastq\.gz/\.dm6\.bam/g');
    oFname2="${oDir}/${of2}"
    
    bowtie2 --no-discordant --no-mixed --very-fast -p 2 -X 1000 -x $BOWTIE2_INDEXES/dm6 -1 $f -2 $inf2 | samtools view -b -S - > $oFname2

done

###### find ../FASTQ/ -name '*r1.fastq.gz'
#### ../FASTQ//04_0LRJ_02ZNUSC_OF-Pool2_IRF2_mm_i4-13_r1.fastq.gz
#### ../FASTQ//02_0LRH_02ZNUSC_YF-Pool2_IRF2_mm_i2-13_r1.fastq.gz
#### ../FASTQ//01_0LRG_02ZNUSC_YF-Pool1_IRF2_mm_i1-13_r1.fastq.gz
#### ../FASTQ//03_0LRI_02ZNUSC_OF-Pool1_IRF2_mm_i3-13_r1.fastq.gz

### 28177840 reads; of these:
###   28177840 (100.00%) were paired; of these:
###     6990332 (24.81%) aligned concordantly 0 times
###     15873766 (56.33%) aligned concordantly exactly 1 time
###     5313742 (18.86%) aligned concordantly >1 times
### 75.19% overall alignment rate
### 28177840 reads; of these:
###   28177840 (100.00%) were paired; of these:
###     25641718 (91.00%) aligned concordantly 0 times
###     2095630 (7.44%) aligned concordantly exactly 1 time
###     440492 (1.56%) aligned concordantly >1 times
### 9.00% overall alignment rate
### 27634551 reads; of these:
###   27634551 (100.00%) were paired; of these:
###     3880005 (14.04%) aligned concordantly 0 times
###     18022838 (65.22%) aligned concordantly exactly 1 time
###     5731708 (20.74%) aligned concordantly >1 times
### 85.96% overall alignment rate
### 27634551 reads; of these:
###   27634551 (100.00%) were paired; of these:
###     26727354 (96.72%) aligned concordantly 0 times
###     741404 (2.68%) aligned concordantly exactly 1 time
###     165793 (0.60%) aligned concordantly >1 times
### 3.28% overall alignment rate
### 28176693 reads; of these:
###   28176693 (100.00%) were paired; of these:
###     5573105 (19.78%) aligned concordantly 0 times
###     17212798 (61.09%) aligned concordantly exactly 1 time
###     5390790 (19.13%) aligned concordantly >1 times
### 80.22% overall alignment rate
### 28176693 reads; of these:
###   28176693 (100.00%) were paired; of these:
###     25765520 (91.44%) aligned concordantly 0 times
###     1994088 (7.08%) aligned concordantly exactly 1 time
###     417085 (1.48%) aligned concordantly >1 times
### 8.56% overall alignment rate
### 29190477 reads; of these:
###   29190477 (100.00%) were paired; of these:
###     6993923 (23.96%) aligned concordantly 0 times
###     16480333 (56.46%) aligned concordantly exactly 1 time
###     5716221 (19.58%) aligned concordantly >1 times
### 76.04% overall alignment rate
### 29190477 reads; of these:
###   29190477 (100.00%) were paired; of these:
###     25967711 (88.96%) aligned concordantly 0 times
###     2682418 (9.19%) aligned concordantly exactly 1 time
###     540348 (1.85%) aligned concordantly >1 times
### 11.04% overall alignment rate