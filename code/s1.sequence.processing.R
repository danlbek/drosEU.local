# DrosEU part A includes QC, read mapping, and re-mapping. The output for this step should be QC reports and summaries, as well as a sorted BAM file for each sample.


## Print versions for all used applications to log file

## Insert fastQC or equivalent here

# Trim raw FASTQ files using cutadapt

    # export PATH=$PATH:scripts/cutadapt-1.8.3/bin
    # cutadapt -q 18 --minimum-length 75 -o trimmed-read1.fq.gz -p trimmed-read2.fq.gz -b ACACTCTTTCCCTACACGACGCTCTTCCGATC -B CAAGCAGAAGACGGCATACGAGAT -O 15 -n 3 read1.fq.gz read2.fq.gz
    
# Map reads to genome using BWA and filter unwanted reads using samtools

    # export PATH=$PATH:scripts/samtools-0.1.19
    # export PATH=$PATH:scripts/bwa-0.7.15
    # bwa mem -M -t 24 reference.fa.gz trimmed-read1.fq.gz trimmed-read2.fq.gz | samtools view -Sbh -q 20 -F 0x100 - > library.bam

# Using Picard, sort BAM files, remove PCR duplicates and add group tags

    # java -Xmx20g -Dsnappy.disable=true -jar scripts/picard-tools-1.109/SortSam.jar I=library.bam O=library-sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT
    # 
    # java -Xmx20g -Dsnappy.disable=true -jar scripts/picard-tools-1.109/MarkDuplicates.jar REMOVE_DUPLICATES=true I=library-sort.bam O=library-dedup.bam M=library-dedup.txt VALIDATION_STRINGENCY=SILENT
    # 
    # java -jar -Xmx10g scripts/picard-tools-1.109/AddOrReplaceReadGroups.jar INPUT=librtary-dedup.bam OUTPUT=library-dedup_rg.bam SORT_ORDER=coordinate RGID=library RGLB=library RGPL=illumina RGSM=sample RGPU=library CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

# Generate InDel target list and re=align those positions using GATK

    # java -Xmx20g -jar scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T IndelRealigner -R reference.fa -I library-dedup_rg.bam -targetIntervals library-dedup_rg.list -o library-dedup_rg_InDel.bam




