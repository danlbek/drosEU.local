# DrosEU part A includes QC, read mapping, and re-mapping. The output for this step should be 
# QC reports and summaries, as well as a sorted BAM file for each sample.

## Load configuration script and any required libraries
source("config.R")

## Print versions for all used applications to log file

############################
## Generate FastQC report ##
############################
## This section generates the FastQC reports for each sample file. Separate reports are 
## generated for each paired end. Should be combined into single MultiQC report.

if (generate.fastqc.reports) {
  fastqc.command <- paste("fastqc -o", report.folder, 
                          paste(paste(data.folder, c(sample.table$p1filename, 
                                                     sample.table$p2filename), sep=""), 
                                collapse=" "), sep=" ")
  system(fastqc.command)
}

#######################
## Trim/filter reads ##
#######################
## This section runs Cutadapt to trim and filter low quality reads.

if (trim.reads) {
  for (i in 1:nrow(sample.table)) {
    # Build system command
    cutadapt.command <- paste("cutadapt --cores ", cutadapt.ncore, 
                              " --quality-cutoff ", cutadapt.minq, 
                              " --minimum-length ", cutadapt.minlength,
                              " -b ", cutadapt.fwdadapter,
                              " -B ", cutadapt.revadapter,
                              " --overlap ", cutadapt.minadapteroverlap,
                              " --times ", cutadapt.adaptercopies,
                              " -o ", data.folder, "/trimmed.", sample.table$p1filename[i],
                              " -p ", data.folder, "/trimmed.", sample.table$p2filename[i],
                              " --json ", report.folder, "/log.", sample.table$sampleID[i], ".cutadapt.json ",
                              data.folder, sample.table$p1filename[i], " ",
                              data.folder, sample.table$p2filename[i], sep="")
    # Run command
    system(cutadapt.command)
                              
  }
}
    
###############
## Map reads ##
###############
# Map reads to genome using BWA. SAM file is passed to Samtools without saving to disk.
# Samtools filters out low quality reads and saves the file in BAM format.

# Bwa-mem2 may be a better option. It claims to be substantially faster.
# The 0x100 flag removes secondary alignments
if (map.reads) {
  if (create.bwa.index){
    bwa.index.command <- paste("bwa index ", genome.folder, "/", reference.genome.file.name, sep="")
    system(bwa.index.command)
  }
  for (i in 1:nrow(sample.table)) {
    bwa.command <- paste("bwa mem -M -t ", bwa.ncore, " ", genome.folder, "/", reference.genome.file.name,
                         " ", data.folder, "/trimmed.", sample.table$p1filename[i],
                         " ", data.folder, "/trimmed.", sample.table$p2filename[i],
                         " | samtools view -Sbh --threads ", bwa.ncore-1, " -q ", bwa.minq, " -F 0x100 -o ", 
                         data.folder, sample.table$sampleID[i], ".bam", sep="")
                         
    system(bwa.command)
  }
}

#########################
## BAM file processing ##
#########################
# Using Picard, sort BAM files, remove PCR duplicates and add group tags. DrosEU uses old Picard 
# version. I've updated it here but should watch for problems. The Kamiak version of Picard is an
# intermediate version.
# This stores all intermediate BAM files. This should be modified if the intermediate files are not 
# used downstream.

if (process.bam.to.indel.targets) {
  for (i in 1:nrow(sample.table)) {
    
    # Sort BAM file. 
    sortsam.command <- paste("java ", java.option.string, " -jar ", picard.folder, "/picard.jar SortSam ",
                             "I=", normalizePath(data.folder), "/", sample.table$sampleID[i], ".bam ", 
                             "O=", normalizePath(data.folder), "/", "sorted.", sample.table$sampleID[i], ".bam ", 
                             "SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT", sep="")
    system(sortsam.command)

    # Remove PCR duplicates
    dedup.command <- paste("java ", java.option.string, " -jar ", picard.folder, "/picard.jar MarkDuplicates ",
                           "I=", normalizePath(data.folder), "/", "sorted.", sample.table$sampleID[i], ".bam ", 
                           "O=", normalizePath(data.folder), "/", "dedup.", sample.table$sampleID[i], ".bam ", 
                           "M=", normalizePath(report.folder), "/", "dedup.metrics.", sample.table$sampleID[i], ".bam ",
                           "VALIDATION_STRINGENCY=SILENT", sep="")
    system(dedup.command)
    
    # Add read group tags. RGID (ID), RGLB (library), RGSM (sample), and RGPU (barcode) are all set to sample name. 
    rgtag.command <- paste("java ", java.option.string, " -jar ", picard.folder, 
                           "/picard.jar AddOrReplaceReadGroups ",
                           "I=", normalizePath(data.folder), "/", "dedup.", sample.table$sampleID[i], ".bam ", 
                           "O=", normalizePath(data.folder), "/", "dedup.rg.", sample.table$sampleID[i], ".bam ", 
                           "SORT_ORDER=coordinate",
                           " RGID=", sample.table$sampleID[i], 
                           " RGLB=", sample.table$sampleID[i], 
                           " RGPL=illumina",
                           " RGSM=", sample.table$sampleID[i],
                           " RGPU=", sample.table$sampleID[i], 
                           " CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT", sep="")
    system(rgtag.command)
    
    # # Generate InDel target list and re=align those positions using GATK3
    # # These functions are eliminated in GATK4. GATK3 is no longer supported. Does it make sense to keep 
    # # this analysis in place? Skipping for now.
    # gatk3.target.command <- paste("java ", java.option.string, " -jar ", gatk.folder, 
    #                               "/GenomeAnalysisTK.jar -T RealignerTargetCreator",
    #                               " -R ", genome.folder, reference.genome.file.name,
    #                               " -I ", data.folder, "dedup.rg.", sample.table$sampleID[i], ".bam",
    #                               " -o ", data.folder, "dedup.rg.", sample.table$sampleID[i], ".list",
    #                               sep="")
    # system(gatk3.target.command)
    # 
    # gatk3.realign.command <- paste("java ", java.option.string, " -jar ", gatk.folder, 
    #                                "/GenomeAnalysisTK.jar -T IndelRealigner",
    #                                "-R ", genome.folder, reference.genome.file.name,
    #                                " -I ", normalizePath(data.folder), "/", "dedup.rg.", sample.table$sampleID[i], ".bam",
    #                                " -targetIntervals ", normalizePath(data.folder), "/", "dedup.rg.", sample.table$sampleID[i], ".list",
    #                                " -o ", normalizePath(data.folder), "/", "dedup.rg.indel.", sample.table$sampleID[i], ".bam",
    #                                sep="")
    # system(gatk3.realign.command)
    
  }
}



