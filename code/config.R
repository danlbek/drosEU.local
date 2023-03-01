
## Main configuration script for the drosEU.local pipeline


#######################
### Project Details ###
#######################

# Give the project any name and description. These will be used in the reports. 
# It might be good practice to also use the project.name as the project folder name.
project.name <- "test"
project.description <- ""

project.folder <- "/test/"
# Location of project code files
code.folder <- paste(project.folder, "code/", sep="")
# Location of external scripts and binary files not loaded using Kamiak modules
script.folder <- paste(code.folder, "scripts/", sep="")
# Location of picard jar file (generally within the scripts folder, but could be centrally located)
picard.folder <- script.folder

# Location of project data files (intermediate files will be stored here as well)
data.folder <- paste(project.folder, "data/", sep="")
# Location of reference genomes
genome.folder <- paste(data.folder, "genome/", sep="")
# Location to store result files
result.folder <- paste(project.folder, "results/", sep="")
# Location to store reports and logs
report.folder <- paste(project.folder, "reports/", sep="")

# The reference.genome.name is user defined
reference.genome.name <- "dm6"
reference.genome.file.name <- "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"


######################
### Analysis Flags ###
######################
# Generate FastQC reports for sequencing quality verification
generate.fastqc.reports <- TRUE
# Filter and trim reads to remove adapters and low quality bases
trim.reads <- TRUE
# Map the reads to the reference genome using Bwa-mem2
map.reads <- TRUE
# Process the BAM files using Picard-tools and GTK
process.bam.to.indel.targets <- TRUE

decontaminate.libraries <- FALSE


######################
### Sample Details ###
######################
# This table links a sample name with the files associated with that sample. This may be
# better implemented as a sample spreadsheet that is imported by the pipeline.

sample.table <- data.frame(
  sampleID = c("LD7", "LD11", "LD18", "LD24", "LD29", "LD34", 
               "ZE10", "ZE14", "ZE17", "ZE22", "ZE30", "ZE31"),
  p1filename = c("18057XD-04-12_S0_L001_R1_001.fastq.gz", "18057XD-04-13_S0_L001_R1_001.fastq.gz",
                 "18057XD-04-14_S0_L001_R1_001.fastq.gz", "18057XD-04-15_S0_L001_R1_001.fastq.gz",
                 "18057XD-04-16_S0_L001_R1_001.fastq.gz", "18057XD-04-17_S0_L001_R1_001.fastq.gz", 
                 "18057XD-04-06_S0_L001_R1_001.fastq.gz", "18057XD-04-07_S0_L001_R1_001.fastq.gz",
                 "18057XD-04-08_S0_L001_R1_001.fastq.gz", "18057XD-04-09_S0_L001_R1_001.fastq.gz",
                 "18057XD-04-10_S0_L001_R1_001.fastq.gz", "18057XD-04-11_S0_L001_R1_001.fastq.gz"),
  p2filename =c("18057XD-04-12_S0_L001_R2_001.fastq.gz", "18057XD-04-13_S0_L001_R2_001.fastq.gz",
                "18057XD-04-14_S0_L001_R2_001.fastq.gz", "18057XD-04-15_S0_L001_R2_001.fastq.gz",
                "18057XD-04-16_S0_L001_R2_001.fastq.gz", "18057XD-04-17_S0_L001_R2_001.fastq.gz", 
                "18057XD-04-06_S0_L001_R2_001.fastq.gz", "18057XD-04-07_S0_L001_R2_001.fastq.gz",
                "18057XD-04-08_S0_L001_R2_001.fastq.gz", "18057XD-04-09_S0_L001_R2_001.fastq.gz",
                "18057XD-04-10_S0_L001_R2_001.fastq.gz", "18057XD-04-11_S0_L001_R2_001.fastq.gz"),
  metadata = c(rep("LD", 6), rep("ZE", 6)),
  stringsAsFactors = F
)



######################
### QC and Mapping ###
######################

## Cutadapt parameters. There are many other options and settings available through Cutadapt.
## Only the parameters used in the DrosEU pipeline are shown here, but added flexibility can
## be added as needed. Most defaults are taken from the DrosEU pipeline.

# Number of cores for parallel processing. Default is set to use all available
cutadapt.ncore <- 0
# Quality threshold for trimming low quality ends from reads
cutadapt.minq <- 18
# Filter to remove reads that are too short 
cutadapt.minlength <- 75
# Forward adapter to be trimmed
cutadapt.fwdadapter <- "ACACTCTTTCCCTACACGACGCTCTTCCGATC"
# Reverse adapter to be trimmed
cutadapt.revadapter <- "CAAGCAGAAGACGGCATACGAGAT"
# Allow for partial adapter fragments to be identified and removed by setting this minimum 
# overlap of the read and adapter.
cutadapt.minadapteroverlap <- 15
# Default assumption is at most one adapter will be present but this is modified here
cutadapt.adaptercopies <- 3

## Mapping parameters
# Number of cores for parallel processing.
bwa.ncore <- 1
# Minimum MAPQ score for mapped reads. All others will be discarded.
bwa.minq <- 20

## Java options - used by Picard-tools and GATK
# Xmx20g is a java option limiting the "heap space" (which is a bit less than the total memory used). 
# Picard documentation recommends using 2g, but the DrosEU documentations uses 20g or 10g depending 
# on the tool.
# -Dsnappy.disable=true is apparently a system property. I assume it is disabling snappy compression, 
# but I have yet to find documentation. It may not be necessary on Kamiak but is kept here since I 
# don't think it will break anything.
java.option.string <- "-Xmx10g -Dsnappy.disable=true"
##################################
### Decontamination parameters ###
##################################

contamination.genome.name <- "D. simulans"
contamination.genome.location <- ""

