
## Main configuration script for the drosEU.local pipeline


#######################
### Project Details ###
#######################

# Give the project any name and description. These will be used in the reports. 
# It might be good practice to also use the project.name as the project folder name.
project.name <- "test"
project.description <- "Small synthetic test set"

project.folder <- "~/test/"
# Location of project code files
code.folder <- paste(project.folder, "code/", sep="")
# Location of external scripts and binary files not loaded using Kamiak modules
script.folder <- paste(code.folder, "tools/", sep="")
## Both Picard and GATK are available as Kamiak modules. I'm using a local version here since
## the DrosEU pipeline requires obselete functions only available in GATK 3.8.
# Location of picard jar file
picard.folder <- paste(script.folder, "picard/", sep="")
# Location of gatk jar file
gatk.folder <- paste(script.folder, "gatk3.8/", sep="")

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
generate.fastqc.reports <- FALSE
# Filter and trim reads to remove adapters and low quality bases
trim.reads <- FALSE
# Index reference genome (required for BWA, but generated index can be copied to other projects)
create.bwa.index <- FALSE
# Map the reads to the reference genome using Bwa-mem2
map.reads <- FALSE
# Process the BAM files using Picard-tools and GTK
process.bam.to.indel.targets <- TRUE

decontaminate.libraries <- FALSE


######################
### Sample Details ###
######################
# This table links a sample name with the files associated with that sample. This may be
# better implemented as a sample spreadsheet that is imported by the pipeline.

# sample.table <- data.frame(
#   sampleID = c("LD7", "LD11", "LD18", "LD24", "LD29", "LD34", 
#                "ZE10", "ZE14", "ZE17", "ZE22", "ZE30", "ZE31"),
#   p1filename = c("18057XD-04-12_S0_L001_R1_001.fastq.gz", "18057XD-04-13_S0_L001_R1_001.fastq.gz",
#                  "18057XD-04-14_S0_L001_R1_001.fastq.gz", "18057XD-04-15_S0_L001_R1_001.fastq.gz",
#                  "18057XD-04-16_S0_L001_R1_001.fastq.gz", "18057XD-04-17_S0_L001_R1_001.fastq.gz", 
#                  "18057XD-04-06_S0_L001_R1_001.fastq.gz", "18057XD-04-07_S0_L001_R1_001.fastq.gz",
#                  "18057XD-04-08_S0_L001_R1_001.fastq.gz", "18057XD-04-09_S0_L001_R1_001.fastq.gz",
#                  "18057XD-04-10_S0_L001_R1_001.fastq.gz", "18057XD-04-11_S0_L001_R1_001.fastq.gz"),
#   p2filename =c("18057XD-04-12_S0_L001_R2_001.fastq.gz", "18057XD-04-13_S0_L001_R2_001.fastq.gz",
#                 "18057XD-04-14_S0_L001_R2_001.fastq.gz", "18057XD-04-15_S0_L001_R2_001.fastq.gz",
#                 "18057XD-04-16_S0_L001_R2_001.fastq.gz", "18057XD-04-17_S0_L001_R2_001.fastq.gz", 
#                 "18057XD-04-06_S0_L001_R2_001.fastq.gz", "18057XD-04-07_S0_L001_R2_001.fastq.gz",
#                 "18057XD-04-08_S0_L001_R2_001.fastq.gz", "18057XD-04-09_S0_L001_R2_001.fastq.gz",
#                 "18057XD-04-10_S0_L001_R2_001.fastq.gz", "18057XD-04-11_S0_L001_R2_001.fastq.gz"),
#   metadata = c(rep("LD", 6), rep("ZE", 6)),
#   stringsAsFactors = F
# )

sample.table <- data.frame(
  sampleID = c("g1.rss1", "g1.rss2", "g1.rss3", "g2.rss1", "g2.rss2", "g2.rss3"),
  p1filename = c("test.g1.rss1.R1.fastq.gz", "test.g1.rss2.R1.fastq.gz", "test.g1.rss3.R1.fastq.gz",
                 "test.g2.rss1.R1.fastq.gz", "test.g2.rss2.R1.fastq.gz", "test.g2.rss3.R1.fastq.gz"),
  p2filename =c("test.g1.rss1.R2.fastq.gz", "test.g1.rss2.R2.fastq.gz", "test.g1.rss3.R2.fastq.gz",
                "test.g2.rss1.R2.fastq.gz", "test.g2.rss2.R2.fastq.gz", "test.g2.rss3.R2.fastq.gz"),
  metadata = c(rep("g1", 3), rep("g2", 3)),
  stringsAsFactors = F
)


######################
### QC and Mapping ###
######################

## Cutadapt parameters. There are many other options and settings available through Cutadapt.
## Only the parameters used in the DrosEU pipeline are shown here, but flexibility can be 
## added as needed. Most defaults are taken from the DrosEU pipeline.

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
bwa.ncore <- 4
# Minimum MAPQ score for mapped reads. All others will be discarded.
bwa.minq <- 20

## Java options - used by Picard-tools and GATK
# Xmx20g is a java option limiting the "heap space" (which is a bit less than the total memory used). 
# Picard documentation recommends using 2g, but the DrosEU documentations uses 20g or 10g depending 
# on the tool.
# -Dsnappy.disable=true is apparently a system property. I assume it is disabling snappy compression, 
# but I have yet to find documentation. It may not be necessary on Kamiak. I've removed it here.
java.option.string <- "-Xmx10g"

##################################
### Decontamination parameters ###
##################################

contamination.genome.name <- "Drosophila simulans"
contamination.genome.file.name <- "GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz"

##############################
### SNP calling parameters ###
##############################

## Pileup file generation
pileup.min.mapping.quality <- 20
pileup.min.base.quality <- 20

## PoolSNP.sh SNP calling
poolsnp.names <- c()
poolsnp.maxcov <- 0.99
poolsnp.mincov <- 10
poolsnp.mincount <- 10
poolsnp.minfreq <- 0.001
poolsnp.missfrac <- 0.2
poolsnp.jobs <- 4 # This should be combined with other available thread parameters by default

## DetectIndels.py InDel masking
detectindel.minimum.count <- 20
detectindel.mask <- 5

## RepeatMasker
repeatmasker.pa <- 20

## snpEFF
snpeff.ud <- 2000

## SubsampleSync.py
subsamplesync.targetcov <- 40
subsamplesync.mincov <- 10

## TrueWindows.py
truewindows.window <- 200000
truewindows.step <- 200000

## PoolGen_var.py
poolgenvar.poolsize <- c(rep(80,35), 66, rep(80,8), 70, rep(80,3))
poolgenvar.mincount <- 2
poolgenvar.window <- 200000
poolgenvar.step <- 200000
poolgenvar.minsitesfrac <- 0.75

## IntronicSnps.py
intronicsnps.targetlength <- 60

## FilterByRecomRateNInversion.py
filterinversions.D <- 1000000 
filterinversions.r <- 3 

## FST.py
fst.poolsize <- poolgenvar.poolsize ## Both these should be specified in a more general way.
fst.minimumcount <- 2 
fst.minimumcov <- 10 

## CombineFST.py
combinefst.stat <- 0







