
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
# Location of project data files (intermediate files will be stored here as well)
data.folder <- paste(project.folder, "data/", sep="")
# Location to store result files
result.folder <- paste(project.folder, "results/", sep="")
# Location to store reports and logs
report.folder <- paste(project.folder, "reports/", sep="")

# The reference.genome.name is user defined
reference.genome.name <- ""
reference.genome.location <- ""


######################
### Analysis Flags ###
######################
generate.fastqc.reports <- TRUE
trim.reads <- TRUE
map.reads <- TRUE

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

##################################
### Decontamination parameters ###
##################################

contamination.genome.name <- "D. simulans"
contamination.genome.location <- ""

