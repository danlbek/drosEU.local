# drosEU.local

## Overview

This is a (currently incomplete) reformatting of the DrosEU pipeline available on [GitHub](https://github.com/capoony/DrosEU_pipeline) and described [in this paper](https://academic.oup.com/mbe/article/37/9/2661/5837682). The goal is to perform similar analysis methods in a more automated manner on WSU's Kamiak cluster. The initial wrapper implementation will be in R, calling external tools as necessary.

## Structure

The *config.R* script holds all user-definable analysis options and parameters.

The full analysis is separated into several scripts to allow for breakpoints after computationally intensive procedures. 

The *s1.sequence.processing.R* script requires the *config.R* configuration file. It reads raw FASTQ files, performs basic QC, and maps the quality reads to a reference genome using BWA. The final output will be both QC summary reports as well as sorted BAM files.

The *s2.decontamination.R* script requires the *config.R* configuration file and the output files from *s1.sequence.processing.R*. It removes any potential *D. simulans* contamination.

The *s3.snp.analysis.R* script requires the *config.R* configuration file and the output files from *s1.sequence.processing.R*. It generates the Pileup file combining all study samples and calls SNPs. It then filters SNPs based on several characteristics. Finally, the script annotates any identified SNPs.

The *s4.additional.analyses.R* script requires the *config.R* configuration file and the output files from *s1.sequence.processing.R* and *s3.snp.analysis.R*. It performs several analyses including calculation of Watterson's Theta, Tajima's D, Tajima's pi, pairwise Fst, isolation by distance, and other values.

The *s5.generate.reports.R* script requires the *config.R* configuration file and the output from the previous pipeline steps. It generates a summary of the pipeline results, including QC metrics, a description of the analyses performed with parameter values, and relevant figures and tables.

## Progress

### Mapping and QC
| Step   | Draft | Tested | Kamiak | Example | Notes |
|--------|:-----:|:------:|:------:|:-------:|:-----:|
|FastQC  |      x|       x|        |         |       |
|Cutadept|      x|       x|        |         |       |
|BWA     |      x|       x|        |         |       |
|P-sort  |      x|       x|        |         |       |
|P-filter|      x|       x|        |         |       |
|P-tag   |      x|       x|        |         |       |
|InDels  |      x| skipped|        |         |       |

### Decontamination
| Step   | Draft | Tested | Kamiak | Example | Notes |
|--------|:-----:|:------:|:------:|:-------:|:-----:|
|reformat|      x| skipped|        |         |       |
|bam2fq  |      x| skipped|        |         |       |
|BWA     |      x| skipped|        |         |       |
|fixbam  |      x| skipped|        |         |       |

### SNP calling
| Step   | Draft | Tested | Kamiak | Example | Notes |
|--------|:-----:|:------:|:------:|:-------:|:-----:|
|mpileup |      x|       x|        |         |       |
|PoolSNP |      x|        |        |         |       |
|Indels  |      x|        |        |         |       |
|RM TFF  |      x|        |        |         |       |
|annotate|      x|        |        |         |       |

### Additional analyses
| Step   | Draft | Tested | Kamiak | Example | Notes |
|--------|:-----:|:------:|:------:|:-------:|:-----:|
|pi/tht/D|      x|        |        |         |       |
|Fst     |      x|        |        |         |       |
|IBD     |      x|        |        |         |       |
|A freq  |      x|        |        |         |       |
|PCA     |       |        |        |         |       |
|pop strc|       |        |        |         |       |

### Reports
| Step   | Draft | Tested | Kamiak | Example | Notes |
|--------|:-----:|:------:|:------:|:-------:|:-----:|
|QC      |       |        |        |         |       |
|summary |       |        |        |         |       |
|results |       |        |        |         |       |



