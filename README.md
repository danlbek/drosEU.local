# drosEU.local

## Overview

This is a (currently incomplete) reformatting of the DrosEU pipeline available on [GitHub](https://github.com/capoony/DrosEU_pipeline) and described [in this paper](https://academic.oup.com/mbe/article/37/9/2661/5837682). The goal is to perform similar analysis methods in a more automated manner on WSU's Kamiak cluster. The initial wrapper implementation will be in R, calling external tools as necessary.

## Structure

The *config.R* script holds all user-definable analysis options and parameters.

The full analysis is separated into several scripts to allow for breakpoints after computationally intensive procedures. 

The *s1.sequence.processing.R* script requires the *config.R* configuration file. It reads raw FASTQ files, performs basic QC, and maps the quality reads to a reference genome using BWA. The final output will be both QC summary reports as well as sorted BAM files.


## Current TODO list

Outline pipeline and get working mapping step by 2/24. Run on Kamiak test data over weekend. Subsample test data as needed for speed.




