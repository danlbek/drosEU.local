# DrosEU part B


## This script should be optional. It looks reasonably quick, so may be combined with s3.snp.analysis.R in the future

# # add "sim_" to FASTA headers
# sed 's/>/>sim_/g' reference/sim_genome.fa | gzip -c > reference/sim_genome_prefix.fa.gz
# 
# # combine with D. melanogaster reference
# zcat reference/sim_genome_prefix.fa.gz | cat reference.fa - | gzip -c > reference/combined.fa.gz
# 
# # extract high confidence mapped reads from BAM file with bam2fastq
# 
# export PATH=$PATH:scripts/bam2fastq-1.1.0
# bam2fastq -s -o reads/library# library-dedup_rg_InDel.bam
# 
# # competitive mapping of extracted reads against combined reference genomes
# export PATH=$PATH:scripts/bwa-0.7.15
# bwa mem -Mt 20 reference/combined.fa.gz reads/library\_1.gz reads/library\_2.gz > library_deSim.sam
# 
# # deconvolute the reads in the original BAM file
# python2.7 scripts/FixBAM.py --contaminated library-dedup_rg_InDel.bam --prefix sim_ --detect library_deSim.sam --output library_deSim