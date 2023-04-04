# DrosEU part C merges all BAM files in a MPILEUP file and calls SNPs

## Load configuration script and any required libraries
source("config.R")

#################################
## Create combined pileup file ##
#################################

# Samtools is used to combine the BAM files into a Pileup file. The size of the Pileup file may be concerning.
bamfilelist <- paste(normalizePath(data.folder), "/", "dedup.rg.", sample.table$sampleID, ".bam", sep="")
writeLines(bamfilelist, con = paste(data.folder, "/bamfilelist.txt", sep=""))
mpileup.command <- paste("samtools mpileup -B -f ", genome.folder, "/", reference.genome.file.name,
                         " -b ", data.folder, "/bamfilelist.txt -q ",pileup.min.mapping.quality, 
                         " -Q ", pileup.min.base.quality, " | gzip > ", data.folder, "/all.pileup.gz", 
                         sep="")
system(mpileup.command)
  
###############
## Call SNPs ##
###############

poolsnp.command <- paste("bash ", script.folder, "/PoolSNP/PoolSNP.sh mpileup=all.pileup.gz",
                         " reference=", genome.folder, "/", reference.genome.file.name, 
                         " names=", paste(sample.table$sampleID, collapse=","),
                         " max-cov=", poolsnp.maxcov, 
                         " min-cov=", poolsnp.mincov,
                         " min-count=", poolsnp.mincount,
                         " min-freq=", poolsnp.minfreq,
                         " miss-frac=", poolsnp.missfrac,
                         " jobs=", poolsnp.jobs,
                         " output=all.snps", sep="")
system(poolsnp.command)

###############################
## Mask InDel adjacent sites ##
###############################
detectindel.command <- paste("bash ", script.folder, "/DetectIndels.py --mpileup all.pileup.gz",
                         " --minimum-count ", detectindel.minimum.count,
                         " --mask ", detectindel.mask, 
                         " | gzip > Indel-positions.txt.gz", sep="")
system(detectindel.command)

## The DrosEU pipeline downloads and partially modifies TE libraries here.
## They use a python script named "adjust-id.py". I do not see that script included
## in the repository. I've downloaded the TE libraries and modified them using the
## included "adjust-id.R" script based on the description of the missing python script.
## The next steps assume the TE libraries are modified and present in the genome folder.

#####################################
## Mask repeats using RepeatMasker ##
#####################################
repeatmasker.command <- paste(script.folder, "/RepeatMasker/RepeatMasker",
                              " -pa ", repeatmasker.pa,
                             " --lib ", genome.folder, "/dmel-all-transposon-r6.10.fixed-id.fasta.gz",
                             " --gff --qq --no_is --nolow",
                             " dmel.all.chromosome.r6.10.fasta", sep="")
system(repeatmasker.command)

#######################################
## Filter SNPs around InDels and TEs ##
#######################################
## This script uses Python 2.7. This should be updated to python3.
########### Fix file names during testing. Paths and exact names are unclear ############
filter.snps.command <- paste("python2.7 ", script.folder, "/FilterPosFromVCF.py",
                              " --indel Indel-positions.txt.gz",
                              " --te dmel-all-chromosome-r6.10.fasta.out.gff",
                              " --vcf SNPs.vcf.gz", 
                              " | gzip > SNPs.clean.vcf.gz", sep="")
system(filter.snps.command)

###############################
## Annotate SNPs with snpEff ##
###############################

annotate.snps.command <- paste("java ", java.option.string, " -jar ", script.folder, "/snpEff/snpEff.jar",
                               " -ud ", snpeff.ud,
                               " BDGP6.82 -stats SNPs_clean.html SNPs.clean.vcf.gz",
                               " | gzip > SNPs.clean-ann.vcf.gz", sep="")
                               
system(annotate.snps.command)





