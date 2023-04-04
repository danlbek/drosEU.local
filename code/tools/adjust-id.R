## This script adjusts the TE libraries used in s3.snp.analysis.R.
## It is intended to be run manually, since there are several hard-coded parameters
## that will need to be modified when the TE libraries change. This script replaces
## the missing adjust-id.py script from the DrosEU pipeline.

## The TE libraries are downloaded using curl
#curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.10_FB2016_02/fasta/dmel-all-transposon-r6.10.fasta.gz
#curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.10_FB2016_02/fasta/dmel-all-chromosome-r6.10.fasta.gz

library(Biostrings)

## I'm not sure if the chromosome file needs to be modified. It may remain unused.
# Read fasta file
all.chr <- readDNAStringSet(filepath="../data/genome/dmel-all-chromosome-r6.10.fasta.gz")
# Replace names with only the chromosome name
names(all.chr) <- sapply(strsplit(names(all.chr), split=" "), function(i) i[1])
# Write to fasta
writeXStringSet(x=all.chr, filepath="../data/genome/dmel-all-chromosome-r6.10.fixed-id.fasta.gz",
                compress=TRUE, format="fasta")

# Do the same thing for the transposon file
all.tran <- readDNAStringSet(filepath="../data/genome/dmel-all-transposon-r6.10.fasta.gz")
names(all.tran) <- sapply(strsplit(names(all.tran), split=" "), function(i) i[1])
writeXStringSet(x=all.tran, filepath="../data/genome/dmel-all-transposon-r6.10.fixed-id.fasta.gz",
                compress=TRUE, format="fasta")

