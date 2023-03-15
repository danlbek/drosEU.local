# DrosEU part B

## This script should be optional.

if (decontaminate.libraries){
  # This step renames all decontamination reference sequences with the prefix "sim". I'm not sure it is necessary, but
  # it is quick enough not to matter
  add.sim.command <- paste("zcat ", genome.folder, "/", contamination.genome.file.name, 
                           " | sed 's/>/>sim_/g' | gzip -c > ", genome.folder, 
                           "/dc_genome_prefix.fa.gz", sep="")
  system(add.sim.command)
  
  
  # Combine decontamination reference sequences with study reference
  combine.references.command <- paste("zcat ", genome.folder, "/dc_genome_prefix.fa.gz | cat ",
                                      genome.folder, "/", reference.genome.file.name,
                                      " - | gzip -c > ", genome.folder, "combined.fa.gz", sep="")
  system(combine.references.command)
  
  # Create BWA-MEM index for combined reference
  bwa.index.comb.command <- paste("bwa index ", genome.folder, "/combined.fa.gz", sep="")
  system(bwa.index.comb.command)
  
  for (i in 1:nrow(sample.table)) {
    
    # The DrosEU pipeline uses bam2fastq here. That script is not maintained and can be replaced by
    # the SamToFastq tool in the Picard toolset. I don't see a good reason to use the original.
  
    samtofastq.command <- paste("java ", java.option.string, " -jar ", picard.folder, "/picard.jar SamToFastq ",
                                "I=", normalizePath(data.folder), "/", "dedup.rg.", sample.table$sampleID[i], ".bam ", 
                                "FASTQ=", sample.table$sampleID[i], ".R1.fromBAM.fq ",
                                "SECOND_END_FASTQ=", sample.table$sampleID[i], ".R2.fromBAM.fq ",
                                "UNPAIRED_FASTQ=", sample.table$sampleID[i], ".RNP.fromBAM.fq", sep="")
    system(samtofastq.command)
    
    # Map reads from BAM file to combined reference genome. This seems out of order. Why not map reads
    # to combined reference initially, then remove all reads mapping to the decontamination genome?
    # This should be further explored, since the mapping is computationally expensive and this step
    # seems to nearly double the overall run time.
    # Other problems, the current implementation ignores non-paired reads. The output is a SAM file. I've
    # left the file version alone, since the next step may require it. Unused intermediate files should be
    # identified and deleted.
    bwa.comb.command <- paste("bwa mem -M -t ", bwa.ncore, " ", genome.folder, "/combined.fa.gz",
                              " ", data.folder, "/", sample.table$sampleID[i], ".R1.fromBAM.fq",
                              " ", data.folder, "/", sample.table$sampleID[i], ".R2.fromBAM.fq",
                              " > ", data.folder, "/", sample.table$sampleID[i], ".deSim.sam", sep="")
    system(bwa.comb.command)
    
    # Deconvolute the reads in the original BAM file. Possible alternatives should be explored. There may be
    # better alternatives in SAMTools, Picard, or GATK. At the very least it should be updated to Python3.
    
    fix.bam.command <- paste("python2.7 ", script.folder, "/FixBAM.py --contaminated ",
                             normalizePath(data.folder), "/", "dedup.rg.", sample.table$sampleID[i], ".bam ",
                             "--prefix sim_ --detect ", data.folder, "/", sample.table$sampleID[i], ".deSim.sam ",
                             "--output ", data.folder, "/", sample.table$sampleID[i], ".decon ", sep="")
    system(fix.bam.command)
  }
}

