# DrosEU part D calculates Tajima's pi, Watterson's Theta and Tajima's D
# DrosEU part E infers demographic patterns

## Load configuration script and any required libraries
source("config.R")

## DrosEU part D starts here


####################################
## Change format from VCF to SYNC ##
#################################### 

vcf2sync.command <- paste("python ", script.folder, "/VCF2sync.py",
                          " --vcf SNPs.clean-ann.vcf.gz",
                          " | gzip > SNPs.sync.gz", sep="")
system(vcf2sync.command)

####################
## Re-sample SNPS ##
####################

subsamplesync.command <- paste("python ", script.folder, "/SubsampleSync.py",
                               " --sync SNPs.sync.gz",
                               " --target-cov ", subsamplesync.targetcov,
                               " --min-cov ", subsamplesync.mincov,
                               " | gzip > SNPs.40x.sync.gz", sep="")
system(subsamplesync.command)

###################################
## Calculate "true" window-sizes ##
###################################

truewindows.command <- paste("python ", script.folder, "/TrueWindows.py",
                             " --badcov SNP_BS.txt.gz",
                             " --indel Indel-positions.txt.gz ",
                             " --te te.gff ",
                             " --window ", truewindows.window, 
                             " --step ", truewindows.step,
                             " --output truewindows", sep="")
system(truewindows.command)

###################################################
## Calculate Tajima's pi/D and Watterson's Theta ##
###################################################

PoolGen_var.command <- paste("python ", script.folder, "/PoolGen_var.py",
                             " --input SNPs.40x.sync.gz",
                             " --pool-size ", paste(poolgenvar.poolsize, collapse=","),
                             " --min-count ", poolgenvar.mincount, 
                             " --window ", poolgenvar.window,
                             " --step ", poolgenvar.step, 
                             " --sitecount truewindows-200000-200000.txt",
                             " --min-sites-frac ", poolgenvar.minsitesfrac,
                             " --output Popgen", sep="")
system(PoolGen_var.command)


## DrosEU part E starts here

##################################
## Identify SNPs inside introns ##
##################################

intronicsnps.command <- paste("python ", script.folder, "/IntronicSnps.py",
                              " --gff dmel-all-filtered-r6.09.gff.gz",
                              " --sync SNPs.sync.gz",
                              " --target-length ", intronicsnps.targetlength,
                              " > intron60_all.sync", sep="")

system(intronicsnps.command)

#######################
## Filter SNPs again ##
#######################
## Remove SNPs within a user defined distance of chromosomal inversions 
## and locations with low recombination rates

filterinversions.command <- paste("python ", script.folder, "/FilterByRecomRateNInversion.py",
                                  " --inv ", genome.folder, "/inversions_breakpoints_v5v6.txt",
                                  " --RecRates ", genome.folder, "/DrosEU-SNPs.recomb.gz",
                                  " --input intron60_all.sync",
                                  " --D", filterinversions.D,
                                  " --r", filterinversions.r,
                                  " > intron60.sync", sep="")
                                  
system(filterinversions.command)

############################################
## Pairwise FST (Weir and Cockerham 1984) ##
############################################

fst.command <- paste("python ", script.folder, "/FST.py",
                     " --pool-size ", fst.poolsize,
                     " --input ", "intron60.sync",
                     " --minimum-count ", fst.minimumcount, 
                     " --minimum-cov ", fst.minimumcov,
                     " | gzip > intron.fst.gz", sep="")
system(fst.command)
                     
#################################
## Average FST across all loci ##
#################################

combinefst.command <- paste("python ", script.folder, "/CombineFST.py",
                            " --diff ", "intron.fst.gz",
                            " --stat ", combinefst.stat,
                            " > intron_average.fst",  sep="")
system(combinefst.command)

#################################
## Isolation by distance (IBD) ##
#################################

ibd.command <- paste("python ", script.folder, "/IBD.py",
                     " --fst ", "intron_average.fst",
                     " --meta ", genome.folder, "/MetaData.txt",
                     " --output ", "IBD_EU", sep="")
system(ibd.command)

##############################
## Major allele frequencies ##
##############################

sync2AF.command <- paste("python ", script.folder, "/sync2AF.py",
                         " --input intron60.sync",
                         " --output intron60-af", sep="")
system(sync2AF.command)






