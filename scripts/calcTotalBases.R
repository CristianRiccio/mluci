# Calculate the total length of the M. luci subreads

# Read in the output of seqtk comp
seqtkComp <- read.delim('output/reads/dna/pacbio/all.subreads.seqtkComp.txt',
                        header = FALSE)

# Name the columns of the seqtkComp data frame
colnames(seqtkComp) <- c('chr', 'length', 'A', 'C', 'G', 'T', 'two', 'three', 'four', 'CpG', 'tv', 'ts', 'CpGts')

sum(seqtkComp$length)
