###coombine mean values with new genotype specific column
WT.means <- read.delim("WT_means.txt")
MPIMT.means <- read.delim("MPIMT_means.txt")

WT.means$geno_pair <- paste("WT_", WT.means$interacting_pair, sep = "")
WT.means$geno_pair

MPIMT.means$geno_pair <- paste("MPIMT_", MPIMT.means$interacting_pair, sep = "")
MPIMT.means$geno_pair

joint.means <- rbind(WT.means, MPIMT.means)
joint.means$geno_pair

write_excel_csv(joint.means, "all_means.txt")

###combine p value tables with new genotype specific column
WT.pvals <- read.delim("WT_pvalues.txt")
MPIMT.pvals <- read.delim("MPIMT_pvalues.txt")

WT.pvals$geno_pair <- paste("WT_", WT.pvals$interacting_pair, sep = "")
WT.pvals$geno_pair

MPIMT.pvals$geno_pair <- paste("MPIMT_", MPIMT.pvals$interacting_pair, sep = "")
MPIMT.pvals$geno_pair

joint.pvals <- rbind(WT.pvals, MPIMT.pvals)
joint.pvals$geno_pair

write_excel_csv(joint.pvals, "all_pvalues.txt")



