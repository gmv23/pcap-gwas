#!/usr/bin/env Rscript

#Read arguments
args <- commandArgs(trailingOnly=TRUE)
pos_file <- args[1]
geno_file <- args[2]
out_file <- args[3]
het_cutoff <- as.numeric(args[4])

#Get SNP names
pos <- read.table(pos_file)
snps <- apply(pos, 1, function(x) paste("S", x[1], "_", x[2], sep=""))

#Find excess hets
geno <- read.table(geno_file)
geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA
het_rate <- apply(geno,2,function(x) sum(x==1, na.rm=T)/sum(!is.na(x)))
excess <- which(het_rate > het_cutoff)

#Print file with SNP IDs of excess hets
excess_snp_ids <- snps[excess]
write.table(excess_snp_ids, out_file, row.names=F, col.names=F, quote=F)
