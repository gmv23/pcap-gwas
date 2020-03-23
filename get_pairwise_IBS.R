#!/usr/bin/env Rscript

# This script takes a VCF file and creates a pairwise IBS similarity matrix 
# for use in identifying and filtering out clones

################################## LOAD 012 FILES ###############################

snps <- read.table("./capsici_diversity_FILTER.012.pos")
head(snps)

indvs <- read.table("./capsici_diversity_FILTER.012.indv", stringsAsFactors = F)
indvs <-unlist(indvs$V1)

geno <- read.table("./capsici_diversity_FILTER.012")
print("geno loaded")
geno <- geno[,-1]
geno <- t(geno)

geno[geno==-1] <- NA
print('finished replacing NAs')

################################ IBS FUNCTION ###################################

ibs <- function(x,y){

  alleles_ibs <- 2 - abs(x-y)
  return(sum(alleles_ibs, na.rm = T)/(2*sum(!is.na(alleles_ibs))))
  
}

#################### CALCULATE IBS FOR EACH PAIRWISE ISOLATE COMBINATION ###########

d <- ncol(geno)

IBS_matrix <- matrix(nrow=d, ncol=d)

print("got to loop")

for(i in 1:(d-1)){
	for (j in (i +1):d){
		IBS_matrix[i,j] <- ibs(geno[,i], geno[,j])
	}
	print(i)
}

rownames(IBS_matrix) <- indvs
colnames(IBS_matrix) <- indvs

write.csv(IBS_matrix, "IBS_matrix.csv")
