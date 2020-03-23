# This script is for looking at genotypic discrepancies between isolates of clonal lineage 16
# (3 isolates that are clonal but differ in mating type)

#################################   Load and clean data   ################################

setwd("~/Documents/work/Smart_lab/P_capsici/isolate_collection/paper/gwas/")

#Load scripts for plotting in sliding windows
source("../scripts/plot_genome.R")
source("../scripts/window_smooth.R")

#Reread data prior to removing clones
pcap <- read.csv("../phenotypes_and_clones/tables/isolate_metadata.csv")
indvs <- read.table("capsici_diversity_FILTER.012.indv", stringsAsFactors = F)
geno <- read.table("capsici_diversity_FILTER.012")
snps <- read.table("capsici_diversity_FILTER.012.pos")
depths <- read.table("../ploidy/read_depths.txt", header = T)
mt <- read.csv("tables/MT_GWAS.csv")

#Clean up data
indvs <-unlist(indvs$V1)
indvs <- sapply(indvs, function(x) unlist(strsplit(x,":"))[1])
geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#Put data in same order
pcap <- pcap[match(indvs, pcap$SampleSZ),]

#################################   Subset relevant data and calculate some numbers   ################################

#Subset data out for just clonal lineage 16
pcap.16 <- pcap[pcap$UniqueGenotype==16 & !is.na(pcap$UniqueGenotype==16),]
samples.16 <- pcap.16$SampleSZ
geno.16 <- cbind(snps, t(geno[match(samples.16, indvs),]))
colnames(geno.16) <- c("CHROM", "BP", as.character(samples.16))
geno.16$CHROM <- as.factor(geno.16$CHROM)
#Look at discordant SNPs between individuals
geno.16$EH05A_v_EH26A <- geno.16$`13EH05A` == geno.16$`13EH26A`
geno.16$EH05A_v_EH76A <- geno.16$`13EH05A` == geno.16$`13EH76A`
geno.16$EH76A_v_EH26A <- geno.16$`13EH76A` == geno.16$`13EH26A`

mt_snps <- mt[mt$Score.pval < .05/nrow(mt),]
geno_mt <- which(geno.16$CHROM %in% mt_snps$chr & geno.16$BP > min(mt_snps$pos) & geno.16$BP < max(mt_snps$pos))
geno_sig <- which(geno.16$CHROM %in% mt_snps$chr & geno.16$BP %in% mt_snps$pos & geno.16$BP %in% mt_snps$pos)

discordant <- function(x){
  sum(!x, na.rm=T)/sum(!is.na(x))
}
genomewide_discordant <- apply(geno.16[,6:8], 2, discordant); genomewide_discordant
mt_discordant <- apply(geno.16[geno_mt,6:8], 2, discordant); mt_discordant
mt_sig_discordant <- apply(geno.16[geno_sig,6:8], 2, discordant); mt_sig_discordant
geno_16_sig <- geno.16[geno_sig,]
geno_16_sig_nonmissing <- geno_16_sig[apply(geno_16_sig, 1, function(x) all(!is.na(x))),]
nrow(geno_16_sig_nonmissing)

#################################   Subset relevant data and calculate some numbers   ################################

#Create variables used in plotting function
scaffold_sizes <- aggregate(snps$V2~snps$V1, FUN=max)
colnames(scaffold_sizes) = c("SCAFFOLD", "MAX_SNP")
window_size <- 50000
step_size <- 50000
scaffold_sizes <- scaffold_sizes[scaffold_sizes$MAX_SNP>window_size,]
scaffold_sizes$SCAFFOLD <- as.factor(scaffold_sizes$SCAFFOLD)
geno.16 <- geno.16[geno.16$CHROM %in% scaffold_sizes$SCAFFOLD,]

#Get percentage of discordant calls in windows
geno.16.discordant <- window_smooth(stats=geno.16,data_columns=6:8, FUN=discordant, window_size=window_size, step_size=step_size, 
                                    smooth_by="distance", scaffold_sizes=scaffold_sizes)
#Get SNP counts in windows
snp.counts <- window_smooth(stats=geno.16, data_columns=6:8, FUN = function(x) sum(!is.na(x)),
                            window_size = window_size, step_size = step_size,
                            scaffold_sizes = scaffold_sizes, smooth_by="distance")
colnames(snp.counts)[3:5] <- colnames(geno.16)[6:8]

#Get rid of any smoothed frequency where fewer than 15 snps contribute in window
minsnpcount <- 20
geno.16.discordant.filter <- geno.16.discordant
for(i in 3:5){
  small_windows <- which(snp.counts[,i] < minsnpcount)
  geno.16.discordant.filter[small_windows,i] <- NA
}
nas <- apply(geno.16.discordant.filter[,3:5], 1, function(x) any(is.na(x)))
geno.16.discordant.filter <- geno.16.discordant.filter[!nas,]

#Make plot showing discordant call rate in first 50 scaffolds
geno.16.discordant.filter$CHROM <- droplevels(geno.16.discordant.filter$CHROM)
scaffold_sizes <- scaffold_sizes[scaffold_sizes$SCAFFOLD %in% levels(geno.16.discordant.filter$CHROM),]
scaffold_sizes$SCAFFOLD <- droplevels(scaffold_sizes$SCAFFOLD)
plot_stat(geno.16.discordant.filter[geno.16.discordant.filter$CHROM %in% 1:50,],
          scaffold_sizes=scaffold_sizes[scaffold_sizes$SCAFFOLD %in% 1:50,],
          number_frames=1, tick_spacing = 100000, horizontal_lines=c(),
          stat_colors=adjustcolor(c("blue", "red", "orange"), alpha.f=0.6), stat_lwds=rep(3,3),
          plot_base_name="plots/discordant_snps",
          plot_arguments = list(ylab = "% Discordant Genotypes"),
          legend_arguments = list(x="topright",
                                  legend = c("13EH05A (A2) vs 13EH26A (A1)",
                                             "13EH05A (A2) vs 13EH76A (A2)",
                                             "13EH26A (A1) vs 13EH76A (A2)"),
                                  fill = c("blue", "red", "orange"),
                                  cex=1.5,
                                  bty="n")) ###Note scaffold 45 is not plotted due to no SNPs

##############################   Read depth in LOH region --- haplotype loss or copy neutral?   ##############################
depths.16 <- depths[,which(indvs %in% samples.16) + 2]
depths.16.mt <- depths.16[snps$V1==4 & abs(snps$V2 - 579765) < 400000,]
apply(depths.16,2,mean,na.rm=T)
apply(depths.16.mt,2,mean,na.rm=T)
boxplot(depths.16[,1], depths.16.mt[,1],
        depths.16[,2], depths.16.mt[,2],
        depths.16[,3], depths.16.mt[,3],
        outline=F)



