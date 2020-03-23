#The goal of this script is to make two figures with 4 panels each
#Showing, for both MT and mefenoxam sensitivity,
# NJ tree colored by phenotype
# Histogram/barplot of phenotype distribution
# qqPlot
# Manhattan plot
#In addition:
# Make figure showing allelic effect of most significant SNP for mefenoxam
# LD figures
setwd("~/Documents/work/Smart_lab/P_capsici/isolate_collection/paper/gwas/")
library(ape)
library(qqman)
library(RColorBrewer)

#################################   Load and clean data   ################################

mef <- read.csv("tables/MEF_GWAS.csv")
mt <- read.csv("tables/MT_GWAS.csv")
pcap <- read.csv("../phenotypes_and_clones/tables/isolate_metadata.csv")
tree <- read.nexus("../pop_structure/capsici_diversity_PG.nex")
meta <- read.csv("../phenotypes_and_clones/isolate_plotting_metadata.csv", stringsAsFactors = F)
geno <- read.table("../pop_structure/capsici_diversity_PG.012")
indvs <- read.table("../pop_structure/capsici_diversity_PG.012.indv", stringsAsFactors = F)
snps <- read.table("../pop_structure/capsici_diversity_PG.012.pos")
pheno <- read.csv("../../pheno/Pcap_isolates.csv")
ld <- read.table("meanDecayFig.bin", header = F)

#Clean up data
indvs <-unlist(indvs$V1)
indvs <- sapply(indvs, function(x) unlist(strsplit(x,":"))[1])
geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA
colnames(ld) <- c("Dist", "Meanr2", "MeanD", "Sumr2", "SumD", "NumberPairs")

#Put meta data, tree tip labels, geno, and pheno in same order
pcap <- pcap[match(tree$tip.label, pcap$SampleSZ),]
geno <- geno[match(tree$tip.label, indvs),]
indvs <- indvs[match(tree$tip.label, indvs)]
meta <- meta[match(tree$tip.label, meta$SampleSZ),]
pheno$SampleSZ <- paste(pheno$Sample, pheno$SZ_for_analysis, sep="")
pheno <- pheno[match(tree$tip.label, pheno$SampleSZ),]

#Reorder levels of mef sensitivity
pcap$Mefenoxam_sensitivity = factor(pcap$Mefenoxam_sensitivity,levels(pcap$Mefenoxam_sensitivity)[c(3,1,2)])

#################################   GWAS figure Plotting function   ################################

make_plot <- function(pheno, tree, gwas_results, cols, name, ytop, ybom, x1,x2,x3){
  pdf(paste(name, "_GWAS_figure.pdf", sep=""), width=7, height=4.5)
  old.par <- par(no.readonly = T)
  par(oma=c(2,2,2,2))
  
  #Arrange layout
  plot_layout <- rbind(c(1,2,3),
                       c(4,4,4))
  layout(plot_layout)
  
  #Par settings for plot 1: NJ tree
  par(mar=c(4,4,1,4))
  #Plot 1: NJ tree
  tree_colors <- cols[match(pheno, levels(pheno))]
  tree_colors[is.na(tree_colors)] <- 'gray'
  plot(tree, type="unrooted", show.tip.label = F, 
       no.margin=T, direction="downwards")
  for(i in 1:length(tree$tip.label)){
    tiplabels(tip=i,pch=16,col='black', cex=1.5)
    tiplabels(tip=i,pch=16,col=tree_colors[i], cex=1)
  }
  
  #Par settings for plot 2: Histogram
  par(mar=c(4,5,1,4))
  #Plot 2: Histogram
  barplot(table(pheno), col=cols, ylab = "Counts")
  
  par(xpd=NA)
  text(c(x1,x2,x3, x1),c(ytop, ytop, ytop, ybom),c("A","B","C","D"), cex=1.5)
  
  #Par settings for plot 3: qqplot
  par(mar=c(4,4,1,4), xpd=F)
  #Plot 3: qqplot
  qq(gwas_results$P)
  
  #Par settings for plot 4: Manhattan plot
  par(mar=c(4,4,2,4), xpd=T)
  #Plot
  manhattan(gwas_results, 
            genomewideline=F, 
            suggestiveline = F,
            xlab = "Scaffold")
  par(xpd=F)
  abline(h=-log10(.05/nrow(gwas_results)), col='red')
  #End
  par(old.par)
  dev.off()
  
}

##############################   Make GWAS results plot for MT and Mef   ################################

setwd("./plots/")

#MT
make_plot(pheno = pcap$MT,
          tree = tree,
          gwas_results = data.frame("CHR" = mt$chr,
                                    "BP" = mt$pos,
                                    "P" = mt$Score.pval),
          cols = brewer.pal(3, "PuOr")[c(1,3)],
          name = "MT",
          ytop=72, ybom=-35,
          x1=-6.8,x2=-1.2,x3=3.3)

make_plot(pheno = pcap$Mefenoxam_sensitivity,
          tree = tree,
          gwas_results = data.frame("CHR" = mef$chr,
                                    "BP" = mef$pos,
                                    "P" = mef$Score.pval),
          cols = brewer.pal(3, "YlOrRd"),
          name = "Mef",
          ytop=112, ybom=-58,
          x1=-10.75, x2=-2.75,x3=4.75)

#########################   Make plot showing allelic effect of most significant SNP on mef resistance   ################################

sig_SNP <- which(mef$Score.pval==min(mef$Score.pval))
sig_geno <- geno[,sig_SNP]

pheno$Mef5 <- 100*.5*((pheno$Mef_5ug_diameter1 - 10) / (pheno$Mef_0ug_diameter1 - 10) + 
                        (pheno$Mef_5ug_diameter2 - 10) / (pheno$Mef_0ug_diameter2 - 10))
pheno$Mef100 <- 100*.5*((pheno$Mef_100ug_diameter1 - 10) / (pheno$Mef_0ug_diameter1 - 10) + 
                          (pheno$Mef_100ug_diameter2 - 10) / (pheno$Mef_0ug_diameter2 - 10))

pdf("Mefenoxam_effect.pdf", width=3.25, height=5)
layout(plot_layout)
old.par <- par(no.readonly = T)
par(mar=c(1,4,0,0), oma=c(4,2,1,1), mfrow=c(2,1), xpd=NA)
boxplot(pheno$Mef5~sig_geno, 
        outline=F, border='gray0', boxwex=0.5, staplewex=0,
        xaxt="n", xlab = "",
        ylab = expression(atop("% Relative Growth", 
                               paste("(5 ", mu, "g/ml mefenoxam)"))),        cex.axis=0.7,
        cex.lab=0.8)
points(x = sig_geno+rnorm(length(sig_geno), 0, 0.025) +1, 
       y = pheno$Mef5,
       col=meta$color, pch=meta$pch, cex=1, lwd=1.5)
par(mar=c(1,4,0,0))
boxplot(pheno$Mef100~sig_geno, 
        outline=F, border='gray0', boxwex=0.5, staplewex=0,
        xaxt="n", 
        xlab = "Genotype",
        ylab = expression(atop("% Relative Growth", 
                         paste("(100 ", mu, "g/ml mefenoxam)"))),
        cex.axis=0.7,
        cex.lab=0.8)    
points(x = sig_geno+rnorm(length(sig_geno), 0, 0.025) +1, 
       y = pheno$Mef100,
       col=meta$color, pch=meta$pch, cex=1, lwd=1.5)
axis(1, at=c(1,2,3), labels = c("AA", "AC", "CC"), cex.axis=0.75)

par(xpd=NA)
text(c(-0.75,-0.75),c(225,110), c("A","B"), cex=1.25)
par(old.par)
dev.off()

#############################   Make LD decay plot   ################################

pdf("ld_decay.pdf")
old.par <- par(no.readonly = T)
par(mar=c(5,5,3,5))
plot(ld$Dist, ld$Meanr2, 
     type='l', xaxt='n', xlim=c(0,500000),
     xlab = 'Distance (kb)',
     ylab = expression(paste("Mean ", "r"^"2")))
axis(1, at=seq(0,500000, by = 100000), 
     labels = seq(0,500, by=100))
abline(h=0.04, col='red', lty=2)
abline(h=0.10, col='red', lty=2)
par(mar=c(5,4,4,6), xpd=NA)
text(560000, 0.10, labels = expression(paste("r"^"2", " = 0.10")))
text(560000, 0.04, labels = expression(paste("r"^"2", " = 0.04")))
dev.off()

#########################   Show LD to peak SNPs in significant regions   ################################

####### First MT
#400 Kb on either side of peak SNP
scaff4 <- mt[mt$chr==4 & mt$pos > 179765 & mt$pos < 979765,]

#get -log10 pvalues
p.scaff4 <- -log10(scaff4$Score.pval)

#Get LD to peak SNP
geno.scaff4 <- geno[,mt$chr==4 & mt$pos > 179765 & mt$pos < 979765]
peak <- which(scaff4$Score.pval == min(scaff4$Score.pval))
geno.peak <- geno.scaff4[,peak]
lds <- rep(NA, nrow(scaff4))
for(i in 1:length(lds)){
  lds[i] <- cor(geno.scaff4[,i], geno.peak, use = "complete.obs")^2
}

#rescale r2 distribution to -log10 distribution for plotting
lds.rescale <- (max(p.scaff4)-min(p.scaff4))/(max(lds)-min(lds)) * (lds-max(lds)) + max(p.scaff4)

pdf("MT_LD.pdf")
old.par <- par(no.readonly=T)
par(mar=c(5,5,3,5))
plot(0, type="n", xlim = range(scaff4$pos), ylim=c(0, max(p.scaff4)), xaxt="n",
     ylab = expression(paste("-log"[10], "(", italic("p"), ")")),
     xlab = "Physical position (Kb) on scaffold 4")
for(r in 1:nrow(scaff4)){
  lines(x=rep(scaff4$pos[r],2), y = c(0,p.scaff4[r]), col='gray')
}
points(scaff4$pos, lds.rescale, pch=2)
points(scaff4$pos[peak], lds.rescale[peak], pch=17, cex=0.9, col='orange')
axis(side=4, at = seq(0,max(p.scaff4), length.out=6), labels = seq(0,max(lds), length.out=6))
axis(side=1, at = seq(100000,1000000,by=100000), labels = seq(100000,1000000,by=100000)/1000)
mtext(expression(italic("r")^2), side=4, line=3)
par(old.par)
dev.off()


####### Now Mef
#400 Kb on either side of peak SNP
scaff62 <- mef[mef$chr==62,]

#get -log10 pvalues
p.scaff62 <- -log10(scaff62$Score.pval)

#Get LD to peak SNP
geno.scaff62 <- geno[,mt$chr==62]
peak <- which(scaff62$Score.pval == min(scaff62$Score.pval))
geno.peak <- geno.scaff62[,peak]
lds <- rep(NA, nrow(scaff62))
for(i in 1:length(lds)){
  lds[i] <- cor(geno.scaff62[,i], geno.peak, use = "complete.obs")^2
}

#rescale r2 distribution to -log10 distribution for plotting
lds.rescale <- (max(p.scaff62)-min(p.scaff62))/(max(lds)-min(lds)) * (lds-max(lds)) + max(p.scaff62)

pdf("Mef_LD.pdf")
old.par <- par(no.readonly=T)
par(mar=c(5,5,3,5))
plot(0, type="n", xlim = range(scaff62$pos), ylim=c(0, max(p.scaff62)), xaxt="n",
     ylab = expression(paste("-log"[10], "(", italic("p"), ")")),
     xlab = "Physical position (Kb) on scaffold 62")
for(r in 1:nrow(scaff62)){
  lines(x=rep(scaff62$pos[r],2), y = c(0,p.scaff62[r]), col='gray')
}
points(scaff62$pos, lds.rescale, pch=2)
points(scaff62$pos[peak], lds.rescale[peak], pch=17, cex=0.9, col='orange')

axis(side=4, at = seq(0,max(p.scaff62), length.out=6), labels = seq(0,max(lds), length.out=6))
axis(side=1, at = seq(0, 400000,by=50000), labels = seq(0,400000,by=50000)/1000)
mtext(expression(italic("r")^2), side=4, line=3)
par(old.par)
dev.off()

