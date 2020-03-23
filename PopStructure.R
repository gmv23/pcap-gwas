setwd("~/Documents/work/Smart_lab/P_capsici/isolate_collection/paper/pop_structure/")

library(pcaMethods)
library(ape)
library(RColorBrewer)
library(StAMPP)

##############################################         Clean and load data       ######################################################

#Load data
geno <- read.table("capsici_diversity_PG.012")
indvs <- read.table("capsici_diversity_PG.012.indv", stringsAsFactors = F)
snps <- read.table("capsici_diversity_PG.012.pos")
pcap <- read.csv("../phenotypes_and_clones/tables/isolate_metadata.csv", stringsAsFactors = F)
meta <- read.csv("../phenotypes_and_clones/isolate_plotting_metadata.csv", stringsAsFactors = F)
field_info <- read.csv("../phenotypes_and_clones/field_plotting_metadata.csv", stringsAsFactors = F)
tree <- read.nexus("capsici_diversity_PG.nex")

#Clean up data
indvs <-unlist(indvs$V1)
indvs <- sapply(indvs, function(x) unlist(strsplit(x,":"))[1])

geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#Subset metadata for just cc set and get in same order as geno
pcap <- pcap[match(indvs, pcap$SampleSZ),]
meta <- meta[match(indvs, meta$SampleSZ),]

######################################################     PCA       #################################################################  

#Get PCs
geno.pca <- pca(geno, method = "nipals", nPcs = 4, scale = 'uv')
pve <- round(geno.pca@R2*100,2)

##### Make plot of PC1v2 and PC3v4
pdf("plots/PCA.pdf", width=6, height=6)
plot_layout <- cbind(c(1,1, 2, 2),
                     c(3, 3, 3, 3))
layout(plot_layout)
#Par settings for plots 1 and 2 (PCA plots)
old.par <- par(no.readonly = T)
par(mar=c(4,1,2,1), oma=c(4,4,4,4), xpd=NA)

#PC1 vs PC2
plot(geno.pca@scores[,1],
     geno.pca@scores[,2],
     pch=meta$pch,cex=1,col=meta$color,
     xlab = paste("PC1:", pve[1], "%"),
     ylab = paste("PC2:", pve[2], "%"))

#PC3 vs PC4
plot(geno.pca@scores[,3],
     geno.pca@scores[,4],
     pch=meta$pch,cex=1,col=meta$color,
     xlab = paste("PC1:", pve[3], "%"),
     ylab = paste("PC2:", pve[4], "%"))

#Plot 3: legend
populations <- meta$Field
plot(0,"n", bty="n", xlim=c(1,10), ylim=c(1,10), xaxt='n', yaxt='n', xlab='', ylab='')
legend_starts <- c(11, 9, 6, 4, 2.5, 1.5)
site_counts <- table(populations)
for(i in 1:length(unique(field_info$Region))){
  region <- unique(field_info$Region)[i]
  sites <- field_info$Field[field_info$Region==region]
  colors <- field_info$color[field_info$Region==region]
  pch <- field_info$pch[field_info$Region==region]
  legend(legend_starts[i],legend=paste(sites, ", n=", site_counts[match(sites, names(site_counts))], sep=""),
         col = colors,
         pch = pch,
         bty='n', title = as.expression(bquote(bold(.(region)))),
         title.adj = 0, cex=0.9)
}
par(old.par)
dev.off()

######################################################     NJ Tree       #################################################################  

tip_pops <- populations[match(tree$tip.label, pcap$SampleSZ)]
tip_colors <- field_info$color[match(tip_pops, field_info$Field)]
tip_pch <- field_info$pch[match(tip_pops, field_info$Field)]
pdf("plots/labeled_tree.pdf")
plot(tree, type="unrooted", show.tip.label = T, lab4ut = "axial",
     no.margin=T, direction="downwards", cex=0.75)
dev.off()

####################################       PCA and tree in one plot       ####################################################################  

#### New layout
pdf("plots/PCA_and_tree.pdf", width=7, height=4.75)
plot_layout <- rbind(c(1, 2),
                     c(1, 2),
                     c(3, 3))

layout(plot_layout)
#Par settings for plot 1: PC 1 vs 2
old.par <- par(no.readonly = T)
par(mar=c(4,2,2,2), oma=c(1,4,1,2), xpd=NA)

#Plot 1: PC 1 vs 2
plot(geno.pca@scores[,1],
     geno.pca@scores[,2],
     pch=meta$pch,cex=1.5,col=meta$color, lwd=1.5,
     xlab = paste("PC1:", pve[1], "%"),
     ylab = paste("PC2:", pve[2], "%"))

#Par settings for plot 2: tree
par(mar=c(4, 1, 2, 2))
#Plot 2: Tree
plot(tree, type="unrooted", show.tip.label = F, 
     no.margin=T, direction="downwards")
tiplabels(pch=tip_pch,col=tip_colors, cex=1.5, lwd=1.5)

#Par settings for plot 3: legend
par(mar=c(3,1,2.5,2))
#Empty plot for legend
plot(0,"n", bty="n", xlim=c(1,10), ylim=c(1,10), xaxt='n', yaxt='n', xlab='', ylab='')

#Get coordinates of where each region is oging to go
x_starts <- c(0, 3, 6, 6, 9, 9)
y_starts <- c(13.5, 13.5, 13.5, 3, 13.5, 8)
site_counts <- table(populations)

#Loop through each region
region_order <- c("CD", "CNY", "LI", "WNY", "NY unknown", "Non NY")
for(i in 1:length(region_order)){
  region <- region_order[i]
  sites <- field_info$Field[field_info$Region==region]
  colors <- field_info$color[field_info$Region==region]
  pch <- field_info$pch[field_info$Region==region]
  legend(x_starts[i], y_starts[i], legend=paste(sites, ", n=", site_counts[match(sites, names(site_counts))], sep=""),
         col = colors,
         pch = pch,
         bty='n', title = as.expression(bquote(bold(.(region)))),
         title.adj = 0, cex=0.95)
}

text(c(0,6),c(50,50), c("A", "B"), cex=2)

par(old.par)
dev.off()

####################################       FST between fields       ####################################################################  

#Get fields with more than 6 isolates after clone-correction
populations.sizes <- table(populations)
large_pops <- names(populations.sizes)[populations.sizes>=6 & names(populations.sizes)!="NonNY"]
fields <- populations
fields[!fields %in% large_pops] <- NA

#Turn geno into allele frequency format
geno.stamp <- geno/2

#Make column names SNP names
colnames(geno.stamp) <- paste("S_",snps$V1, "_", snps$V2, sep="")

#Add isolate, population, ploidy, and format information 
geno.stamp <- cbind("Isolate" = pcap$SampleSZ,
                    "Subpop" = fields,
                    "Ploidy" = 2,
                    "Format" = 'freq',
                    geno.stamp)
geno.stamp <- as.data.frame(geno.stamp, stringsAsFactors = F)
geno.stamp <- geno.stamp[!is.na(geno.stamp$Subpop),]

#Make object for StAMPP functions
geno.stamp <- stamppConvert(geno.stamp, "r")

#Weir and Cockerham's Fst, 100 bootstraps, 95% CI
fst.stamp <- stamppFst(geno.stamp,100,95)
fst.mat <- fst.stamp$Fsts

#Re-order matrix
order_matrix <- function(mat, ord){
  
  mat.sort <- matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
  rownames(mat.sort) <- ord
  colnames(mat.sort) <- ord
  
  for(i in 1:(nrow(mat.sort)-1)){
    ord.i <- ord[i]
    for(j in (i+1):ncol(mat.sort)){
      ord.j <- ord[j]
      mat.values <- c(mat[ord.i,ord.j], mat[ord.j,ord.i])
      mat.value <- mat.values[which(!is.na(mat.values))]
      mat.sort[i,j] <- mat.value
    }
  }
  return(mat.sort) 
}
fst.mat <- order_matrix(fst.mat,sort(rownames(fst.mat)))
field.sizes <- populations.sizes[match(rownames(fst.mat), names(populations.sizes))]
fst.mat <- cbind("Number isolates" = field.sizes, fst.mat)
write.table(fst.mat, "tables/pairwise_fst_ont13_and_17.txt", row.names = T, quote=F, sep="\t")

#### Now do it again but collapse ontario 2013 and 2017 into one population

#Get fields with more than 6 isolates after clone-correction
populations.sizes <- table(populations)
large_pops <- names(populations.sizes)[populations.sizes>=6 & names(populations.sizes)!="NonNY"]
fields <- populations
fields[!fields %in% large_pops] <- NA
fields[fields %in% c("Ontario #1 (2017)", "Ontario #1 (2013)")] <- "Ontario #1 (2013, 2017)"

#Turn geno into allele frequency format
geno.stamp <- geno/2

#Make column names SNP names
colnames(geno.stamp) <- paste("S_",snps$V1, "_", snps$V2, sep="")

#Add isolate, 'population' (year), ploidy, and format information 
geno.stamp <- cbind("Isolate" = pcap$SampleSZ,
                    "Subpop" = fields,
                    "Ploidy" = 2,
                    "Format" = 'freq',
                    geno.stamp)
geno.stamp <- as.data.frame(geno.stamp, stringsAsFactors = F)
geno.stamp <- geno.stamp[!is.na(geno.stamp$Subpop),]

#Make object for StAMPP functions
geno.stamp <- stamppConvert(geno.stamp, "r")

#Weir and Cockerham's Fst, 100 bootstraps, 95% CI
fst.stamp <- stamppFst(geno.stamp,100,95)
fst.mat <- fst.stamp$Fsts

fst.mat <- order_matrix(fst.mat,sort(rownames(fst.mat)))
field.sizes <- table(fields)
field.sizes <- field.sizes[order(names(field.sizes))]
fst.mat <- cbind("Number isolates" = field.sizes, fst.mat)
fst.mat <- round(fst.mat, 2)
write.table(fst.mat, "tables/pairwise_fst.txt", row.names = T, quote=F, sep="\t")

