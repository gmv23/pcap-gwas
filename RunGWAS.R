setwd("~/Documents/work/Smart_lab/P_capsici/isolate_collection/paper/gwas/")
library(rrBLUP)
library(qqman)
library(GENESIS)
library(pcaMethods)
library(GWASTools)
library(SNPRelate)

#################################   Load and clean data   ################################

#Reformat VCF file as GDS file
snpgdsVCF2GDS(vcf.fn="capsici_diversity_PG.recode.vcf", out.fn="capsici_diversity_PG.gds")

#Load data
geno <- read.table("../pop_structure/capsici_diversity_PG.012")
indvs <- read.table("../pop_structure/capsici_diversity_PG.012.indv", stringsAsFactors = F)
snps <- read.table("../pop_structure/capsici_diversity_PG.012.pos")
phenos <- read.csv("../phenotypes_and_clones/tables/isolate_metadata.csv")
geno.gds <- GdsGenotypeReader("capsici_diversity_PG.gds")
  
#Clean up sample names
indvs <-unlist(indvs$V1)
indvs <- sapply(indvs, function(x) unlist(strsplit(x,":"))[1])

#Get missing values as NAs in genotype data
geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#Clean up snps data frame
colnames(snps) <- c("CHROM", "BP")
snps$MARKER <- paste(snps$CHROM, snps$BP, sep="_")

#Put phenos in sample order
phenos <- phenos[match(indvs, phenos$SampleSZ),]

#################################   Get phenotypes and covariates   ################################

#Cor between 5 and 100 ug/ml mef RG phenotypes
cor(phenos$Mef5, phenos$Mef100, use = "complete.obs")

#Get clean binary phenotypes
phenos$MT_code <- ifelse(phenos$MT=="A1",1,0)
phenos$Mef_code <- ifelse(phenos$Mef5 > 0.10,1,0)

#Get first two PCs of geno matrix
geno.pca <- pca(geno, method = "nipals", nPcs = 4, scale = 'uv')

#Make ScanAnnotationDataFrame object with phenotypes and covariates
scanAnnot.df <- data.frame("scanID" = phenos$SampleSZ,
                           "PC1" = geno.pca@scores[,1],
                           "PC2" = geno.pca@scores[,2],
                           "PC3" = geno.pca@scores[,3],
                           "PC4" = geno.pca@scores[,4],
                           "MT" = phenos$MT_code,
                           "Mef" = phenos$Mef_code)
scanAnnot <- ScanAnnotationDataFrame(scanAnnot.df)

#Make GenotypeData object
genoData <- GenotypeData(geno.gds, scanAnnot = scanAnnot)

#Get genomic relationship matrix
K <- A.mat(geno-1)
rownames(K) <- scanAnnot.df$scanID
colnames(K) <- scanAnnot.df$scanID

#Association of PCs with MT and mefenoxam
t.test(geno.pca@scores[,1] ~ phenos$MT)
t.test(geno.pca@scores[,2] ~ phenos$MT)
t.test(geno.pca@scores[,3] ~ phenos$MT)
t.test(geno.pca@scores[,4] ~ phenos$MT)
t.test(geno.pca@scores[,1] ~ phenos$Mef_code)
t.test(geno.pca@scores[,2] ~ phenos$Mef_code)
t.test(geno.pca@scores[,3] ~ phenos$Mef_code)
t.test(geno.pca@scores[,4] ~ phenos$Mef_code)

###########################   Run mefenoxam association tests   ##########################

#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "Mef",
                           family = binomial)
mod_K <- fitNullModel(scanAnnot, outcome = "Mef",
                      cov.mat = K, family = binomial)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1"),
                        cov.mat = K, family = binomial)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1"),
                       family = binomial)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1","PC2"),
                        cov.mat = K, family = binomial)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1", "PC2"),
                       family = binomial)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = binomial)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1", "PC2", "PC3"),
                       family = binomial)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = binomial)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "Mef", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = binomial)
pdf("plots/Mefenoxam_models.pdf")
old.par <- par(no.readonly = T)
par(xpd=NA)
plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
       mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n", main = "Mefenoxam sensitivity models")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                            "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC) + 30,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
       mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))
par(old.par)
dev.off()

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests for Mefenoxam
assocMEF <- assocTestSingle(gdsobj=genoIterator, null.model=mod_KQ3)
assocMEF$chr <- snps$CHROM #Provide correct chromosome names

###########################   Run mating type association tests   ###########################

#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "MT",
                           family = binomial)
mod_K <- fitNullModel(scanAnnot, outcome = "MT",
                      cov.mat = K, family = binomial)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1"),
                        cov.mat = K, family = binomial)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1"),
                       family = binomial)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1","PC2"),
                        cov.mat = K, family = binomial)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1", "PC2"),
                       family = binomial)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = binomial)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1", "PC2", "PC3"),
                       family = binomial)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = binomial)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "MT", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = binomial)
pdf("plots/MatingType_models.pdf")
old.par <- par(no.readonly = T)
par(xpd=NA)
plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n", main = "Mating type sensitivity models")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC) + 30,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))
par(old.par)
dev.off()

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

assocMT <- assocTestSingle(gdsobj=genoIterator, null.model=mod_simple)
assocMT$chr <- snps$CHROM #Provide correct chromosome names

##########################   Close gds connection and save results   #########################

close(geno.gds)

write.csv(assocMEF, "tables/MEF_GWAS.csv", quote=F, row.names = F)
write.csv(assocMT, "tables/MT_GWAS.csv", quote=F, row.names = F)


##########################   Random calculations   #########################
table(phenos$MT)
mt_snps <- assocMT[assocMT$Score.pval < .05/nrow(assocMT),]
nrow(mt_snps)
range(mt_snps$pos)
diff(range(mt_snps$pos))
mt_snps[which(mt_snps$Score.pval==min(mt_snps$Score.pval)),]

#Look at predictions at peak SNP
peak_geno <- geno[,which(mt_snps$Score.pval==min(mt_snps$Score.pval))]
peak_mt <- data.frame("sample" = phenos$SampleSZ,
                      "mt" = phenos$MT,
                      "geno" = geno[,which(assocMT$Score.pval==min(assocMT$Score.pval))])
peak_mt$match <- (peak_mt$mt == 'A1' & peak_mt$geno == 0) | (peak_mt$mt=='A2' & peak_mt$geno == 1)

#Look at mefenoxam SNPs
mef_snps <- assocMEF[assocMEF$Score.pval < .05/nrow(assocMEF),]
range(mef_snps$pos)
diff(range(mef_snps$pos))
geno_res <- geno[,which(assocMEF$Score.pval == min(assocMEF$Score.pval))]
phenos$Mef5[geno_res==2]
phenos$Mef100[geno_res==2]
geno_res_MAF <- (sum(geno_res==1, na.rm = T) + 2*sum(geno_res==2, na.rm = T)) / sum(!is.na(geno_res))
meta <- read.csv("../phenotypes_and_clones/isolate_plotting_metadata.csv", stringsAsFactors = F)
meta <- meta[match(phenos$SampleSZ, meta$SampleSZ),]
unique(meta$Field[geno_res == 1 | geno_res == 2])
meta$Field[geno_res == 1 | geno_res == 2]

unique(meta$Field[geno_res == 2])
meta$Field[geno_res == 2]

sort(phenos$Mef5[geno_res==0])
sort(phenos$Mef100[geno_res==0])

sort(phenos$Mef5[geno_res==2])
sort(phenos$Mef100[geno_res==2])

phenos$SampleSZ[phenos$Mef5>.75 & geno_res==0]
phenos$SampleSZ[phenos$Mef5==0 & geno_res==2]

#Fst of these 6 SNPs
get_p <- function(x){
  x <- x[!is.na(x)]
  return((2*sum(x==0) + sum(x==1))/(length(x)*2))
}

hudson <- function(x,y){
  p1 <- get_p(x)
  p2 <- get_p(y)
  n1 <- 2*sum(!is.na(x))
  n2 <- 2*sum(!is.na(y))
  num1 <- (p1 - p2)^2
  num2 <- (p1 * (1-p1))/(n1-1)
  num3 <- (p2 * (1-p2))/(n2-1)
  den <- (p1*(1-p2)) + (p2*(1-p1))
  num <- num1 - num2 - num3
  return(num/den)
}

fields <- c("Cayuga #1", "Erie #1", "Erie #2", "Ontario #1", "Suffolk #1")
for(snp in which(assocMEF$Score.pval < .05/nrow(assocMEF))){
  fst_mat <- matrix(NA, nrow=5, ncol=5)
  for(i in 1:4){
    for(j in ((i+1):5)){
      x <- geno[phenos$Field %in% fields[i], snp]
      y <- geno[phenos$Field %in% fields[j], snp]
      fst_mat[i,j] <- hudson(x, y)
    }
  }
  print(assocMEF[snp,])
  print(fst_mat)
}
  
#Lets look at haplotype at mef region
mef_genos <- geno[,assocMEF$Score.pval < .05/nrow(assocMEF)]
mef_haplos <- apply(mef_genos, 1, paste, collapse="")
mef_haplos[grep('NA', mef_haplos)] <- NA
stripchart(phenos$Mef5 ~ mef_haplos, vertical=T, method='jitter')

#Look at mt snps
geno_mt <- geno[,assocMT$Score.pval < .05/nrow(assocMT)]
geno_mt_a1 <- geno_mt[phenos$MT=="A1",]
geno_mt_a2 <- geno_mt[phenos$MT=="A2",]
sum(geno_mt_a1==0, na.rm=T)/sum(geno_mt_a1==2 | geno_mt_a1==0, na.rm=T)

freq_homo <- function(x){
  x <- x[!is.na(x)]
  x <- sum(x %in% c(0,2))/length(x)
}

a1_homo <- apply(geno_mt_a1, 1, freq_homo)
a2_homo <- apply(geno_mt_a2, 1, freq_homo)
summary(a1_homo)
sum(geno_mt_a1==0, na.rm=T)/sum(geno_mt_a1 %in% c(0,2), na.rm=T)
summary(1-a2_homo)

mt_test <- rbind(geno_mt_a1, geno_mt_a2)
rotate <- function(x) t(apply(x,2,rev))
image(rotate(mt_test), xaxt='n')
axis(side = 1, 
     at = seq(0,1,length.out = sum(assocMT$Score.pval < .05/nrow(assocMT))), 
     labels = snps$BP[assocMT$Score.pval < .05/nrow(assocMT)], 
     las=2, cex.axis=0.75)

