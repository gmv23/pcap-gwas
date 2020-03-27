# Goal of this script is to first use allele balances to determine possible triploid LGs
# Then look at read depths and confirm that triploid LGs differ from diploid LGs for given isolate
# Then look at some trends with the ploidy determinations across all isolates

setwd("~/Documents/work/Smart_lab/P_capsici/isolate_collection/paper/ploidy/")
library(RColorBrewer)

############################################   LOAD AND CLEAN DATA   ######################################
depths <- read.table("read_depths.txt", header = T)
abs <- read.table("allele_balances.txt", header = T)
pcap <- read.csv("~/Documents/work/Smart_lab/P_capsici/isolate_collection/paper/phenotypes_and_clones/tables/isolate_metadata.csv")
lg_order <- read.csv("scaffold_assignments.csv", header=T, stringsAsFactors = F) #Table S12 from Lamour et al 2012
subscaffolds <- read.csv("subscaffold_assignments.csv", header=T, stringsAsFactors = F) #Table S11 from Lamour et al 2012

lg_order <- lg_order[-nrow(lg_order),] #Get rid of totals row 

#Function to divide string like 'Sc45.2' or 'LG01.01' into first or second number
get_value <- function(x, position){
  x.number.start <- regexpr("[0-9]", x) #Find where numbers start
  x.clean <- substr(x,x.number.start,nchar(x)) #Remove beginning letters
  x.split <- strsplit(x.clean, ".", fixed = T)[[1]] #Split on "."
  return(as.integer(x.split[position])) #Return first or second number
}

#Function to get marker (2nd part) from string like "1_529343"
get_marker <- function(x){
  start_pos <- regexpr("_", x)
  marker <- as.integer(substr(x, start_pos + 1, nchar(x)))
  return(marker)
}

#Separate Scaffold_Block into Scaffold and Block
lg_order$Scaffold <- sapply(lg_order$Scaffold_Block, get_value, position = 1)
lg_order$Block <- sapply(lg_order$Scaffold_Block, get_value, position = 2)

subscaffolds$Scaffold <- sapply(subscaffolds$Scaffold_Block, get_value, position = 1)
subscaffolds$Block <- sapply(subscaffolds$Scaffold_Block, get_value, position = 2)

#Clean up marker positions
subscaffolds$Lowest_Marker <- sapply(subscaffolds$Lowest_Marker, get_marker)
subscaffolds$Highest_Marker <- sapply(subscaffolds$Highest_Marker, get_marker)

##################################   COMBINE GENETIC MAP DATASETS INTO ONE DATA FRAME   #######################

print("Pulling information from genetic map tables to assign scaffolds to linkage groups")

block_assignments <- matrix(NA, nrow=nrow(lg_order), ncol=6)
colnames(block_assignments) <- c("Scaffold", "Block", "Low", "High", "LG", "LG_order")

#Loop through lg_order by row
#Take out values or match values in subscaffolds data frame to populate block_assignments data frame
for(i in 1:nrow(block_assignments)){
  row_data <- unlist(lg_order[i,]) 
  #Pull out values from each row of lg_order
  scaffold <- as.integer(row_data[5])
  block <- as.integer(row_data[6])
  lg <- as.integer(get_value(row_data[1], 1))
  lg.pos <- as.integer(get_value(row_data[1], 2))
  #If scaffold isn't divided into subscaffold blocks, "low" and "high" markers cover entire scaffold
  if(is.na(block)){ 
    low <- 1
    high <- as.integer(row_data[3])
    #Otherwise find matching marker positions of subscaffold from subscaffolds data frame
  }else{
    low <- subscaffolds$Lowest_Marker[subscaffolds$Scaffold==scaffold & subscaffolds$Block==block]
    high <- subscaffolds$Highest_Marker[subscaffolds$Scaffold==scaffold & subscaffolds$Block==block]
  }
  block_assignments[i,] <- c(scaffold, block, low, high, lg, lg.pos) #Fill in row in block_assignments df
}

block_assignments <- as.data.frame(block_assignments)

#Look at block_assignments in order by scaffold and scaffold block
head(block_assignments[order(block_assignments$Scaffold, block_assignments$Block),])

############################################   ASSIGN POSITIONS IN VCF FILE TO LGs    ##################################

print("Assign markers to linkage groups")

# Insert LG column to depths data frame after chromosome and position
depths <- data.frame(depths[1:2], "LG"=rep(NA, nrow(depths)), depths[3:ncol(depths)])
abs <- data.frame(abs[1:2], "LG"=rep(NA, nrow(abs)), abs[3:ncol(abs)])

#Make sure depths and abs data frame contain same markers
any(depths[,1:2] != abs[,1:2])

#Loop through rows in depths data frame
#Use CHROM and POS values to find matching LG from block_assignments data frame
for(i in 1:nrow(depths)){
  chrom <- depths[i,1]
  pos <- depths[i,2]
  #How many times does scaffold appear in block_assignments
  #0 means not assigned to lg, 1 means not divided into subscaffolds, 2+ means divided into subscaffold blocks
  scaffold_block_number <- sum(block_assignments$Scaffold==chrom) 
  if(scaffold_block_number==0){ #If scaffold not assigned a LG, return NA
    lg <- NA
  }else if(scaffold_block_number==1){ #If no subdivision of scaffold, match LG
    lg <- block_assignments$LG[block_assignments$Scaffold==chrom]
  }else if(scaffold_block_number > 0){ #If subdivision of scaffold, find appropriate scaffold block
    lg <- NA #NA if POS is not in any of the scaffold blocks
    chrom_blocks <- block_assignments[block_assignments$Scaffold==chrom,] #Pull out info on scaffold blocks
    for(j in 1:nrow(chrom_blocks)){ #Loop through scaffold blocks
      if(pos >= chrom_blocks$Low[j] & pos <= chrom_blocks$High[j]){
        lg <- chrom_blocks$LG[j] #If marker position is in between low and high, pull out matching LG
      }
    }
  }
  depths$LG[i] <- lg #Populate LG column of depths data frame
  abs$LG[i] <- lg #Populate LG column of abs data frame
}

###############################################   DRAW HISTOGRAMS FOR ALL ISOLATES AND LGS    ######################################

#Function to pull out abbreviated isolate name
get_name <- function(x){
  name <- unlist(strsplit(x, "[:.]"))[1]
  if(substr(name,1,1) == "X"){ #Remove X if from column name
    return(substr(name,2,nchar(name)))
  }else{
    return(name)
  }
}

#For each isolate, make 2 plots:
#One including all depth ratios
#One separating histograms by linkage group

for(i in 4:ncol(abs)){ #Loop through individuals
  isolate_name <- get_name(colnames(abs)[i])
  
  #First by linkage group
  pdf(paste("./plots/all_isolates_by_lg/", isolate_name, "_by_lg.pdf", sep=""), height=11.5, width=8)
  
  #Save par settings before switching to 6x3 layout
  old.par <- par(no.readonly = T)
  par(mfrow=c(6,3))
  
  for(j in 1:18){ #Loop through linkage groups
    #Count number of markers
    n <- sum(abs$LG==j & !is.na(abs[,i]), na.rm = T)
    if(n > 0){ #Don't draw any histograms for lgs with 0 markers
      hist(abs[abs$LG==j,i],
           main = paste("LG:", j, "n:", n),
           xlab = paste(isolate_name, "Depth ratio"),
           col = "gray")
      #Draw lines at 0.25, 0.33, 0.5, 0.66, 0.75
      abline(v = 0.333, lty=3)
      abline(v=0.666, lty=3)
      abline(v = 0.25, lty=2)
      abline(v=0.75, lty=2)
      abline(v=0.5, lty=1)
    }
  }
  
  par(old.par) #Go back to old par settings
  dev.off()
  
  #Now do all markers for isolate
  pdf(paste("./plots/all_isolates/", isolate_name, ".pdf", sep=""))
  
  n <- sum(is.na(abs[,i]), na.rm = T)
  
  hist(abs[,i],
       main = paste(isolate_name, "n:", n),
       xlab = "Depth ratio",
       col = "gray")
  abline(v = 0.333, lty=3)
  abline(v=0.666, lty=3)
  abline(v = 0.25, lty=2)
  abline(v=0.75, lty=2)
  abline(v=0.5, lty=1)
  
  dev.off()
  
}
###############################################   ASSIGN PLOIDY LEVEL TO LGs    ######################################

#This function calls even if peak around 0.5, odd if peak around 0.33 or 0.66
get_ploidy <- function(x){
  x <- x[!is.na(x)]
  evens <- x >= 0.47 & x <= 0.53
  odds1 <- x >= 0.30 & x <= 0.36
  odds2 <- x >= 0.63 & x <= 0.69
  if(sum(evens) > sum(odds1) & sum(evens) > sum(odds2)){
    return('2n')
  }else if(sum(odds1) > sum(evens) | sum(odds2) > sum(evens)){
    return('3n')
  }else{
    return(NA)
  }
}

#Given a vector of allele balanaces at heterozygous sites, this function will bootstrap over the sites and call ploidies
bootstrap_ploidies <- function(x, n, alpha, min_snp){
  x <- x[!is.na(x)]
  l <- length(x)
  if(l < min_snp){
    return(NA)
  }else{
    ploidies <- rep(NA, n)
    for(i in 1:n){
      x.boot <- sample(x,l,replace=T)
      ploidies[i] <- get_ploidy(x.boot)
    }
    counts <- table(ploidies)
    counts <- counts[order(counts, decreasing=T)]
    if(counts[1] > (n-(alpha*n))){
      return(names(counts)[1])
    }else{
      return('?')
    }
  }
}

#Call ploidies on samples
n.indvs <- ncol(abs)-3
ploidy_calls <- matrix(NA, nrow=n.indvs, ncol=18)
set.seed(531)
for(i in 1:n.indvs){
  col.i <- i + 3
  for(j in 1:18){
    x <- abs[abs$LG==j, col.i]
    ploidy_calls[i, j] <- bootstrap_ploidies(x, 1000, .05, 30)
  }
}
isolate_names <- sapply(colnames(depths)[4:ncol(depths)], get_name)
rownames(ploidy_calls) <- isolate_names
colnames(ploidy_calls) <- 1:18

###############################################   MAKE BOXPLOTS OF READ DEPTHS    ######################################

#Make plots of read depths by LG
for(i in 4:ncol(depths)){ #Loop through individuals
  isolate_name <- get_name(colnames(depths)[i])
  #First by linkage group
  pdf(paste("./plots/total_depths/", isolate_name, ".pdf", sep=""))
  boxplot(depths[,i] ~ depths$LG)
  dev.off()
}

###################################   VERIFY PLOIDY CALLS USING READ DEPTHS AND INFER BASE PLOIDY LEVELS   ######################################

#Function to pull out comparisons between a given linkage group and the "even" lgs for that given isolate
test_depths <- function(reads, lgs, ploidies, test_lg, downsample = T){
  reads <- reads[!is.na(lgs)]
  lgs <- lgs[!is.na(lgs)]
  lgs.unique <- unique(lgs)
  
  #Downsample to number of snps on lg with least number of snps
  if(downsample==T){
    lg.counts <- table(lgs)
    min.snp <- min(lg.counts)
    reads.down <- rep(NA, min.snp*length(lgs.unique))
    lgs.down <- rep(lgs.unique, each=min.snp)
    for(i in 1:length(lgs.unique)){
      lg <- lgs.unique[i]
      reads.down[((i-1)*min.snp + 1):(i*min.snp)] <- sample(reads[lgs==lg], size=min.snp, replace=F)
    }
    lgs <- lgs.down
    reads <- reads.down
  }
  
  #Fit model and run pairwise comparisons
  reads.lm <- lm(reads~lgs)
  reads.tukey <- TukeyHSD(aov(reads.lm))
  reads.tukey <- as.data.frame(reads.tukey$lgs)
  #Make new columns for lgs being compared
  reads.tukey <- cbind(reads.tukey, t(sapply(rownames(reads.tukey), function(x) unlist(strsplit(x, "-")))))
  colnames(reads.tukey)[5:6] <- c("lg1", "lg2")
  reads.tukey$lg1 <- as.integer(as.character(reads.tukey$lg1))
  reads.tukey$lg2 <- as.integer(as.character(reads.tukey$lg2))
  
  #Pull out comparisons for test isolate and make all comparisons in same direction 
  comp.test <- reads.tukey[reads.tukey$lg1==test_lg | reads.tukey$lg2==test_lg,]
  for(i in 1:nrow(comp.test)){
    if(comp.test$lg1[i]!=test_lg){
      comp.test[i,5:6] <- rev(comp.test[i,5:6])
      comp.test[i,1] <- comp.test[i,1]*-1 #switch direction when switching comparison order
    }
  }
  #determine which are against "even" ploidies
  lg_ploidies <- ploidies[match(comp.test$lg2,1:18)]
  test_to_even <- comp.test[which(lg_ploidies=="2n"),c(1,4,5,6)]
  return(test_to_even)
}

#Loop through and check all the ploidy calls
ploidy_checks <- ploidy_calls
base_ploidies <- rep(NA, nrow(ploidy_calls))
lgs <- as.factor(depths$LG)
set.seed(1207)
for(i in 1:nrow(ploidy_calls)){
  for(j in 1:ncol(ploidy_calls)){
    if(!is.na(ploidy_calls[i,j]) & ploidy_calls[i,j]=="3n"){
      test_isolate <- rownames(ploidy_calls)[i]
      test_reads <- depths[,which(isolate_names==test_isolate) + 3]
      ploidy_compare <- test_depths(test_reads, lgs, ploidy_calls[i,], test_lg = j)
      if(nrow(ploidy_compare) > 0){
        if(any(ploidy_compare$`p adj` > 0.05)){
          ploidy_checks[i,j] <- '?'
        }else{
          if(all(ploidy_compare$diff > 0)){
            base_ploidies[i] <- "diploid"
          }else if(all(ploidy_compare$diff < 0)){
            base_ploidies[i] <- "tetraploid"
          }else{
            base_ploidies[i] <- "?"
          }
        }
      }
    }
  }
}

write.csv(ploidy_checks, "tables/ploidy_calls.csv", quote=F, row.names = T)

###################### Summarize results ####################

#Look at ploidy calls for all isolates
chrom_counts <- matrix(NA, ncol=4, nrow=nrow(ploidy_checks))
colnames(chrom_counts) <- c('NA', '2n', '3n', '?')
rownames(chrom_counts) <- rownames(ploidy_checks)
for(i in 1:nrow(ploidy_checks)){
  chrom_counts[i,1] <- sum(is.na(ploidy_checks[i,]))
  chrom_counts[i,2] <- sum(ploidy_checks[i,] == '2n', na.rm=T)
  chrom_counts[i,3] <- sum(ploidy_checks[i,] == '3n', na.rm=T)
  chrom_counts[i,4] <- sum(ploidy_checks[i,] == '?', na.rm=T)
}

#Look at 'called' linkage groups by isolate
summary(chrom_counts[,2] + chrom_counts[,3])
sum((chrom_counts[,2] + chrom_counts[,3]) == 18)

#Compare number unknown calls to avg read depth
percent_known <- apply(chrom_counts,1, function(x) sum(x[2:3])/18)
avg_depths <- apply(depths[4:ncol(depths)], 2, mean)
cor(percent_known, avg_depths)

chrom_counts <- chrom_counts[order(chrom_counts[,3], chrom_counts[,2], chrom_counts[,4]),]

pdf("plots/ploidy_levels.pdf")
old.par <- par(no.readonly = T)
par(mar=c(8,1,1,1))
barplot(t(chrom_counts), beside=F, col=c('gray', 'blue', 'red'),
        names.arg = rownames(chrom_counts), las=2,
        main = "Chromosomal copy numbers by isolate",
        cex.names = 0.2)
par(old.par)
dev.off()

#Look at ploidy levels by linkage group
chrom_sums <- matrix(NA, ncol=4, nrow=18)
colnames(chrom_sums) <- c('NA', '2n', '3n', '?')
rownames(chrom_sums) <- 1:18
for(i in 1:18){
  chrom_sums[i,1] <- sum(is.na(ploidy_checks[,i]))
  chrom_sums[i,2] <- sum(ploidy_checks[,i] == '2n', na.rm=T)
  chrom_sums[i,3] <- sum(ploidy_checks[,i] == '3n', na.rm=T)
  chrom_sums[i,4] <- sum(ploidy_checks[,i] == '?', na.rm=T)
}
pdf("plots/lg_sums.pdf", height=4, width=6)
par(mar=c(6,4,4,6), xpd=NA)
barplot(t(chrom_sums), beside=F, col=c('gray', 'blue', 'red', 'purple'),
        xlab = "Linkage Group",
        ylab = "Counts",
        cex.names = 0.9, las=2)
legend(22,150, 
       c("NA", "2n", "3n", "?"), 
       fill=c('gray', 'blue', 'red', 'purple'),
       bty='n')
par(old.par)
dev.off()

chrom_sums.q <- t(chrom_sums[,-1])
chisq.test(chrom_sums.q)
apply(chrom_sums.q,2,function(x) x[1]/sum(x))

#Look at aneuploid levels
number_odds <- rep(0, 19)
names(number_odds) <- 0:18
odd_counts <- table(chrom_counts[,3])
number_odds[match(names(odd_counts), names(number_odds))] <- odd_counts
pdf("plots/aneuploid_distribution.pdf", height=4, width=6)
barplot(odd_counts, 
        xlab = "Number triploid chromosomes",
        ylab = "Isolate counts")
dev.off()

## Ploidy variation within clonal groups
rotate <- function(x) t(apply(x,2,rev))
for(ug in 1:41){
  clones <- pcap$SampleSZ[pcap$UniqueGenotype == ug]
  ploidy_checks.ug <- ploidy_checks[rownames(ploidy_checks) %in% clones,]
  ploidy_checks.ug[ploidy_checks.ug == "2n"] <- 0
  ploidy_checks.ug[ploidy_checks.ug == "3n"] <- 1
  ploidy_checks.ug[ploidy_checks.ug == "?"] <- NA
  class(ploidy_checks.ug) <- "numeric"
  pdf(paste("plots/clonal_lineages/ug", ug, ".pdf", sep=""))
  image(rotate(ploidy_checks.ug), col=c("blue", "red"), xaxt='n', yaxt='n',
        main = ug, zlim=c(0,1))
  dev.off()
}

#Enrichment for aneuploidy in certain clonal groups compared to others?
ug_sums <- matrix(NA, nrow=41, ncol=2)
for(ug in 1:41){
  clones <- pcap$SampleSZ[pcap$UniqueGenotype == ug]
  ploidy_checks.ug <- ploidy_checks[rownames(ploidy_checks) %in% clones,]
  aneuploid_counts <- sum(apply(ploidy_checks.ug,1,function(x) any(x=="3n", na.rm = T) & any(x=="2n", na.rm = T)))
  euploid_counts <- nrow(ploidy_checks.ug) - aneuploid_counts
  ug_sums[ug,] <- c(aneuploid_counts, euploid_counts)
}
fisher.test(ug_sums)

#Make plot showing example featuring isolate with clear patterns
pdf("plots/ploidy_fig_13EH38A.pdf", height=5.5, width=7)

#Choose colors representing NA, even, and odd
ploidy_types <- c(NA, "2n", "3n")
ploidy_colors <- c("gray", col=brewer.pal(3, "RdBu")[c(1,3)])

#Set up the plot layout
plot_layout <- rbind(c(1:6),
                     c(7:12),
                     c(13:18),
                     c(19),
                     c(19),
                     c(19))
layout(plot_layout)

#Par settings for plot 1: Ploidy example
old.par <- par(no.readonly = T)
par(mar=c(2,1,1,0), oma=c(4,6,3,4), xpd=NA)

#Plot 1: Ploidy histogram for a representative isolate
isolate <- "13EH38A"
calls <- ploidy_checks[isolate,]
depth_column <- which(rownames(ploidy_checks) == isolate) + 3
for(i in 1:18){
  n <- sum(abs$LG==i & !is.na(abs[,depth_column]), na.rm = T)
  ylabel <- ""
  xlabel <- ""
  xaxt <- "n"
  yaxt <- "n"
  if(i %in% c(1,7,13)){
    yaxt <- "s"
  }
  if(i %in% c(13:18)){
    xaxt <- "s"
  }
  if(i==7){
    ylabel <- "Density"
  }
  if(i==15){
    xlabel <- "                     Allele balance"
  }
  hist(abs[abs$LG==i,depth_column], 
       col = ploidy_colors[match(calls[i], ploidy_types)],
       main = bquote(paste(bold("LG "), bold(.(as.character(i))),
                                  bold(" ("),bolditalic("n"),bold("="),
                                  bold(.(as.character(n))),bold(")"), sep='')),
       xlab = xlabel,
       ylab = ylabel,
       cex.main=1,
       cex.lab=1.5,
       cex.axis=1.1,
       xlim=c(0,1),
       ylim=c(0,4),
       xaxt=xaxt,
       yaxt=yaxt,
       freq=F)
  axis(side=1,labels=F)
  axis(side=2,labels=F)
}

#Par settings for plot 2
par(mar=c(4,1,4,0), xpd=NA)

#Plot 2: total read depths for representative isolate
boxplot(depths[,depth_column]~depths$LG, 
        outline=F, 
        col=ploidy_colors[match(calls, ploidy_types)],
        xlab = "Linkage group",
        ylab = "Read depth",
        cex.lab=1.5,
        cex.axis=1.1)

legend(5.1,-38, 
       c("NA", expression(paste("2", italic("n"))), expression(paste("3", italic("n")))), 
       fill=ploidy_colors,
       x.intersp = 0.4,
       text.width=1.5,
       ncol=3,
       bty='n', cex=1.7)

text(c(-2,-2), c(260, 85), c("A", "B"), cex=2)
par(old.par)
dev.off()

#Make plot showing example with one isolate ----- unclear patterns
pdf("plots/ploidy_fig_12889MIA.pdf", height=5.5, width=7)

#Choose colors representing NA, even, and odd
ploidy_types <- c('?', "2n", "3n")
ploidy_colors <- c("gray", col=brewer.pal(3, "RdBu")[c(1,3)])

#Set up the plot layout
plot_layout <- rbind(c(1:6),
                     c(7:12),
                     c(13:18),
                     c(19),
                     c(19),
                     c(19))
layout(plot_layout)

#Par settings for plot 1: Ploidy example
old.par <- par(no.readonly = T)
par(mar=c(2,1,1,0), oma=c(4,6,3,4), xpd=NA)

#Plot 1: Ploidy histogram for a representative isolate
isolate <- "12889MIA"
calls <- ploidy_checks[isolate,]
depth_column <- which(rownames(ploidy_checks) == isolate) + 3
for(i in 1:18){
  n <- sum(abs$LG==i & !is.na(abs[,depth_column]), na.rm = T)
  ylabel <- ""
  xlabel <- ""
  xaxt <- "n"
  yaxt <- "n"
  if(i %in% c(1,7,13)){
    yaxt <- "s"
  }
  if(i %in% c(13:18)){
    xaxt <- "s"
  }
  if(i==7){
    ylabel <- "Density"
  }
  if(i==15){
    xlabel <- "                     Allele balance"
  }
  hist(abs[abs$LG==i,depth_column], 
       col = ploidy_colors[match(calls[i], ploidy_types)],
       main = bquote(paste(bold("LG "), bold(.(as.character(i))),
                                  bold(" ("),bolditalic("n"),bold("="),
                                  bold(.(as.character(n))),bold(")"), sep='')),
       xlab = xlabel,
       ylab = ylabel,
       cex.main=1,
       cex.lab=1.5,
       cex.axis=1.1,
       xlim=c(0,1),
       ylim=c(0,4),
       xaxt=xaxt,
       yaxt=yaxt,
       freq=F)
  axis(side=1,labels=F)
  axis(side=2,labels=F)
}

#Par settings for plot 2
par(mar=c(4,1,4,0), xpd=NA)

#Plot 2: total read depths for representative isolate
boxplot(depths[,depth_column]~depths$LG, 
        outline=F, 
        col=ploidy_colors[match(calls, ploidy_types)],
        xlab = "Linkage group",
        ylab = "Read depth",
        cex.lab=1.5,
        cex.axis=1.1)

legend(5.1,-15, 
       c("?", expression(paste("2", italic("n"))), expression(paste("3", italic("n")))), 
       fill=ploidy_colors,
       x.intersp = 0.4,
       text.width=1.5,
       ncol=3,
       bty='n', cex=1.7)

text(c(-2,-2), c(105, 34), c("A", "B"), cex=2)
par(old.par)
dev.off()

########################################## Random calculations ################################
sum(!is.na(ploidy_checks)) #Non NA 
sum(ploidy_checks!="2n",na.rm=T) #Odd
odd_counts

aneuploid <- apply(ploidy_checks,1, function(x) any(x=="3n", na.rm=T))
sum(aneuploid)
odd_counts
pcap <- pcap[match(rownames(ploidy_checks), pcap$SampleSZ),]
aggregate(aneuploid~pcap$UniqueGenotype, FUN=sum)

#Fisher test for enrichment of triploid
chrom_test <- chrom_sums[,c("2n", "3n")]
fisher.test(chrom_test, simulate.p.value=T, B=100000)

########### Let's make chrom sums table
chrom_table <- chrom_sums[,c('2n', '3n', '?', 'NA')]
chrom_table <- as.data.frame(chrom_table)
chrom_table$Percent3n <- round(chrom_table$`3n`/(chrom_table$`2n` + chrom_table$`3n`),2)
chrom_table <- rbind(chrom_table, apply(chrom_table,2,sum))
chrom_table <- cbind("LG"=c(1:18, "Sum"), chrom_table)
chrom_table[nrow(chrom_table),6] <- NA
write.csv(chrom_table, "tables/lg_sums.csv", row.names = F, quote = F)


### quickly look at median read depths both across individuals and sites
depths.filt <- depths[,-c(1:3)]
depths.filt[depths.filt < 5] <- NA #All sites with fewer than 5 reads were set to NA in VCF file
site_means <- apply(depths.filt, 1, mean, na.rm=T)
sample_means <- apply(depths.filt, 2, mean, na.rm=T)
summary(site_means)
summary(sample_means)
