# This script is to summarize trends from clone correction of P. capsici samples
# Output:
     # A table summarizing number of isolates and clones found among different fields
     # A map of NY showing sample locations
     # Analysis of phenotypic consistency by clonal lineage

################################################ Import data and packages #############################################

setwd("~/Documents/work/Smart_lab/P_capsici/isolate_collection/paper/phenotypes_and_clones/")

library(RColorBrewer)
library(sp)
library(raster)
library(lme4)

pcap <- read.csv("tables/isolate_metadata.csv", na.strings = "NA", stringsAsFactors = F)
clones <- read.table("Clone_assignments.txt", header=T, stringsAsFactors = F)
indvs.cc <- read.table("capsici_diversity_CC.012.indv", stringsAsFactors = F)
indvs.cc <- unlist(indvs.cc$V1)

#Clean up all sample names --- get rid of extra fields output by TASSEL
shorten_name <- function(x){
  return(unlist(strsplit(x,":"))[1])
}
clones$Sample <- sapply(clones$Sample, shorten_name)
indvs.cc <-sapply(indvs.cc, shorten_name)

################################################ Isolate number table and plotting metadata #############################################

#Make table assigning counties to geographical regions of NY
regions <- rbind(c("Cayuga", "CNY"),
                 c("Columbia", "CD"),
                 c("Erie", "WNY"),
                 c("Herkimer", "CD"),
                 c("Monroe", "CNY"),
                 c("Ontario", "CNY"),
                 c("Rensselaer", "CD"),
                 c("Schenectady", "CD"),
                 c("Suffolk", "LI"),
                 c("Tioga", "CNY"),
                 c("Tompkins", "CNY"))
regions <- as.data.frame(regions)
colnames(regions) <- c("County", "Region")

#Get 'populations' designation: combination of year and field
populations <- paste(pcap$Field, " (", pcap$year, ")", sep="")
populations[populations == "Other NY (NA)"] <- "NY unknown"
populations[populations == "NonNY (NA)"] <- pcap$State[populations == "NonNY (NA)"]

#Now make table assigning fields to counties, regions, and plotting colors/symbols
field_info <- data.frame("Field" = unique(populations),
                         "County" = NA,
                         "Region" = NA)
for(i in 1:nrow(field_info)){
  field <- as.character(field_info$Field[i])
  county <- NA
  if(field == "NY unknown"){
    region <- "NY unknown"
  }else if(field %in% c("OH", "MI", "SC", "FL", "CA", "NM")){
    region <- "Non NY"
  }else{
    county <- unlist(strsplit(field," #", fixed=T))[1]
    region <- as.character(regions$Region)[regions$County==county]
  }
  field_info$Region[i] <- region
  field_info$County[i] <- county
}

#Put in alphabetical order with non NY stuff at end
region_order <- c("CD", "CNY", "LI", "WNY", "NY unknown", "Non NY")
field_unordered <- field_info
start <- 1
for(region in region_order){
  region.size <- sum(field_unordered$Region==region)
  field_info[start:(start+region.size-1),] <- field_unordered[field_unordered$Region==region,]
  start <- start + region.size
}

#Add color and pch
field_info$color <- NA
color_types <- brewer.pal(7,"Set1")[-6]
field_info$color <- color_types[match(field_info$Region, unique(field_info$Region))]
#Order sites within fields and add pch info
field_info$pch <- NA
for(region in unique(field_info$Region)){
  n.sites <- sum(field_info$Region==region)
  field_info.region <- field_info[field_info$Region==region,]
  field_info[field_info$Region==region,] <- field_info.region[order(field_info.region$Field),]
  field_info$pch[field_info$Region==region] <- 1:n.sites
}

#Save table assigning isolates to fields, counties, regions, and plotting info
isolate_metadata <- cbind("SampleSZ" = pcap$SampleSZ, 
                          field_info[match(populations, field_info$Field),])
write.csv(isolate_metadata, "isolate_plotting_metadata.csv", quote=F, row.names = F)
write.csv(field_info, "field_plotting_metadata.csv", quote=F, row.names = F)

#Get total isolates per field-year and total unique genotypes per field-year
populations[populations %in% c("NM", "CA", "SC", "MI", "FL", "OH")] <- "Non NY"
population_counts <- table(populations)
field_sums <- data.frame("Field" = gsub("(.*) \\(.*", "\\1", rownames(population_counts), perl=T),
                         "Year" = as.integer(gsub(".*\\((.*)\\)", "\\1", rownames(population_counts), perl=T)),
                         "NumberIsolates" = as.integer(population_counts))
field_counties <- sapply(as.character(field_sums$Field), function(x) unlist(strsplit(x, " "))[1])
field_sums$Region <- regions$Region[match(field_counties, regions$County)]

#List MTs and Mefs for all isolates 
MT_counts <- as.data.frame.matrix(table(populations, pcap[,c("MT")]))
Mef_counts <- as.data.frame.matrix(table(populations, pcap[,c("Mefenoxam_sensitivity")]))
populationsCC <- populations[pcap$InCCDataset]
populationCC_counts <- table(populationsCC)
names(population_counts) == names(populationCC_counts)

field_sums <- data.frame(field_sums[c("Field", "Region", "Year", "NumberIsolates")], 
                         MT_counts[,c("A1", "A2")], 
                         Mef_counts[,c("S", "IS", "R")])
field_sums$NumberUniqueGenotypes <- as.integer(populationCC_counts)
field_sums <- field_sums[order(field_sums$Region, field_sums$Field),]

#Get totals
isolate_totals <- apply(field_sums[,4:ncol(field_sums)], 2, sum)
field_sums$Field <- as.character(field_sums$Field)
field_sums <- rbind(field_sums, c("Total", NA, NA, isolate_totals))

write.csv(field_sums, "tables/isolate_sums.csv", quote=F, row.names = F)

################################################### Sample map ######################################################

#Note this includes all isolates including those that didnt pass genotyping

usa <- getData('GADM', country='US', level=2)
ny <- usa[usa$NAME_1=="New York",]
ny.clean <- ny[ny@data$TYPE_2=='County',]
counties.unique <- unique(pcap$County)
counties.unique <- counties.unique[!is.na(counties.unique)]
pdf("plots/fieldmap.pdf", width = 3.25, height= 2.141176)
old.par <- par(no.readonly = T)
par(ann = FALSE,
    bg = "white",
    bty = "n",
    mai = c(0,0,0,0),
    mgp = c(0, 0, 0),
    oma = c(0,0,0,0),
    omd = c(1,1,1,1),
    omi = c(0,0,0,0),
    usr = c(-79.5, -71, 40, 45.6),
    pin = c(3.25,2.141176),
    plt = c(0,1,0,1),
    pty = "m",
    xaxs = 'i',
    xaxt = 'n',
    xpd = FALSE,
    yaxs = 'i',
    yaxt = 'n')
plot(ny.clean, xlim = c(-79.5, -71), ylim=c(40,45.6))
set.seed(823)
for(county in counties.unique){
  county.polygon <- ny[ny$NAME_2==county,]
  county.color <- unique(field_info$color[field_info$County==county])
  county.color <- county.color[!is.na(county.color)]
  for(field in unique(pcap$Field[pcap$County == county & !is.na(pcap$County)])){
    if(field != "Other NY"){
      field.point <- spsample(county.polygon,1,type='random')
      for(year in unique(pcap$year[pcap$Field == field])){
        n.isolates <- sum(pcap$year == year & pcap$Field == field, na.rm = T)
        if(year < 2017){
          point.pch = 16
        }else{
          point.pch = 17
        }
        points(field.point, cex=((n.isolates-.95)^.20)+0.4, 
               col = 'black',
               pch = point.pch)
        points(field.point, cex=(n.isolates-.95)^.20, 
               col = county.color,
               pch = point.pch)
      }
    }
  }
}

legend(-72.8, 45.35, 
       pt.cex=1.5,
       pch = c(16,17),
       legend = c("Prior to 2017", "2017-2018"),
       cex=0.7, bty='n', y.intersp=1.5)
legend(-72.8, 44.2, 
       legend = unique(field_info$Region)[1:4], 
       fill = unique(field_info$color)[1:4],
       cex=0.7, bty='n', y.intersp=1.5)
par(old.par)
dev.off()

############################################# Random calculations ###############################################

#How many clonal lineages are there?
ug_sizes <- table(pcap$UniqueGenotype)
clineages <- as.integer(names(ug_sizes[ug_sizes>1]))
length(clineages)
cl_sizes <- ug_sizes[names(ug_sizes) %in% clineages]
sum(cl_sizes)
max(cl_sizes)
median(cl_sizes)

max_cl <- max(clineages)

#Number of each MT
table(pcap$MT[pcap$State == "NY"])

#Mefenoxam by location
mef_table <- table(pcap$Field[pcap$State == "NY"], pcap$Mefenoxam_sensitivity[pcap$State == "NY"])

#MT by clonal lineage and how many discordant
MT_variation <- table(pcap$UniqueGenotype, pcap$MT)
which(apply(MT_variation, 1, function(x) sum(x!=0))>1)
MT_variation[16,]

#Mef by clonal lineage
mef_variation <- table(pcap$UniqueGenotype, pcap$Mefenoxam_sensitivity)
discordant_mef <- mef_variation[apply(mef_variation, 1, function(x) sum(x!=0))>1,]
discordant_mef

#% variation in mef due to differences between clonal lineage
clonal <- pcap[pcap$UniqueGenotype %in% clineages & !is.na(pcap$Mef5) & !is.na(pcap$Mef100),]
clonal$UniqueGenotype <- as.factor(clonal$UniqueGenotype)

#Mef5
y <- clonal$Mef5
X <- cbind(rep(1,length(y)))
Z <- model.matrix(~0+clonal$UniqueGenotype)
mef5.mm <- mixed.solve(y=y, X=X, Z=Z)
mef5.mm$Vu/(mef5.mm$Vu+mef5.mm$Ve)

#Mef100
y <- clonal$Mef100
mef100.mm <- mixed.solve(y=y, X=X, Z=Z)
mef100.mm$Vu/(mef100.mm$Vu+mef100.mm$Ve)

mef5.lmer <- lmer(clonal$Mef5~(1|clonal$UniqueGenotype))
mef5.var <- as.data.frame(VarCorr(mef5.lmer))
mef5.var$vcov[1]/sum(mef5.var$vcov)
mef100.lmer <- lmer(clonal$Mef100~(1|clonal$UniqueGenotype))
mef100.var <- as.data.frame(VarCorr(mef100.lmer))
mef100.var$vcov[1]/sum(mef100.var$vcov)

#Clone corrected
pcap.cc <- pcap[pcap$InCCDataset,]
table(pcap.cc$MT)
table(pcap.cc$Mefenoxam_sensitivity)

#Testing for deviations from 1:1
pcap.cc$Fieldyear <- paste(pcap.cc$Field, pcap.cc$year)
MT_sums <- table(pcap.cc$Fieldyear, pcap.cc$MT)
apply(MT_sums, 1, binom.test, alternative = "two.sided")
Mef_sums <- table(pcap.cc$Fieldyear, pcap.cc$Mefenoxam_sensitivity)



