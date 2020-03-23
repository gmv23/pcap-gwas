# Clone correction of P capsici isolates from around NY and US
# This will return file with clonal group designations of all isolates
# And then sample the isolate from each clonal group with the least missing data
# Then it will return a list of the samples corresponding the cc dataset

library(igraph)

################################# READ IBS MATRIX, GENOTYPE FILES, and MISSING DATA FILE ##############################

#Load IBS matrix and make individuals row and column names
IBS_matrix <- read.csv("IBS_matrix.csv", header = T)
row.names(IBS_matrix) <- IBS_matrix$X
IBS_matrix$X <- NULL
IBS_matrix <- as.matrix(IBS_matrix)

#Load individual info
indvs <- read.table("capsici_diversity_FILTER.012.indv", stringsAsFactors = F)
indvs <-unlist(indvs$V1)

#Load missing data info
miss <- read.table("ind_missing.imiss", header=T)

#Load phenos data to see which samples are alive in storage
phenos <- read.csv("/workdir/gmv23/pcap_diversity/phenos/Pcap_isolates.csv", stringsAsFactors = FALSE)

############################### ASSIGN INDIVIDUALS TO CLONAL GROUPS ################################

#Turn high IBS cells of matrix to 1
modify_matrix <- function(x){
  if(is.na(x) | x<.95){
    return(0)
  }else{
    
    return(1)
  }
}
clone_or_not <- structure(sapply(IBS_matrix, modify_matrix), dim=dim(IBS_matrix))

# Create network -> Each isolate is a node and there is an edge between isolates that are clones
g <- graph_from_adjacency_matrix(clone_or_not, "undirected")

# Clusters are isolates that only have edges between themselves and not the rest of the network (ie clones)
g.clusters <- clusters(graph = g)

#### Create table of clonal group assignments

# Make list of cluster size corresponding to each member of network (used later)
cluster_sizes <- rep(NA, length(indvs))
for(i in 1:length(cluster_sizes)){
  member <- g.clusters$membership[i]
  size <- sum(g.clusters$membership == member)
  cluster_sizes[i] <- size
}

# Prepare table and variables for loop
clonal_groups <- 1:(g.clusters$no)
clone_assignments <- matrix(ncol=2)
colnames(clone_assignments) <- c("Sample", "Clonal_group")
counter <- 0

# Assign individuals to clonal groups starting with largest group
for(i in 1:length(unique(g.clusters$csize))){ #loop through all unique cluster sizes
  # Start with largest cluster size
  current_size <- sort(unique(g.clusters$csize), decreasing=T)[i] 
  # how many groups of this size are there
  same_size_clonal_groups <- unique(g.clusters$membership[cluster_sizes == current_size]) 
  #loop through groups of that size
  for(j in 1:length(same_size_clonal_groups)){ 
    counter <- counter +1
    old_clonal_group_id <- same_size_clonal_groups[j] #Assignment to group from g.clusters$membership
    new_clonal_group_assignment <- clonal_groups[counter] #New assignment going from largest to smallest
    clone_assignments <- rbind(clone_assignments, cbind(
      indvs[which(g.clusters$membership == old_clonal_group_id)],
      new_clonal_group_assignment))
  }
}
clone_assignments <- clone_assignments[-1,]
clone_assignments <- as.data.frame(clone_assignments, stringsAsFactors = F)
clone_assignments$Clonal_group <- as.integer(clone_assignments$Clonal_group)

write.table(clone_assignments, "Clone_assignments.txt", row.names = F, quote=F, sep="\t")

########### Clone correct -- choosing individual with least missing data per clonal group #########
########### and ensuring individual is still alive in storage if possible #########

#Get just part of name that will match with pheno file
shorten_name <- function(x){
  return(unlist(strsplit(x,":"))[1])
}

pheno_names <- paste(phenos$Sample, phenos$SZ_for_analysis, sep="")

indvs_cc <- rep(NA, max(clone_assignments$Clonal_group))
cc_samples_list <- matrix(rep(NA, 2*max(clone_assignments$Clonal_group)), ncol=2)
colnames(cc_samples_list) <- c("Sample", "In_storage")
for(i in 1:max(clone_assignments$Clonal_group)){
  samples.i <- clone_assignments$Sample[clone_assignments$Clonal_group == i]
  samples.short.i <- sapply(samples.i, shorten_name)
  samples.exist <- as.logical(phenos$In_storage[match(samples.short.i, pheno_names)])
  if(any(samples.exist)){
    samples.i <- samples.i[samples.exist]
    samples.short.i <- samples.short.i[samples.exist]
    cc_samples_list[i,2] <- 1
  }else{
    cc_samples_list[i,2] <- 0
  }
  samples.missing <- miss$F_MISS[match(samples.i, miss$INDV)]
  best_sample <- which(samples.missing == min(samples.missing))
  indvs_cc[i] <- samples.i[best_sample]
  cc_samples_list[i,1] <- samples.short.i[best_sample]
}

write.table(indvs_cc, "capsici_diversity_CC.012.indv", quote = F, row.names = F, col.names = F)
write.csv(cc_samples_list, "capsici_diversity_cc_samplelist.csv", quote=F, row.names = F)

