#Sliding window plotting

window_smooth <- function(
  stats, #data frame with columns for CHROM, POS, and stats to smooth
  chrom = 1, #column corresponding to CHROM
  pos = 2, #column corresponding to POS
  data_columns = 3:ncol(stats), #columns corresponding to data to plot
  smooth_stat = 2, #column corresponding to stat to use to divide windows (ie physical position, number SNPs, etc)
  smooth_by = c("distance", "total", "sum"), #Windows based on distance or number of X statistic per window
  FUN = mean, #function to apply on data in window
  window_size = 1000000, 
  step_size = 100000,
  scaffold_sizes = NULL #First column scaffold name, second column size of scaffold only for smooth_by=distance
){
  
  FUN <- match.fun(FUN)
  out <- data.frame()
  lgs <- stats[,chrom]  
  lgs <- droplevels(lgs) #Get rid of unused levels
  lgs.unique <- unique(lgs)
  
  ############ First Get 'pseudo-position' vector depending on type of smooth_by #############

  if(smooth_by == "distance"){
    if(is.null(scaffold_sizes)){
      stop("Need scaffold_sizes matrix if smooth_by == distance")
    }
    pseudopos <- stats[,pos] 
    
  }else if(smooth_by == "sum" | smooth_by == "total"){
    pseudopos <- rep(NA, nrow(stats))
    lg.snpsums <- as.vector(table(lgs))
    lg.snpcumsums <- cumsum(lg.snpsums)
    lg.snpstarts <- c(1, (lg.snpcumsums+1)[-length(lg.snpcumsums)])
    for(i in 1:length(lgs.unique)){
      
      if(smooth_by == "sum"){
        pseudopos[lg.snpstarts[i]:lg.snpcumsums[i]] <- 1:lg.snpsums[i]
        
      }else{ #smooth_by == "total"
        lg.smoothstats <- stats[lg.snpstarts[i]:lg.snpcumsums[i], smooth_stat]
        pseudopos[lg.snpstarts[i]:lg.snpcumsums[i]] <- cumsum(as.numeric(lg.smoothstats))
        
      }
    }
  }else{
      stop("Need a valid smooth_by type")
  }

############ Now start applying smoothing function in windows ###########
   
  # Do each chromosome separately
  for(lg in lgs.unique){
    stats.lg <- stats[lgs == lg,]
    pseudopos.lg <- pseudopos[lgs==lg]
    pos.lg <- stats.lg[,pos]
    
    if(smooth_by == "distance"){
      pos.max <- scaffold_sizes[which(scaffold_sizes[,1]==lg),2]
    }else{
      pos.max <- max(pseudopos.lg)
    }
    
    #Get starts and stops of windows
    start_positions <- seq(0, step_size*ceiling(pos.max/step_size)-window_size, by = step_size)
    if(pos.max %% step_size == 0){
      stop_positions <-  seq(window_size, pos.max , by = step_size)
    }else{
      stop_positions <-  c(seq(window_size, pos.max , by = step_size), pos.max)
    }
    
    #Get window centers based on physical position
    centers.n <- length(start_positions)
    if(smooth_by == "distance"){
      centers <- (stop_positions + start_positions)/2
    }else{
      centers <- rep(NA, centers.n)
    }
    smoothed_data <- matrix(NA, nrow=centers.n, ncol=length(data_columns))
    
    for(i in 1:centers.n){
      window_data <- stats.lg[pseudopos.lg > start_positions[i] & pseudopos.lg <= stop_positions[i], data_columns]
      pos.window <- pos.lg[pseudopos.lg > start_positions[i] & pseudopos.lg <= stop_positions[i]]
      if(smooth_by != "distance"){
        centers[i] <- mean(pos.window) #Take mean of SNPs contributing to window if window not based on physical position
      }
      if(class(window_data) == "data.frame" | class(window_data) == "matrix"){
        smoothed_data[i,] <- apply(window_data,2,FUN)
      }else{
        smoothed_data[i,] <- FUN(window_data)
      }
    }
    chrom_results <- data.frame("CHROM" = rep(lg, centers.n),
                                "POS" = centers,
                                smoothed_data)
    colnames(chrom_results)[3:ncol(chrom_results)] <- colnames(stats)[data_columns]
    out <- rbind(out, chrom_results)
  }
  return(out)
}

