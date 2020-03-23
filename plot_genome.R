plot_stat <- function(
  stats, #data frame with columns for CHROM, POS, and plotting stats
  chrom = 1, #column corresponding to CHROM
  pos = 2, #column corresponding to POS,
  plotting_stats = 3:ncol(stats), #columns corresponding to data to plot
  scaffold_sizes, #First column scaffold name, second column size of scaffold
  
  #Fixed plotting arguments
  plot_base_name = "sliding_window",
  number_frames = 1,
  frame_stretch = 0,
  width_to_height = 3,
  horizontal_lines = c(0.5),
  stat_colors = rep("black", length(plotting_stats)), #Vector of colors to use
  stat_ltys = rep(1, length(plotting_stats)),
  stat_lwds = rep(1, length(plotting_stats)),
  stat_types = rep('l', length(plotting_stats)),
  tick_spacing = 5000000,
  #Passed arguments to plot, lines, and legend
  plot_arguments = list(),
  chrom_tick_arguments = list(),
  lines_arguments = list(),
  legend_arguments = list(),
  chrom_arguments = list(),
  horizontal_line_arguments = list(),
  tick_arguments = list()
){
  
  ##########    FIRST SET UP DATA FOR PLOTTING    ##########
  
  #Get information out of scaffold sizes matrix
  lgs.unique <- scaffold_sizes[,1]
  lgs.length <- scaffold_sizes[,2]
  
  #Make cumulative positions
  lgs <- stats[,chrom]
  lgs.cumsum <- c(0,cumsum(scaffold_sizes[,2]))
  cumpos <- rep(NA, nrow(stats))
  for(i in 1:length(lgs.unique)){
    scaff <- lgs.unique[i]
    add <- lgs.cumsum[i]
    cumpos[stats[,chrom]==scaff] <- stats[,pos][stats[,chrom]==scaff] + add
  }
  
  #Make coordinates of 5 MB intervals for tick mark positions 
  tick_names <- c()
  tick_cumpos <- c()
  for(i in 1:length(lgs.unique)){
    scaff <- lgs.unique[i]
    pos.lg <- stats[,pos][stats[,chrom]==scaff]
    lg.tick_names <- seq(0,max(pos.lg), by = tick_spacing)[-1]
    lg.tick_cumpos <- lg.tick_names + lgs.cumsum[i]
    tick_names <- c(tick_names, lg.tick_names)
    tick_cumpos <- c(tick_cumpos, lg.tick_cumpos)
  }
  
  #Get which chromosomes go in each frame
  cumsum.max <- max(lgs.cumsum)
  frame.length <- cumsum.max/number_frames
  frames.lgs <- sapply(1:number_frames,function(x) return(sum(lgs.cumsum/frame.length < x)))
  frames.lgstarts <- c(1, frames.lgs[-number_frames] + 1)
  #Separate x axis positions for separate frames
  frame.starts <- lgs.cumsum[c(1, frames.lgs[-length(frames.lgs)] + 1)] + 1
  frame.ends <- lgs.cumsum[frames.lgs + 1]
  
  #Get relative vertical positions of frames in windows
  frame.vsize <- 1/number_frames
  frame.vstops <- seq(1,0, by=-frame.vsize)[-(number_frames+1)] + frame_stretch                                               
  frame.vstarts <- frame.vstops - frame.vsize - 2*frame_stretch
  frame.vstops[1] <- 1
  frame.vstarts[number_frames] <- 0
  
  #Make matrix with coordinates of chromosome centers for tick marks
  lgs.centers <- matrix(NA, nrow=length(lgs.unique), ncol=2)
  for(i in 1:nrow(lgs.centers)){
    lgs.centers[i,1] <- as.character(lgs.unique)[i]
    lgs.centers[i,2] <- (lgs.cumsum[i] + lgs.cumsum[i+1])/2
  }
  
  ##########    NOW START ADDING PLOTTING ELEMENTS    ##########
  
  #Make plot
  pdf(paste(plot_base_name, ".pdf", sep=""), height=6, width=6*width_to_height)
  old.par <- par(no.readonly = T)
  par(mar=c(4.5,6,4,4))
  
  for(f in 1:number_frames){ 
    par(fig=c(0,1,frame.vstarts[f], frame.vstops[f]), new=T)
    lgs.window <- lgs.unique[frames.lgstarts[f]:frames.lgs[f]]
    #Set up empty plot
    total_plot_arguments <- list(x=0, type="n", xaxt="n",
                                 xlim=c(frame.starts[f], frame.ends[f]),
                                 ylim = c(min(stats[lgs %in% lgs.window, plotting_stats]),
                                          max(stats[lgs %in% lgs.window, plotting_stats])),
                                 xlab = "Genomic position",
                                 cex.lab=1.3,
                                 cex.main = 1.5)
    total_plot_arguments[names(plot_arguments)] <- plot_arguments
    
    #Suppress x label on top frames and main label on lower frames
    if(f < number_frames){
      total_plot_arguments["xlab"] <- ""  
    }else if(f > 1){
      total_plot_arguments["main"] <- NULL
    }
    do.call(plot, total_plot_arguments)
    
    
    #Add axis tick marks for chromosome centers
    total_chrom_tick_arguments <- list(side = 1, at = lgs.centers[,2],
                                       labels = lgs.centers[,1],
                                       tick = F, line=0.5, cex.axis=1)
    total_chrom_tick_arguments[names(chrom_tick_arguments)] <- chrom_tick_arguments
    do.call(axis, total_chrom_tick_arguments)
    #Add axis tick marks for position on chromosome
    total_tick_arguments <- list(side = 1, at = tick_cumpos,
                                 labels = paste(tick_names/1000000, "Mb"),
                                 tick = T, cex.axis=0.6, mgp=c(3,0.25,0))
    total_tick_arguments[names(tick_arguments)] <- tick_arguments
    do.call(axis, total_tick_arguments)
    
    #Add horizontal lines
    for(h in 1:length(horizontal_lines)){
      total_horizontal_line_arguments = list(h = horizontal_lines[h], lty=2, col="snow4", lwd=0.5)
      total_horizontal_line_arguments[names(horizontal_line_arguments)] <- horizontal_line_arguments
      do.call(abline, total_horizontal_line_arguments)
    }
    
    #Loop through stats columns
    counter <- 0
    for(i in plotting_stats){
      counter <- counter + 1
      #Loop through chromosomes
      for(lg in lgs.window){
        y <- stats[lgs==lg, i]
        x <- cumpos[lgs==lg]
        x <- x[!is.na(y)]
        y <- stats[lgs==lg, i]
        y <- y[!is.na(y)]
        total_lines_arguments <- list(x=x,
                                      y=y,
                                      col = stat_colors[counter],
                                      lty = stat_ltys[counter],
                                      lwd = stat_lwds[counter],
                                      type = stat_types[counter])
        #total_lines_arguments[names(lines_arguments)] <- lines_arguments
        do.call(lines, total_lines_arguments)
      }
    }
    
    #Separate chromosomes with vertical lines
    for(i in 1:length(lgs.cumsum)){
      total_chrom_arguments <- list(v=lgs.cumsum[i], col="gray", lty=2, lwd=1.5)
      # total_chrom_arguments[names(chrom_arguments)] <- chrom_arguments
      do.call(abline, total_chrom_arguments)
    }
    
  }
  
  #Add legend
  total_legend_arguments = list(x = "bottomleft", legend = "", bty = "n")
  total_legend_arguments[names(legend_arguments)] <- legend_arguments
  do.call(legend, total_legend_arguments)
  
  par(old.par)
  
  dev.off()
  
}

