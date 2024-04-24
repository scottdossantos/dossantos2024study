# This function is a truncated version of the aldex.effect() function from the
# ALDEx2 package, maintained by Dr. Gregory Gloor. Source code for aldex.effect
# available at: https://github.com/ggloor/ALDEx_bioc/blob/master/R/clr_effect.r

# This version of aldex.effect is used in 'preg_vs_nonpreg.R' to calculate the 
# median CLR values for 31,654 genes (across 128 Monte-Carlo instances) in 341
# vaginal metatranscriptomes from the corresponding aldex.clr object, without
# also calculating within & between group median CLR differences and effect
# sizes for all genes. Using the complete aldex.effect() function on this data
# causes RStudio to throw a fatal error (likely due to memory exhaustion) on a
# 2023 MacBook Pro (M2 chip, 16 GB RAM). This truncated function was used as a
# workaround. As a result, arguments 'CI', 'glm.conds' and 'paired.test' are 
# now functionally irrelevant and DO NOT change the output of the function.

trunc.aldex.effect <- function(clr, verbose = TRUE, include.sample.summary = FALSE, useMC = FALSE,
                         CI = FALSE, glm.conds = NULL, paired.test = FALSE) {
  
  
  # returns the median differences in abundance between 2 conditions
  # returns the median effect size and proportion of effect that overlaps 0
  # data is returned in a data frame
  # requires multicore
  # this uses Rfast
  
  # Use clr conditions slot instead of input
  if (is.vector(clr@conds)) {
    conditions <- clr@conds
  } else if (is.factor(clr@conds)) {
    if (length(levels(clr@conds) == 2)) {
      conditions <- clr@conds
    }
  } else if (is.matrix(clr@conds)){
    if(is.null(glm.conds)) stop("please provide a binary condition vector")
    conditions <- glm.conds
  } else {
    stop("please check that the conditions parameter for aldex.clr is correct.")
  }
  
  is.multicore = FALSE
  
  if ("BiocParallel" %in% rownames(installed.packages()) & useMC==TRUE){
    message("multicore environment is OK -- using the BiocParallel package")
    #require(BiocParallel)
    is.multicore = TRUE
  }
  else {
    if (verbose == TRUE) message("operating in serial mode")
  }
  
  nr <- numFeatures(clr) # number of features
  rn <- getFeatureNames(clr) # feature names
  # ---------------------------------------------------------------------
  
  # sanity check to ensure only two conditons passed to this function
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )
  
  if ( length( conditions ) !=  numConditions(clr) ) stop("mismatch btw 'length(conditions)' and 'ncol(reads)'")
  
  if ( length( levels ) != 2 ) stop("only two condition levels are currently supported")
  
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  
  for ( l in levels( conditions ) ) {
    levels[[l]] <- which( conditions == l )
    if ( length( levels[[l]] ) < 2 ) stop("condition level '",l,"' has less than two replicates")
  }
  
  # end sanity check
  if (verbose == TRUE) message("sanity check complete")
  
  # Summarize the relative abundance (rab) win and all groups
  
  rab <- vector( "list", 3 )
  names(rab) <- c( "all", "win", "spl" )
  rab$win <- list()
  
  #this is the median value across all monte carlo replicates
  # for loops replaced with do.call
  cl2p <- NULL
  cl2p <- do.call(cbind, getMonteCarloInstances(clr))
  # this is a 2X speedup
  rab$all <- Rfast::rowMedians(cl2p)
  names(rab$all) <- rownames(cl2p)
  rm(cl2p)
  
  if (verbose == TRUE) message("rab.all  complete")
  
  for(level in levels(conditions)){
    cl2p <- NULL
    cl2p <- do.call(cbind, getMonteCarloInstances(clr)[levels[[level]]] )
    rab$win[[level]] <-  Rfast::rowMedians(cl2p)
    rm(cl2p)
  }
  
  if (verbose == TRUE) message("rab.win  complete")
  
  if (is.multicore == TRUE)  rab$spl <- bplapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )
  #RMV if (is.multicore == FALSE) rab$spl <- lapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )
  if (is.multicore == FALSE) rab$spl <- lapply( getMonteCarloInstances(clr), function(m) { Rfast::rowMedians(m) } )
  if (verbose == TRUE) message("rab of samples complete")
  
  # ---------------------------------------------------------------------
  
  # make and fill in the data table
  # i know this is inefficient, but it works and is not a bottleneck
  if(CI == FALSE) {
    rv <- list(
      rab = rab
    )
  } else {
    rv <- list(
      rab = rab
    )
  }
  
  if (verbose == TRUE) message("summarizing output")
  
  y.rv <- data.frame(rv$rab$all)
  colnames(y.rv) <- c("rab.all")
  for(i in names(rv$rab$win)){
    nm <- paste("rab.win", i, sep=".")
    y.rv[,nm] <- data.frame(rv$rab$win[[i]])
  }
  if (include.sample.summary == TRUE){
    for(i in names(rv$rab$spl)){
      nm <- paste("rab.sample", i, sep=".")
      if (is.multicore == TRUE) y.rv[,nm] <- data.frame(t(rv$rab$spl[[i]]))
      if (is.multicore == FALSE) y.rv[,nm] <- data.frame(rv$rab$spl[[i]])
    }
    
  }
  
  return(y.rv)
  
}