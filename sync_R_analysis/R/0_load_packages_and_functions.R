
setwd()

# Libraries used for ALL sections ==========

# libbrarys 

library(tidyverse)
library(dplyr)
library(fda)
library(fdaoutlier)
library(dtw)
library(reshape2)
library(data.table)
library(matrixStats)
library(ggplot2)
library(readxl)
library(igraph)
library(biomaRt)
library(patchwork)
library(lomb)
library(stats)
library(pracma)
library(qtl2ggplot)
library(tidyr)
library(dplyr)
library(qtl)

# Functions for use throughout =========

standard_error <- function(x) {
  n <- sum(!is.na(x))  # Number of non-NA values
  sd <- sd(x, na.rm = TRUE)  # Standard deviation excluding NA values
  se <- sd / sqrt(n)  # Standard error
  return(se)
}



#calculate zscore =========

get_z_score_matrix <- function(df_in, include_na = FALSE) {
  # Input should be a data frame containing only numeric variables
  # Use include_na = TRUE to include rows with 0 standard deviation

  # Calculate z scores - i.e. transform each row by subtracting mean
  # then dividing by the standard deviation
  df_z <- (df_in - rowMeans(df_in)) / (rowSds(as.matrix(df_in)))[row(df_in)]

  # remove values with Na - i.e. anything where the standard deviation was 0
  if (!include_na) {
    df_z <- df_z %>% drop_na()
  }

  return(df_z)
}

#filter curves for fda input ==========

filter_curves <- function(x, time, from = min(time, na.rm = TRUE), to = max(time, na.rm = TRUE), min_lum = 0, periodogram_length = NULL) {
  
  if (!any(apply(x, 2, is.numeric)))
    stop("values are NOT numeric. This is not permitted")
  
  if (!is.matrix(time))
    stop("Time is NOT matrix array. This is not permitted")
  
  length_df <- length(x)
  rownumber <- function(time, from, to ) {
    minrow <- which(abs(time - from) == min(abs(time - from)))
    maxrow <- which(abs(time - to) == min(abs(time - to)))
    return(c(minrow,maxrow))
  }
  
  #create data frame for curves and time vector with number of rows based of to and form
  Curves_Timecut <- data.frame(x[rownumber(time = time, from = from, to = to)[1]:rownumber(time = time, from = from, to = to)[2],1:length_df])
  
  time_cut <- time[rownumber(time = time,from = from,to = to)[1]:rownumber(time = time, from = from, to = to)[2],]
  
  # create df of curves that exceed the minimum luminescence threshold
  Curves_set_min <- Curves_Timecut[sapply(Curves_Timecut, function(x) min(x, na.rm = T) > min_lum)]
  
  #trend and then detrend the curves and keep column names the same
  tred <- apply(Curves_set_min, MARGIN = 2, FUN = function(X) stats::lm(X~time_cut, na.action = na.exclude))
  detrend <- as.data.frame(sapply(X = seq(1:length(Curves_set_min)), FUN = function(g){ residuals(tred[[g]])}))
  colnames(detrend) <- colnames(Curves_set_min)
  
  #depending on length specified to perfom periodogram
  
  # Calculate length of the periodogram
  if (is.null(periodogram_length)) {
    periodogram_from <- from
    periodogram_to <- to
  } else {
    if (periodogram_length > to - from) {
      stop("Periodogram length exceeds the length between 'from' and 'to'.")
    }
    if (periodogram_length < 18) {
      stop("Periodogram length must be greater than 18")
    }
    periodogram_from <- max(time, na.rm = TRUE) - periodogram_length
    periodogram_to <- max(time, na.rm = TRUE)
  }
  
  
  #create dataframe of time from the last number of hours specified by periodogram length
  time_periodogram <- time_cut[rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[1]:rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[2]]
  
  
  #create dataframe of the last hours of the curves for periodogram analysis specified by periodogram length
  detrend_periodogram_df <- data.frame(detrend[rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[1]:rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[2],])
  
  #perform ls periodogram then find the peaks of these and remove any curves where this is below 0.001 threshold
  pero <- apply(detrend_periodogram_df, MARGIN = 2, FUN = function(x) lomb::lsp(x, times = time_periodogram, type =  "period",from = 18, to = 30,ofac = 44, plot = FALSE ))
  
  peaks <- lapply(X = seq(1:length(Curves_set_min)), FUN = function(g) {pracma::findpeaks(pero[[g]]$power, nups = 1, ndowns = 1, zero = "0")})
  
  peaks <- lapply(X = seq(1:length(Curves_set_min)), FUN = function(g) {max(peaks[[g]][,1])})
  
  #this removes curves that do not meet the periodogram threshold
  Curves_filterd <- detrend[sapply(X = seq(1:length(detrend)), FUN = function(g) {peaks[[g]] > pero[[g]]$sig.level})]
  finaldf<-cbind(time_cut,Curves_filterd)
  return(finaldf)
  
}

# function for smoothing curves ==================

smooth_fun <- function(x, time, outlier = TRUE, shape = TRUE) {
  
  if (!any(apply(x, 2, is.numeric)))
    stop("values are NOT numeric. This is not permitted")
  
  if (!is.matrix(time))
    stop("Time is NOT matrix array. This is not permitted")
  # NORMALISE each curve
  norm_all <- apply(x, MARGIN = 2, FUN = function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))
  #create list containing each column and the time column with NAs removed
  everycurve <- list()
  for (i in 1:ncol(norm_all)) {
    
    #combine time column with every every column individually
    combine <- cbind(time,norm_all[,i])
    #append list with each time and curve combo with NAs removed.
    everycurve[[i]] <- na.omit(combine)
  }
  #now make a basis functions and functional parameter for each one....
  sp_fdobj <- list()
  for (i in 1:ncol(norm_all)) {
    #range and number of basis for each bassis functtions
    range <- c(min(everycurve[[i]][,1]), max(everycurve[[i]][,1]))
    basis_num <- length(everycurve[[i]][,1]) + 2
    # basis function
    sp_BASIS_bspline <- fda::create.bspline.basis(rangeval = range ,norder = 4,nbasis = basis_num)
    #Define Functional Parameter with Roughness Penalty
    sp_fdobj[[i]] <- fda::fdPar(sp_BASIS_bspline,Lfdobj = 2,lambda = 1)
  }
  # for each curve make the smoothbasis
  sp_smoothbasis <- list()
  for (i in 1:ncol(norm_all)) {
    sp_smoothbasis[[i]] <- fda::smooth.basis(argvals = everycurve[[i]][,1], y = everycurve[[i]][,2], fdParobj = sp_fdobj[[i]])
  }
  #make a new time vector that will cover all the different times
  #the max of the minamum ZT time and the minamum of the maximum time
  min_time <- max(sapply(everycurve, function(x) min(x[,1])))
  max_time <- min(sapply(everycurve, function(x) max(x[,1])))
  new_time <- seq(min_time,max_time,length.out = 300)
  
  #evaluate the curves all under this new time
  eval_curves <- sapply(sp_smoothbasis, function(x) fda::eval.fd(new_time,x$fd))
  colnames(eval_curves) <- colnames(norm_all)
  
  #now create smooth functions for all cuvres together!
  Total_curves <- data.matrix(eval_curves)
  BASIS_bspline <- fda::create.bspline.basis(rangeval = c(min(new_time),max(new_time)),norder = 4,nbasis = 302)
  fdobj_outlier = fda::fdPar(BASIS_bspline,Lfdobj = 2,lambda = 50)
  sp_totalsmooth_outlier <- fda::smooth.basis(argvals = new_time, y = Total_curves, fdParobj = fdobj_outlier)
  
  
  v_sp_totalsmooth_outlier <- fda::eval.fd(new_time,sp_totalsmooth_outlier$fd, outlier)
  #transpose and do tvs mss outlier detection....
  
  if (shape) {
    sp_transposed_curves <- data.matrix(t(v_sp_totalsmooth_outlier))
    
    
    tvoutlier <- fdaoutlier::tvdmss(dts = sp_transposed_curves)
    
    #tvoutlier$shape_outliers
    #tvoutlier$magnitude_outliers
    if (!is.null(tvoutlier$shape_outliers)) {
      Total_curves <- data.matrix(Total_curves[, -tvoutlier$shape_outliers])
      
    }
  }
  fdobj = fda::fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
  sp_totalsmooth <- fda::smooth.basis(argvals = new_time, y = Total_curves, fdParobj = fdobj)
  #list_both <- list(sp_totalsmooth,new_time,Total_curves)
  return(sp_totalsmooth)
  return(list_both)
  return(new_time)
  return(tvoutlier)
}

#

