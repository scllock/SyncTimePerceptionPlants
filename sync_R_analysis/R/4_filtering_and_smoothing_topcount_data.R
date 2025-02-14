# filtering and smoothing of topcount data =============== 

#set working directory to folder contaning all the SD_LD plates. 
setwd("data/topcount/SD_LD")

#load in raw data 
temp = list.files(pattern = "plate_*")
myfiles = lapply(temp, read.csv)

#make list containing the times for each plate (will be 16) and list for the plates 
timesZT <- lapply(X = seq(1:length(myfiles)), FUN = function(g){ data.matrix(myfiles[[g]][,1])})
plates <- lapply(X = seq(1:length(myfiles)), FUN = function(g){ myfiles[[g]][,-1]})

#Run filter curves function (use 55-153) with periodogram on the last 48 hours
num <- lapply(X = seq(1:length(plates)), FUN = function(i){filter_curves(plates[[i]], time = timesZT[[i]], from = 55, to = 153, min_lum = 200, periodogram_length = 48)})

#remove time as fist column
timesZT <- lapply(X = seq(1:length(num)), FUN = function(g){ data.matrix(num[[g]][,1])})
plates <- lapply(X = seq(1:length(num)), FUN = function(g){ data.matrix(num[[g]][,-1])})

#Run filter curves function (use 55-153) with periodogram on the last 48 hours
smoothing <- lapply(X = seq(1:length(plates)), FUN = function(i){smooth_fun(plates[[i]], time = timesZT[[i]], outlier = 1, shape = TRUE)})

# then evaluate each one over a new time frame (56,152)

new_timezt <- seq(56,152, length.out = 300)
eval_smooth <- lapply(X=seq(1:length(smoothing)), FUN = function(i){eval.fd(new_timezt,smoothing[[i]]$fd,0)})

#combine all curves together and removelist
Total_curves <- data.frame(bind_cols(eval_smooth))

#sort based on column names
Total_curves_SDLD <- Total_curves[,sort(colnames(Total_curves))]


# same for SDSD condition ========

#set working directory to folder contaning all the SD_SD plates. 
setwd("data/topcount/SD_SD")

#load in all plates
temp = list.files(pattern = "plate_*")
myfiles = lapply(temp, read.csv)

#make list containing the times for each plate (will be 16) and list for the plates 
timesZT <- lapply(X = seq(1:length(myfiles)), FUN = function(g){ data.matrix(myfiles[[g]][,1])})
plates <- lapply(X = seq(1:length(myfiles)), FUN = function(g){ myfiles[[g]][,-1]})

#Run through my filter curves function (use 55-153)
num <- lapply(X = seq(1:length(plates)), FUN = function(i){filter_curves(plates[[i]], time = timesZT[[i]], from = 55, to = 153, min_lum = 200,periodogram_length = 48)})

#remove time as fist colum
timesZT <- lapply(X = seq(1:length(num)), FUN = function(g){ data.matrix(num[[g]][,1])})
plates <- lapply(X = seq(1:length(num)), FUN = function(g){ data.matrix(num[[g]][,-1])})

#Run filter curves function (use 55-153) with periodogram on the last 48 hours
smoothing <- lapply(X = seq(1:length(plates)), FUN = function(i){smooth_fun(plates[[i]], time = timesZT[[i]], outlier = 1, shape = TRUE)})

# then evaluate each one over a new time frame (56,152)
new_timezt <- seq(56,152, length.out = 300)
eval_smooth <- lapply(X=seq(1:length(smoothing)), FUN = function(i){eval.fd(new_timezt,smoothing[[i]]$fd,0)})

#combine all curves together and removelist
Total_curves <- data.frame(bind_cols(eval_smooth))

#sort based on column names
Total_curves_SDSD <- Total_curves[,sort(colnames(Total_curves))]


#same for LD_LD plates 
setwd("data/topcount/LD_LD")

#load in all plates
temp = list.files(pattern = "plate_*")
myfiles = lapply(temp, read.csv)

#make list containing the times for each plate (will be 16) and list for the plates 
timesZT <- lapply(X = seq(1:length(myfiles)), FUN = function(g){ data.matrix(myfiles[[g]][,1])})
plates <- lapply(X = seq(1:length(myfiles)), FUN = function(g){ myfiles[[g]][,-1]})

#Run through my filter curves function (use 55-153)
num <- lapply(X = seq(1:length(plates)), FUN = function(i){filter_curves(plates[[i]], time = timesZT[[i]], from = 55, to = 153, min_lum = 200, periodogram_length =  48)})

#remove time as fist colum
timesZT <- lapply(X = seq(1:length(num)), FUN = function(g){ data.matrix(num[[g]][,1])})
plates <- lapply(X = seq(1:length(num)), FUN = function(g){ data.matrix(num[[g]][,-1])})

#Run filter curves function (use 55-153) with periodogram on the last 48 hours
smoothing <- lapply(X = seq(1:length(plates)), FUN = function(i){smooth_fun(plates[[i]], time = timesZT[[i]], outlier = 1, shape = TRUE)})

# then evaluate each one over a new time frame (56,152)
new_timezt <- seq(56,152, length.out = 300)
eval_smooth <- lapply(X=seq(1:length(smoothing)), FUN = function(i){eval.fd(new_timezt,smoothing[[i]]$fd,0)})

#combine all curves together and removelist
Total_curves_LDLD <- data.frame(bind_cols(eval_smooth))

#sort based on column names
Total_curves_LDLD <- Total_curves_LDLD[,sort(colnames(Total_curves_LDLD))]

#for final filtering run fdaoutlier on all conditions combined=============

#combine all the datasets to filter for outliers.
Total_curves_ALL <- cbind(Total_curves_SDLD,Total_curves_SDSD, Total_curves_LDLD)
Total_curves_ALL <- data.matrix(Total_curves_ALL)

#create basis function
BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_timezt),max(new_timezt)),norder = 4,nbasis = 300)

#using high smothing to get overall shape of curves (lambda = 20)
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 20)
sp_totalsmooth <- smooth.basis(argvals = new_timezt, y = Total_curves_ALL, fdParobj = fdobj)

o_sp_totalsmooth <- eval.fd(new_time,sp_totalsmooth$fd,0)
v_sp_totalsmooth <- eval.fd(new_time,sp_totalsmooth$fd,Lfdobj = 1)
a_sp_totalsmooth <- eval.fd(new_time,sp_totalsmooth$fd,Lfdobj = 2)

#transpose and do tvs mss outlier detection using second derivative
sp_transposed_curves <- t(a_sp_totalsmooth)

sp_transposed_curves <- data.matrix(sp_transposed_curves)

tvoutlier <- tvdmss(dts = sp_transposed_curves)

tvoutlier$shape_outliers
tvoutlier$magnitude_outliers


#remove the correct curves and split back into the different groups 
#SD_LD
SD_LD_remove <- tvoutlier$shape_outliers[tvoutlier$shape_outliers < length(Total_curves_SDLD)]

#SD_SD
SD_SD_remove <- tvoutlier$shape_outliers[tvoutlier$shape_outliers > length(Total_curves_SDLD) & tvoutlier$shape_outliers < (length(Total_curves_SDLD) + length(Total_curves_SDSD)) ] - length(Total_curves_SDLD)

#LD_LD
LD_LD_remove <- tvoutlier$shape_outliers[tvoutlier$shape_outliers > (length(Total_curves_SDLD) + length(Total_curves_SDSD))] - (length(Total_curves_SDLD) + length(Total_curves_SDSD))

# final set with the removal of the outliers from each section (rmout mean removal of outliers)
Total_curves_SD_LD_rmout <- Total_curves_SDLD[,- SD_LD_remove]
Total_curves_SD_SD_rmout <- Total_curves_SDSD[,- SD_SD_remove]
Total_curves_LD_LD_rmout <- Total_curves_LDLD[,- LD_LD_remove]

#save each dataframe combined with the new time
forsavingSDLD <- cbind(new_timezt, Total_curves_SD_LD_rmout)

forsavingSDSD <- cbind(new_timezt, Total_curves_SD_SD_rmout)

forsavingLDLD <- cbind(new_timezt, Total_curves_LD_LD_rmout)

#save
write.csv(forsavingSDLD,"Total_curves_SDLD.csv", row.names = F)

write.csv(forsavingSDSD,"Total_curves_SDSD.csv", row.names = F)

write.csv(forsavingLDLD,"Total_curves_LDLD.csv", row.names = F)



# calculate the dataframes for curves in each condition pre and post shift  

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_timezt),max(new_timezt)),norder = 4,nbasis = 300)
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)

#make datamaxtric  
Total_curves_SD_LD_rmout <- data.matrix(Total_curves_SD_LD_rmout)
Total_curves_SD_SD_rmout <- data.matrix(Total_curves_SD_SD_rmout)
Total_curves_LD_LD_rmout <- data.matrix(Total_curves_LD_LD_rmout)


#smooth basis 
sp_totalsmooth_SDLD <- smooth.basis(argvals = new_timezt, y = Total_curves_SD_LD_rmout, fdParobj = fdobj)
sp_totalsmooth_SDSD <- smooth.basis(argvals = new_timezt, y = Total_curves_SD_SD_rmout, fdParobj = fdobj)
sp_totalsmooth_LDLD <- smooth.basis(argvals = new_timezt, y = Total_curves_LD_LD_rmout, fdParobj = fdobj)

#evalulate curves for post and pre shift then combine into one dataframe.

#FIRST ------- evaluate the curves from 48 hours before the missed cue..
time_pre <- seq(56,104, length.out = 100)
SDLD_pre_ori <- eval.fd(time_pre, sp_totalsmooth_SDLD$fd)
SDSD_pre_ori <- eval.fd(time_pre, sp_totalsmooth_SDSD$fd)
LDLD_pre_ori <- eval.fd(time_pre, sp_totalsmooth_LDLD$fd)

#SECOND ------- evaluate the curves from 48 hours AFTER the missed cue..
time_post <- seq(104,152, length.out = 100)
SDLD_post_ori <- eval.fd(evalarg = time_post, sp_totalsmooth_SDLD$fd)
SDSD_post_ori <- eval.fd(evalarg = time_post, sp_totalsmooth_SDSD$fd)
LDLD_post_ori <- eval.fd(evalarg = time_post, sp_totalsmooth_LDLD$fd)

#combine all should be two per genotype....

curves_pre_post_SDLD <- cbind(SDLD_pre_ori,SDLD_post_ori)
curves_pre_post_SDSD <- cbind(SDSD_pre_ori,SDSD_post_ori)
curves_pre_post_LDLD <- cbind(LDLD_pre_ori,LDLD_post_ori)

write.csv(curves_pre_post_SDLD, file = "curves_pre_post_SDLD.csv")
write.csv(curves_pre_post_SDSD, file = "curves_pre_post_SDSD.csv")
write.csv(curves_pre_post_LDLD, file = "curves_pre_post_LDLD.csv")
#

# making genotype means ===========


Total_curves_SD_LD_rmout <- data.frame(Total_curves_SD_LD)
Total_curves_SD_SD_rmout <- data.frame(Total_curves_SD_SD)
Total_curves_LD_LD_rmout <- data.frame(Total_curves_LD_LD)
Total_curves_controls <- cbind(Total_curves_SD_SD, Total_curves_LD_LD)

Genotype_names <- sort(Genotype_names)
genotype_listSDLD <- list()
genotype_listSDSD <- list()
genotype_listLDLD <- list()
genotype_listcontrols <- list()

for (i in Genotype_names) {
  
  genotype_listSDLD[[i]] <- Total_curves_SD_LD_rmout %>%
    dplyr::select(dplyr::contains(i))
  
  genotype_listSDSD[[i]] <- Total_curves_SD_SD_rmout %>%
    dplyr::select(dplyr::contains(i))
  
  genotype_listLDLD[[i]] <- Total_curves_LD_LD_rmout %>%
    dplyr::select(dplyr::contains(i))
  
  genotype_listcontrols[[i]] <- Total_curves_controls %>%
    dplyr::select(dplyr::contains(i))
}

#need each aspect of yhe list to be numeric so...

for (i in 1:length(genotype_listSDLD)) {
  genotype_listSDLD[[i]] <- data.matrix(genotype_listSDLD[[i]])
  genotype_listSDSD[[i]] <- data.matrix(genotype_listSDSD[[i]])
  genotype_listLDLD[[i]] <- data.matrix(genotype_listLDLD[[i]])
  genotype_listcontrols[[i]] <- data.matrix(genotype_listcontrols[[i]])
}

#Now have a list again where each component is a genotype, so now for each genotype in this list must creast a smooth basis and then a functional data object.

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 200)
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)

genotype_smoothSDLD <- lapply(X = seq(1,64), FUN = function(x){smooth.basis(argvals = new_ZT, y = genotype_listSDLD[[x]], fdParobj = fdobj)})

genotype_smoothSDSD <- lapply(X = seq(1,64), FUN = function(x){smooth.basis(argvals = new_ZT, y = genotype_listSDSD[[x]], fdParobj = fdobj)})

genotype_smoothLDLD <- lapply(X = seq(1,64), FUN = function(x){smooth.basis(argvals = new_ZT, y = genotype_listLDLD[[x]], fdParobj = fdobj)})

genotype_smoothcontrols <- lapply(X = seq(1,length(genotype_listcontrols)), FUN = function(x){smooth.basis(argvals = new_ZT, y = genotype_listcontrols[[x]], fdParobj = fdobj)})

#now for each one take the first derivative and find a standard deviation curve

###SDLD
mean_genotype_oriSDLD <- lapply(X = seq(1,64), FUN = function(x){mean.fd(genotype_smoothSDLD[[x]]$fd)})

###SDSD
mean_genotype_oriSDSD <- lapply(X = seq(1,64), FUN = function(x){mean.fd(genotype_smoothSDSD[[x]]$fd)})

###LDLD
mean_genotype_oriLDLD <- lapply(X = seq(1,64), FUN = function(x){mean.fd(genotype_smoothLDLD[[x]]$fd)})

#evalulate each of the mean functions under the same time frame
means_dfSDLD <- sapply(mean_genotype_oriSDLD, FUN = function(i) {fda::eval.fd(new_ZT,i)})
colnames(means_dfSDLD) <- Genotype_names

means_dfSDSD <- sapply(mean_genotype_oriSDSD, FUN = function(i) {fda::eval.fd(new_ZT,i)})
colnames(means_dfSDSD) <- Genotype_names

means_dfLDLD <- sapply(mean_genotype_oriLDLD, FUN = function(i) {fda::eval.fd(new_ZT,i)})
colnames(means_dfLDLD) <- Genotype_names

#save means as a csv file 

write.csv(means_dfSDLD, "genotype_means_SDLD.csv")
write.csv(means_dfSDSD, "genotype_means_SDSD.csv")
write.csv(means_dfLDLD, "genotype_means_LDLD.csv")





