# time warping of full length curves ==============

# load in data with outliers removed 
Total_curves_SD_LD_rmout <- read.csv(file = "data/Total_curves_SDLD.csv")
Total_curves_SD_SD_rmout <- read.csv(file = "data/Total_curves_SDSD.csv")
Total_curves_LD_LD_rmout <- read.csv(file = "data/Total_curves_LDLD.csv")
#remove first index column 
new_ZT <- Total_curves_SD_LD_rmout[,1]
Total_curves_SD_LD_rmout <- Total_curves_SD_LD_rmout[,-1]
Total_curves_SD_SD_rmout <- Total_curves_SD_SD_rmout[,-1]
Total_curves_LD_LD_rmout <- Total_curves_LD_LD_rmout[,-1]

#load in mean curves 
Genotype_meansSDLD <- read.csv("data/genotype_means_SDLD.csv")
Genotype_meansSDSD <- read.csv("data/genotype_means_SDSD.csv")
Genotype_meansLDLD <- read.csv("data/genotype_means_LDLD.csv")
#remove first index column 
Genotype_meansSDLD <- Genotype_meansSDLD[,-1]
Genotype_meansSDSD <- Genotype_meansSDSD[,-1]
Genotype_meansLDLD <- Genotype_meansLDLD[,-1]

#make functional object with the genotype means and full length curves 

#make the my vetor twice as there is the SD and LD portion
my_vector_SDLD <- vector("numeric")
my_vector_SDSD <- vector("numeric")
my_vector_LDLD <- vector("numeric")


Total_curves_SD_LD_rmout <- data.frame(Total_curves_SD_LD_rmout)
Total_curves_SD_SD_rmout <- data.frame(Total_curves_SD_SD_rmout)
Total_curves_LD_LD_rmout <- data.frame(Total_curves_LD_LD_rmout)

Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names <- sort(Genotype_names)
library(dplyr)
for (i in Genotype_names) {
  my_vector_SDLD[i] <- length(Total_curves_SD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_SDSD[i] <- length(Total_curves_SD_SD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_LDLD[i] <- length(Total_curves_LD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
}
#make into functional data objects 

genotype_listSDLD <- list()
genotype_listSDSD <- list()
genotype_listLDLD <- list()
#genotype_listcontrols <- list()

for (i in Genotype_names) {
  
  genotype_listSDLD[[i]] <- Total_curves_SD_LD_rmout %>%
    dplyr::select(dplyr::contains(i))
  
  genotype_listSDSD[[i]] <- Total_curves_SD_SD_rmout %>%
    dplyr::select(dplyr::contains(i))
  
  genotype_listLDLD[[i]] <- Total_curves_LD_LD_rmout %>%
    dplyr::select(dplyr::contains(i))
  
}

#need each aspect of yhe list to be numeric so...

for (i in 1:length(genotype_listSDLD)) {
  genotype_listSDLD[[i]] <- data.matrix(genotype_listSDLD[[i]])
  genotype_listSDSD[[i]] <- data.matrix(genotype_listSDSD[[i]])
  genotype_listLDLD[[i]] <- data.matrix(genotype_listLDLD[[i]])
}

#Now have a list again where each component is a genotype, so now for each genotype in this list must creast a smooth basis and then a functional data object.

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 200)
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)

genotype_smoothSDLD <- lapply(X = seq(1,64), FUN = function(x){smooth.basis(argvals = new_ZT, y = genotype_listSDLD[[x]], fdParobj = fdobj)})

genotype_smoothSDSD <- lapply(X = seq(1,64), FUN = function(x){smooth.basis(argvals = new_ZT, y = genotype_listSDSD[[x]], fdParobj = fdobj)})

genotype_smoothLDLD <- lapply(X = seq(1,64), FUN = function(x){smooth.basis(argvals = new_ZT, y = genotype_listLDLD[[x]], fdParobj = fdobj)})

#do the no derivative for all curves 
genotype_smooth_oriSDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothSDLD[[x]]$fd,0)})
genotype_smooth_oriSDSD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothSDSD[[x]]$fd,0)})
genotype_smooth_oriLDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothLDLD[[x]]$fd,0)})

#for means 
genotype_meansSDLD <- list()
genotype_meansSDSD <- list()
genotype_meansLDLD <- list()

for (i in (1:64)) {
  
  genotype_meansSDLD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansSDLD[,i], fdParobj = fdobj)
  genotype_meansSDSD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansSDSD[,i], fdParobj = fdobj)
  genotype_meansLDLD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansLDLD[,i], fdParobj = fdobj)
}

#no derivative for gneotpye means 

genotype_means_oriSDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansSDLD[[x]]$fd,0)})
genotype_means_oriSDSD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansSDSD[[x]]$fd,0)})
genotype_means_oriLDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansLDLD[[x]]$fd,0)})


#For registering the curves in each gneotype to its mean need to define smoothing parameter for warping functions 
Wbasis    <- create.bspline.basis(rangeval = c(56,152),norder = 4,nbasis = 20)
Wfd0      <- fd(matrix(0,20,1),Wbasis)
#  set up the functional parameter object using only
#      a light amount smoothing
WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=0.1)

#register the functions

regSDLD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_oriSDLD[[x]], yfd= genotype_smooth_oriSDLD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

regSDSD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_oriSDSD[[x]], yfd= genotype_smooth_oriSDSD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

regLDLD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_oriLDLD[[x]], yfd= genotype_smooth_oriLDLD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

#save the output 
save(regSDLD, file = "reg_fulllength_to_means_oriSDLD.RData")
save(regSDSD, file = "reg_fulllength_to_means_oriSDSD.RData")
save(regLDLD, file = "reg_fulllength_to_means_oriLDLD.RData")


#    For velocity curves change derivative to = 1 =====================

genotype_smooth_velSDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothSDLD[[x]]$fd,1)})
genotype_smooth_velSDSD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothSDSD[[x]]$fd,1)})
genotype_smooth_velLDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothLDLD[[x]]$fd,1)})


#for means 
genotype_meansSDLD <- list()
genotype_meansSDSD <- list()
genotype_meansLDLD <- list()

for (i in (1:64)) {
  
  genotype_meansSDLD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansSDLD[,i], fdParobj = fdobj)
  genotype_meansSDSD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansSDSD[,i], fdParobj = fdobj)
  genotype_meansLDLD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansLDLD[,i], fdParobj = fdobj)
}

#no derivative for gneotpye means 
genotype_means_velSDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansSDLD[[x]]$fd,1)})
genotype_means_velSDSD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansSDSD[[x]]$fd,1)})
genotype_means_velLDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansLDLD[[x]]$fd,1)})

#For registering the curves in each gneotype to its mean need to define smoothing parameter for warping functions 
Wbasis    <- create.bspline.basis(rangeval = c(56,152),norder = 4,nbasis = 20)
Wfd0      <- fd(matrix(0,20,1),Wbasis)
#  set up the functional parameter object using only
#      a light amount smoothing
WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=0.1)

#register the functions

regSDLD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_velSDLD[[x]], yfd= genotype_smooth_velSDLD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

regSDSD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_velSDSD[[x]], yfd= genotype_smooth_velSDSD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

regLDLD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_velLDLD[[x]], yfd= genotype_smooth_velLDLD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

#save the output 

save(regSDLD, file = "reg_fulllength_to_means_velSDLD.RData")
save(regSDSD, file = "reg_fulllength_to_means_velSDSD.RData")
save(regLDLD, file = "reg_fulllength_to_means_velLDLD.RData")


#    For acceleration curves change derivative to = 2 =====================


genotype_smooth_accSDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothSDLD[[x]]$fd,2)})
genotype_smooth_accSDSD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothSDSD[[x]]$fd,2)})
genotype_smooth_accLDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_smoothLDLD[[x]]$fd,2)})


#for means 

genotype_meansSDLD <- list()
genotype_meansSDSD <- list()
genotype_meansLDLD <- list()

for (i in (1:64)) {
  
  genotype_meansSDLD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansSDLD[,i], fdParobj = fdobj)
  genotype_meansSDSD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansSDSD[,i], fdParobj = fdobj)
  genotype_meansLDLD[[i]] <- fda::smooth.basis(argvals = new_ZT, y = Genotype_meansLDLD[,i], fdParobj = fdobj)
}

#no derivative for gneotpye means 

genotype_means_accSDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansSDLD[[x]]$fd,2)})
genotype_means_accSDSD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansSDSD[[x]]$fd,2)})
genotype_means_accLDLD <- lapply(X = seq(1,64), FUN = function(x){deriv.fd(genotype_meansLDLD[[x]]$fd,2)})


#For registering the curves in each gneotype to its mean need to define smoothing parameter for warping functions 
Wbasis    <- create.bspline.basis(rangeval = c(56,152),norder = 4,nbasis = 20)
Wfd0      <- fd(matrix(0,20,1),Wbasis)
#  set up the functional parameter object using only
#      a light amount smoothing
WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=0.1)

#register the functions

regSDLD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_accSDLD[[x]], yfd= genotype_smooth_accSDLD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

regSDSD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_accSDSD[[x]], yfd= genotype_smooth_accSDSD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

regLDLD <- lapply(X = seq(1,64), FUN = function(x){register.fd(y0fd=genotype_means_accLDLD[[x]], yfd= genotype_smooth_accLDLD[[x]], WfdParobj= WfdParobj, conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)})

#save the output 

save(regSDLD, file = "reg_fulllength_to_means_accSDLD.RData")
save(regSDSD, file = "reg_fulllength_to_means_accSDSD.RData")
save(regLDLD, file = "reg_fulllength_to_means_accLDLD.RData")


#

