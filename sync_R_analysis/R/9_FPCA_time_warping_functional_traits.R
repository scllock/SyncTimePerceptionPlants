# script for reading in the registered curves and performing FPCA and amp phase decomp to get trait table for warping on full length and pre and post curves ===============

#setwd 

#read in the registered data 

Total_curves_SD_LD_rmout <- read.csv("56_152/full_length_curves/origional/data/Total_curves_SDLD.csv")
Total_curves_SD_SD_rmout <- read.csv("56_152/full_length_curves/origional/data/Total_curves_SDSD.csv")
Total_curves_LD_LD_rmout <- read.csv("56_152/full_length_curves/origional/data/Total_curves_LDLD.csv")

Total_curves_SD_LD_rmout <- Total_curves_SD_LD_rmout[,-1]
Total_curves_SD_SD_rmout <- Total_curves_SD_SD_rmout[,-1]
Total_curves_LD_LD_rmout <- Total_curves_LD_LD_rmout[,-1]

#origional curves
load(file = "56_152/full_length_curves/origional/reg_fulllength_to_means_oriSDLD.RData")
load(file = "56_152/full_length_curves/origional/reg_fulllength_to_means_oriSDSD.RData")
load(file = "56_152/full_length_curves/origional/reg_fulllength_to_means_oriLDLD.RData")

reg_oriSDLD <- regSDLD
reg_oriSDSD <- regSDSD
reg_oriLDLD <- regLDLD

#velocity curves
load(file = "56_152/full_length_curves/velocity/reg_fulllength_to_means_velSDLD.RData")
load(file = "56_152/full_length_curves/velocity//reg_fulllength_to_means_velSDSD.RData")
load(file = "56_152/full_length_curves/velocity//reg_fulllength_to_means_velLDLD.RData")

reg_velSDLD <- regSDLD
reg_velSDSD <- regSDSD
reg_velLDLD <- regLDLD

#acceleration curves
load(file = "56_152/full_length_curves/acceleration/reg_fulllength_to_means_accSDLD.RData")
load(file = "56_152/full_length_curves/acceleration//reg_fulllength_to_means_accSDSD.RData")
load(file = "56_152/full_length_curves/acceleration//reg_fulllength_to_means_accLDLD.RData")

reg_accSDLD <- regSDLD
reg_accSDSD <- regSDSD
reg_accLDLD <- regLDLD


#plot(reg_oriSDLD[[1]]$regfd, col = "black")
#plot(reg_oriSDLD[[1]]$yfd, col = "black")
#plot(reg_vel[[1]]$regfd)
#plot(reg_acc[[1]]$regfd)


Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names <- sort(Genotype_names)



################################ single traits 


#first do the Amp phase decomposition within genotypes SDLD
decomp_oriSDLD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_oriSDLD[[x]]$yfd, yfd = reg_oriSDLD[[x]]$regfd, hfd =  reg_oriSDLD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_oriSDLD) <- Genotype_names
decomp_oriSDLD <- data.frame(t(decomp_oriSDLD))
decomp_oriSDLD <- as.data.frame(lapply(decomp_oriSDLD, function(x) as.numeric(unlist(x))))


decomp_velSDLD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_velSDLD[[x]]$yfd, yfd = reg_velSDLD[[x]]$regfd, hfd =  reg_velSDLD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_velSDLD) <- Genotype_names
decomp_velSDLD <- data.frame(t(decomp_velSDLD))
decomp_velSDLD <- as.data.frame(lapply(decomp_velSDLD, function(x) as.numeric(unlist(x))))

decomp_accSDLD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_accSDLD[[x]]$yfd, yfd = reg_accSDLD[[x]]$regfd, hfd =  reg_accSDLD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_accSDLD) <- Genotype_names
decomp_accSDLD <- data.frame(t(decomp_accSDLD))
decomp_accSDLD <- as.data.frame(lapply(decomp_accSDLD, function(x) as.numeric(unlist(x))))


#########################

#first do the Amp phase decomposition within genotypes SDSD
decomp_oriSDSD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_oriSDSD[[x]]$yfd, yfd = reg_oriSDSD[[x]]$regfd, hfd =  reg_oriSDSD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_oriSDSD) <- Genotype_names
decomp_oriSDSD <- data.frame(t(decomp_oriSDSD))
decomp_oriSDSD <- as.data.frame(lapply(decomp_oriSDSD, function(x) as.numeric(unlist(x))))

decomp_velSDSD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_velSDSD[[x]]$yfd, yfd = reg_velSDSD[[x]]$regfd, hfd =  reg_velSDSD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_velSDSD) <- Genotype_names
decomp_velSDSD <- data.frame(t(decomp_velSDSD))
decomp_velSDSD <- as.data.frame(lapply(decomp_velSDSD, function(x) as.numeric(unlist(x))))

decomp_accSDSD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_accSDSD[[x]]$yfd, yfd = reg_accSDSD[[x]]$regfd, hfd =  reg_accSDSD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_accSDSD) <- Genotype_names
decomp_accSDSD <- data.frame(t(decomp_accSDSD))
decomp_accSDSD <- as.data.frame(lapply(decomp_accSDSD, function(x) as.numeric(unlist(x))))


###################################

#first do the Amp phase decomposition within genotypes LDLD
decomp_oriLDLD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_oriLDLD[[x]]$yfd, yfd = reg_oriLDLD[[x]]$regfd, hfd =  reg_oriLDLD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_oriLDLD) <- Genotype_names
decomp_oriLDLD <- data.frame(t(decomp_oriLDLD))
decomp_oriLDLD <- as.data.frame(lapply(decomp_oriLDLD, function(x) as.numeric(unlist(x))))

decomp_velLDLD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_velLDLD[[x]]$yfd, yfd = reg_velLDLD[[x]]$regfd, hfd =  reg_velLDLD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_velLDLD) <- Genotype_names
decomp_velLDLD <- data.frame(t(decomp_velLDLD))
decomp_velLDLD <- as.data.frame(lapply(decomp_velLDLD, function(x) as.numeric(unlist(x))))

decomp_accLDLD  <- sapply(X = seq(1,64), FUN = function(x){AmpPhaseDecomp(xfd = reg_accLDLD[[x]]$yfd, yfd = reg_accLDLD[[x]]$regfd, hfd =  reg_accLDLD[[x]]$warpfd, rng=c(56,152))})
colnames(decomp_accLDLD) <- Genotype_names
decomp_accLDLD <- data.frame(t(decomp_accLDLD))
decomp_accLDLD <- as.data.frame(lapply(decomp_accLDLD, function(x) as.numeric(unlist(x))))

###########################################################

#look at the difference between the SDLD ones and the SDSD and LDLD 

diff_oriSDLD_SDSD <- decomp_oriSDLD - decomp_oriSDSD
diff_oriSDLD_LDLD <- decomp_oriSDLD - decomp_oriLDLD

diff_velSDLD_SDSD <- decomp_velSDLD - decomp_velSDSD
diff_velSDLD_LDLD <- decomp_velSDLD - decomp_velLDLD

diff_accSDLD_SDSD <- decomp_accSDLD - decomp_accSDSD
diff_accSDLD_LDLD <- decomp_accSDLD - decomp_accLDLD

#

#          LOOK AT FULL LENGTH TO GENOTYPE MEANS 

#want query to the mean of each genotype must first eval the genotype means over 300 time points 

#need to get the mean of each evalulted and also each curve indiviually 


#for means
query_time <- seq(56,152,length.out = 300)


geontype_meansSDLD <- read.csv("56_152/full_length_curves/origional/data/genotype_means_SDLD.csv")
geontype_meansSDSD <- read.csv("56_152/full_length_curves/origional/data/genotype_means_SDSD.csv")
geontype_meansLDLD <- read.csv("56_152/full_length_curves/origional/data/genotype_means_LDLD.csv")

geontype_meansSDLD <- geontype_meansSDLD[,-1]
geontype_meansSDSD <- geontype_meansSDSD[,-1]
geontype_meansLDLD <- geontype_meansLDLD[,-1]

#for each individual curve 
#split the origional dataframe into genotypes 
genotype_listSDLD <- list()
genotype_listSDSD <- list()
genotype_listLDLD <- list()
Genotype_names_sorted <- sort(Genotype_names)
for (i in Genotype_names_sorted) {
  genotype_listSDLD[[i]] <- Total_curves_SD_LD_rmout %>%
    dplyr::select(dplyr::contains(i))
  genotype_listSDSD[[i]] <- Total_curves_SD_SD_rmout %>%
    dplyr::select(dplyr::contains(i))
  genotype_listLDLD[[i]] <- Total_curves_LD_LD_rmout %>%
    dplyr::select(dplyr::contains(i))
}

#need each aspect of yhe list to be numeric so...

for (i in 1:64) {
  genotype_listSDLD[[i]] <- data.matrix(genotype_listSDLD[[i]])
  genotype_listSDSD[[i]] <- data.matrix(genotype_listSDSD[[i]])
  genotype_listLDLD[[i]] <- data.matrix(genotype_listLDLD[[i]])
}

#query <- genotype_list[[j]][, i]
referenceSDLD <- geontype_meansSDLD
referenceSDSD <- geontype_meansSDSD
referenceLDLD <- geontype_meansLDLD



#make a waping input each has to have two columns the query and then the ref (genotype mean)
warpinputSDLD <- list()
warpinputSDSD <- list()
warpinputLDLD <- list()

for (j in 1:length(genotype_listSDLD)) {
  for (i in 1:ncol(genotype_listSDLD[[j]])) {
    warpinputSDLD[[paste(j, i, sep = "_")]] <- cbind(genotype_listSDLD[[j]][, i],referenceSDLD[, j])
  }
}

for (j in 1:length(genotype_listSDSD)) {
  for (i in 1:ncol(genotype_listSDSD[[j]])) {
    warpinputSDSD[[paste(j, i, sep = "_")]] <- cbind(genotype_listSDSD[[j]][, i],referenceSDSD[, j])
  }
}

for (j in 1:length(genotype_listLDLD)) {
  for (i in 1:ncol(genotype_listLDLD[[j]])) {
    warpinputLDLD[[paste(j, i, sep = "_")]] <- cbind(genotype_listLDLD[[j]][, i],referenceLDLD[, j])
  }
}

#now do dynamic time warping for each curve to its genotype mean 
## Find the best match with the canonical recursion formula
alignmentSDLD <- lapply(X = seq(1:length(warpinputSDLD)), FUN = function(i){dtw(warpinputSDLD[[i]][,1],warpinputSDLD[[i]][,2],keep = TRUE, step.pattern = symmetric2)})

#extract the indext for each curves alignment 
indext1SDLD <- lapply(X = seq(1:length(alignmentSDLD)), FUN = function(i){alignmentSDLD[[i]]$index1})
indext2SDLD <- lapply(X = seq(1:length(alignmentSDLD)), FUN = function(i){alignmentSDLD[[i]]$index2})

#get the normalised distace 
norm_distSDLD <- sapply(X = seq(1:length(alignmentSDLD)), FUN = function(i){alignmentSDLD[[i]]$normalizedDistance})
norm_distSDLD <- data.frame(t(norm_distSDLD))
colnames(norm_distSDLD) <- colnames(Total_curves_SD_LD_rmout)
hist(log10(as.numeric(norm_distSDLD)))
hist(as.numeric(norm_distSDLD))

#SDSD
alignmentSDSD <- lapply(X = seq(1:length(warpinputSDSD)), FUN = function(i){dtw(warpinputSDSD[[i]][,1],warpinputSDSD[[i]][,2],keep = TRUE, step.pattern = symmetric2)})

#extract the indext for each curves alignment 
indext1SDSD <- lapply(X = seq(1:length(alignmentSDSD)), FUN = function(i){alignmentSDSD[[i]]$index1})
indext2SDSD <- lapply(X = seq(1:length(alignmentSDSD)), FUN = function(i){alignmentSDSD[[i]]$index2})

#get the normalised distace 
norm_distSDSD <- sapply(X = seq(1:length(alignmentSDSD)), FUN = function(i){alignmentSDSD[[i]]$normalizedDistance})
norm_distSDSD <- data.frame(t(norm_distSDSD))
colnames(norm_distSDSD) <- colnames(Total_curves_SD_SD_rmout)
hist(log10(as.numeric(norm_distSDSD)))
hist(as.numeric(norm_distSDSD))

#LDLD

alignmentLDLD <- lapply(X = seq(1:length(warpinputLDLD)), FUN = function(i){dtw(warpinputLDLD[[i]][,1],warpinputLDLD[[i]][,2],keep = TRUE, step.pattern = symmetric2)})

#extract the indext for each curves alignment 
indext1LDLD <- lapply(X = seq(1:length(alignmentLDLD)), FUN = function(i){alignmentLDLD[[i]]$index1})
indext2LDLD <- lapply(X = seq(1:length(alignmentLDLD)), FUN = function(i){alignmentLDLD[[i]]$index2})

#get the normalised distace 
norm_distLDLD <- sapply(X = seq(1:length(alignmentLDLD)), FUN = function(i){alignmentLDLD[[i]]$normalizedDistance})
norm_distLDLD <- data.frame(t(norm_distLDLD))
colnames(norm_distLDLD) <- colnames(Total_curves_LD_LD_rmout)
hist(log10(as.numeric(norm_distLDLD)))
hist(as.numeric(norm_distLDLD))

#split by genotype 

genotype_listSDLD <- list()
genotype_listSDLD <- list()
genotype_listSDLD <- list()
for (i in Genotype_names_sorted) {
  genotype_listSDLD[[i]] <- norm_distSDLD %>%
    dplyr::select(dplyr::contains(i))
  genotype_listSDSD[[i]] <- norm_distSDSD %>%
    dplyr::select(dplyr::contains(i))
  genotype_listLDLD[[i]] <- norm_distLDLD %>%
    dplyr::select(dplyr::contains(i))
}

# calculate the means and then the sd of each genotpye 

mean_norm_dsitSDLD <- sapply(X = seq(1:64), FUN = function(i){mean(as.numeric(genotype_listSDLD[[i]][1,]))})
mean_norm_dsitSDSD <- sapply(X = seq(1:64), FUN = function(i){mean(as.numeric(genotype_listSDSD[[i]][1,]))})
mean_norm_dsitLDLD <- sapply(X = seq(1:64), FUN = function(i){mean(as.numeric(genotype_listLDLD[[i]][1,]))})

logmean_norm_dsitSDLD <- sapply(X = seq(1:64), FUN = function(i){mean(log10(as.numeric(genotype_listSDLD[[i]][1,])))})
logmean_norm_dsitSDSD <- sapply(X = seq(1:64), FUN = function(i){mean(log10(as.numeric(genotype_listSDSD[[i]][1,])))})
logmean_norm_dsitLDLD <- sapply(X = seq(1:64), FUN = function(i){mean(log10(as.numeric(genotype_listLDLD[[i]][1,])))})


sd_norm_dsitSDLD <- sapply(X = seq(1:64), FUN = function(i){sd(genotype_listSDLD[[i]])})
sd_norm_dsitSDSD <- sapply(X = seq(1:64), FUN = function(i){sd(genotype_listSDSD[[i]])})
sd_norm_dsitLDLD <- sapply(X = seq(1:64), FUN = function(i){sd(genotype_listLDLD[[i]])})

#

#hist(as.numeric(decomp_ori[,3]))



#FUNCTIONAL TRAITS of WARPING CURVES!!!!!!!!!!! 

#first need to have all of the arping curves as one functional object so will evaluate each one and combine 

new_time <- seq(56,152, length.out = 300)
BASIS_bspline <- fda::create.bspline.basis(rangeval = c(min(new_time),max(new_time)),norder = 4,nbasis = 300)
fdobj_outlier <- fda::fdPar(BASIS_bspline,Lfdobj = 2,lambda = 0.01)

#ori
#sdld
genotype_smooth_warp_oriSDLD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_oriSDLD[[x]]$warpfd)})
combined_df_oriSDLD <- bind_cols(genotype_smooth_warp_oriSDLD)
genotype_smooth_warp_oriSDLD <- combined_df_oriSDLD
colnames(genotype_smooth_warp_oriSDLD) <- colnames(Total_curves_SD_LD_rmout)
genotype_smooth_warp_oriSDLD <- data.matrix(genotype_smooth_warp_oriSDLD)

#sdsd
genotype_smooth_warp_oriSDSD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_oriSDSD[[x]]$warpfd)})
combined_df_oriSDSD <- bind_cols(genotype_smooth_warp_oriSDSD)
genotype_smooth_warp_oriSDSD <- combined_df_oriSDSD
colnames(genotype_smooth_warp_oriSDSD) <- colnames(Total_curves_SD_SD_rmout)
genotype_smooth_warp_oriSDSD <- data.matrix(genotype_smooth_warp_oriSDSD)

#ldld
genotype_smooth_warp_oriLDLD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_oriLDLD[[x]]$warpfd)})
combined_df_oriLDLD <- bind_cols(genotype_smooth_warp_oriLDLD)
genotype_smooth_warp_oriLDLD <- combined_df_oriLDLD
colnames(genotype_smooth_warp_oriLDLD) <- colnames(Total_curves_LD_LD_rmout)
genotype_smooth_warp_oriLDLD <- data.matrix(genotype_smooth_warp_oriLDLD)

#combine each of the dataframes to one 

ALL_genotype_smooth_warp_ori <- cbind(genotype_smooth_warp_oriSDLD,genotype_smooth_warp_oriSDSD,genotype_smooth_warp_oriLDLD)

#normalise all 
ALL_genotype_smooth_warp_ori <- apply(ALL_genotype_smooth_warp_ori,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))

smooth_warp_all_ori <- fda::smooth.basis(argvals = new_time, y = ALL_genotype_smooth_warp_ori, fdParobj = fdobj_outlier)
#####################################
#vel


#vel
#sdld
genotype_smooth_warp_velSDLD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_velSDLD[[x]]$warpfd)})
combined_df_velSDLD <- bind_cols(genotype_smooth_warp_velSDLD)
genotype_smooth_warp_velSDLD <- combined_df_velSDLD
colnames(genotype_smooth_warp_velSDLD) <- colnames(Total_curves_SD_LD_rmout)
genotype_smooth_warp_velSDLD <- data.matrix(genotype_smooth_warp_velSDLD)

#sdsd
genotype_smooth_warp_velSDSD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_velSDSD[[x]]$warpfd)})
combined_df_velSDSD <- bind_cols(genotype_smooth_warp_velSDSD)
genotype_smooth_warp_velSDSD <- combined_df_velSDSD
colnames(genotype_smooth_warp_velSDSD) <- colnames(Total_curves_SD_SD_rmout)
genotype_smooth_warp_velSDSD <- data.matrix(genotype_smooth_warp_velSDSD)

#ldld
genotype_smooth_warp_velLDLD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_velLDLD[[x]]$warpfd)})
combined_df_velLDLD <- bind_cols(genotype_smooth_warp_velLDLD)
genotype_smooth_warp_velLDLD <- combined_df_velLDLD
colnames(genotype_smooth_warp_velLDLD) <- colnames(Total_curves_LD_LD_rmout)
genotype_smooth_warp_velLDLD <- data.matrix(genotype_smooth_warp_velLDLD)

#combine each of the dataframes to one 

ALL_genotype_smooth_warp_vel <- cbind(genotype_smooth_warp_velSDLD,genotype_smooth_warp_velSDSD,genotype_smooth_warp_velLDLD)

#normalise all 
ALL_genotype_smooth_warp_vel <- apply(ALL_genotype_smooth_warp_vel,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))

smooth_warp_all_vel <- fda::smooth.basis(argvals = new_time, y = ALL_genotype_smooth_warp_vel, fdParobj = fdobj_outlier)


######################################
#acc

#acc
#sdld
genotype_smooth_warp_accSDLD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_accSDLD[[x]]$warpfd)})
combined_df_accSDLD <- bind_cols(genotype_smooth_warp_accSDLD)
genotype_smooth_warp_accSDLD <- combined_df_accSDLD
colnames(genotype_smooth_warp_accSDLD) <- colnames(Total_curves_SD_LD_rmout)
genotype_smooth_warp_accSDLD <- data.matrix(genotype_smooth_warp_accSDLD)

#sdsd
genotype_smooth_warp_accSDSD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_accSDSD[[x]]$warpfd)})
combined_df_accSDSD <- bind_cols(genotype_smooth_warp_accSDSD)
genotype_smooth_warp_accSDSD <- combined_df_accSDSD
colnames(genotype_smooth_warp_accSDSD) <- colnames(Total_curves_SD_SD_rmout)
genotype_smooth_warp_accSDSD <- data.matrix(genotype_smooth_warp_accSDSD)

#ldld
genotype_smooth_warp_accLDLD <- sapply(X = seq(1,64), FUN = function(x){eval.fd(new_time, reg_accLDLD[[x]]$warpfd)})
combined_df_accLDLD <- bind_cols(genotype_smooth_warp_accLDLD)
genotype_smooth_warp_accLDLD <- combined_df_accLDLD
colnames(genotype_smooth_warp_accLDLD) <- colnames(Total_curves_LD_LD_rmout)
genotype_smooth_warp_accLDLD <- data.matrix(genotype_smooth_warp_accLDLD)

#combine each of the dataframes to one 

ALL_genotype_smooth_warp_acc <- cbind(genotype_smooth_warp_accSDLD,genotype_smooth_warp_accSDSD,genotype_smooth_warp_accLDLD)

#normalise all 
ALL_genotype_smooth_warp_acc <- apply(ALL_genotype_smooth_warp_acc,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))

smooth_warp_all_acc <- fda::smooth.basis(argvals = new_time, y = ALL_genotype_smooth_warp_acc, fdParobj = fdobj_outlier)


##############################################################################################

###                                    TRAITS TO MEASURE 



#1                    FPCA of the origional 


#use the warping functions in the fpca 

#ori
fpca_ori <- pca.fd(smooth_warp_all_ori$fd,3,harmfdPar = fdobj_outlier,centerfns = FALSE)
scores <- fpca_ori$scores
#fpca_ori <- fda::varmx.pca.fd(fpca_ori)
#rownames(scores) <- colnames(Total_curves_SD_LD_rmout)
#this plots each of the principle compments
plot(fpca_ori$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_ori$varprop*100
sum(fpca_ori$varprop*100)



#vel
fpca_vel <- pca.fd(smooth_warp_all_vel$fd,3,harmfdPar = fdobj_outlier,centerfns = FALSE)
#fpca_vel <- fda::varmx.pca.fd(fpca_vel)
scores <- fpca_vel$scores
#rownames(scores) <- colnames(Total_curves_SD_LD_rmout)
#this plots each of the principle compments
plot(fpca_vel$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD velgional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_vel$varprop*100
sum(fpca_vel$varprop*100)


#acc
fpca_acc <- pca.fd(smooth_warp_all_acc$fd,3,harmfdPar = fdobj_outlier,centerfns = FALSE)
#fpca_acc <- fda::varmx.pca.fd(fpca_acc)
scores <- fpca_acc$scores
#rownames(scores) <- colnames(Total_curves_SD_LD_rmout)
#this plots each of the principle compments
plot(fpca_acc$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD accgional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_acc$varprop*100
sum(fpca_acc$varprop*100)


#make the my vetor twice as there is the SD and LD portion
my_vector_SDLD <- vector("numeric")
my_vector_SDSD <- vector("numeric")
my_vector_LDLD <- vector("numeric")

Genotype_names_sorted <- sort(Genotype_names)

for (i in Genotype_names_sorted) {
  my_vector_SDLD[i] <- length(Total_curves_SD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_SDSD[i] <- length(Total_curves_SD_SD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_LDLD[i] <- length(Total_curves_LD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
}


my_vector_SDLD <- my_vector_SDLD[order(names(my_vector_SDLD))]
my_vector_SDSD <- my_vector_SDSD[order(names(my_vector_SDSD))]
my_vector_LDLD <- my_vector_LDLD[order(names(my_vector_LDLD))]

my_vector <- c(my_vector_SDLD,my_vector_SDSD,my_vector_LDLD)
cum <- cumsum(my_vector)
genotype_list <- list()


center_coords_ori <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_ori[i,] <- c(mean(fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}
center_coords_vel <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_vel[i,] <- c(mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}

center_coords_acc <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_acc[i,] <- c(mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}


#plot all curve points and the. ccenters of the SD_LD, SD_SD and LD_LD curves 

plot(fpca_ori$scores,xlab = 'PC Score 1',ylab = 'PC Score 2',col = "grey",
     cex.lab = 1.5,cex.axis = 1.5,cex = 1)
points(center_coords_ori[1:64,1], center_coords_ori[1:64,2], col = "darkblue", pch = 19)
points(center_coords_ori[65:128,1], center_coords_ori[65:128,2], col = "darkred", pch = 19)
points(center_coords_ori[129:192,1], center_coords_ori[129:192,2], col = "darkgreen", pch = 19)
legend("bottomright", c("SD_LD", "SD_SD", "LD_LD"), col = c("darkblue", "darkred", "darkgreen"), pch = 19)



plot(fpca_vel$scores,xlab = 'PC Score 1',ylab = 'PC Score 2',col = "grey",
     cex.lab = 1.5,cex.axis = 1.5,cex = 1)
points(center_coords_vel[1:64,1], center_coords_vel[1:64,2], col = "darkblue", pch = 19)
points(center_coords_vel[65:128,1], center_coords_vel[65:128,2], col = "darkred", pch = 19)
points(center_coords_vel[129:192,1], center_coords_vel[129:192,2], col = "darkgreen", pch = 19)
legend("bottomright", c("SD_LD", "SD_SD", "LD_LD"), col = c("darkblue", "darkred", "darkgreen"), pch = 19)


plot(fpca_acc$scores,xlab = 'PC Score 1',ylab = 'PC Score 2',col = "grey",
     cex.lab = 1.5,cex.axis = 1.5,cex = 1)
points(center_coords_acc[1:64,1], center_coords_acc[1:64,2], col = "darkblue", pch = 19)
points(center_coords_acc[65:128,1], center_coords_acc[65:128,2], col = "darkred", pch = 19)
points(center_coords_acc[129:192,1], center_coords_acc[129:192,2], col = "darkgreen", pch = 19)
legend("bottomright", c("SD_LD", "SD_SD", "LD_LD"), col = c("darkblue", "darkred", "darkgreen"), pch = 19)

########
#calculate the spread looking at distance in pc1, pc2 and pc3 separatly its just the distance calulated by absolute

spread_coords_mean_ori <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_ori[i,] <- c((sum(abs(center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_vel <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_vel[i,] <- c((sum(abs(center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_acc <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_acc[i,] <- c((sum(abs(center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_mean_ori_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_ori_total[i,1] <- sum(sqrt(((center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_vel_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_vel_total[i,1] <- sum(sqrt(((center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_acc_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_acc_total[i,1] <- sum(sqrt(((center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}




####################### now the standard error of the spread 


spread_coords_se_ori <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_ori[i,] <- c((sd(abs(center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_vel <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_vel[i,] <- c((sd(abs(center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_acc <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_acc[i,] <- c((sd(abs(center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_se_ori_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_ori_total[i,1] <- sd(sqrt(((center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_vel_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_vel_total[i,1] <- sd(sqrt(((center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_acc_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_acc_total[i,1] <- sd(sqrt(((center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}



##########################################################################

#########################################################################


#make a table of traits with the ap, phase RSQR as traits. 

table_of_phenotypes <- cbind(decomp_oriSDLD[,1:3], decomp_velSDLD[,1:3], decomp_accSDLD[,1:3],decomp_oriSDSD[,1:3], decomp_velSDSD[,1:3], decomp_accSDSD[,1:3],decomp_oriLDLD[,1:3], decomp_velLDLD[,1:3], decomp_accLDLD[,1:3],diff_oriSDLD_SDSD[,1:3],diff_oriSDLD_LDLD[,1:3],diff_velSDLD_SDSD[,1:3],diff_velSDLD_LDLD[,1:3],diff_accSDLD_SDSD[,1:3],diff_accSDLD_LDLD[,1:3],mean_norm_dsitSDLD,logmean_norm_dsitSDLD,sd_norm_dsitSDLD,mean_norm_dsitSDSD,logmean_norm_dsitSDSD,sd_norm_dsitSDSD,mean_norm_dsitLDLD,logmean_norm_dsitLDLD,sd_norm_dsitLDLD,center_coords_ori[1:64,1],center_coords_ori[1:64,2],center_coords_ori[1:64,3],center_coords_ori[65:128,1],center_coords_ori[65:128,2],center_coords_ori[65:128,3],center_coords_ori[129:192,1],center_coords_ori[129:192,2],center_coords_ori[129:192,3],center_coords_vel[1:64,1],center_coords_vel[1:64,2],center_coords_vel[1:64,3],center_coords_vel[65:128,1],center_coords_vel[65:128,2],center_coords_vel[65:128,3],center_coords_vel[129:192,1],center_coords_vel[129:192,2],center_coords_vel[129:192,3],center_coords_acc[1:64,1],center_coords_acc[1:64,2],center_coords_acc[1:64,3],center_coords_acc[65:128,1],center_coords_acc[65:128,2],center_coords_acc[65:128,3],center_coords_acc[129:192,1],center_coords_acc[129:192,2],center_coords_acc[129:192,3], spread_coords_mean_ori[1:64,1], spread_coords_mean_ori[1:64,2],spread_coords_mean_ori[1:64,3],spread_coords_mean_ori[65:128,1], spread_coords_mean_ori[65:128,2],spread_coords_mean_ori[65:128,3],spread_coords_mean_ori[129:192,1], spread_coords_mean_ori[129:192,2],spread_coords_mean_ori[129:192,3],spread_coords_mean_vel[1:64,1], spread_coords_mean_vel[1:64,2],spread_coords_mean_vel[1:64,3],spread_coords_mean_vel[65:128,1], spread_coords_mean_vel[65:128,2],spread_coords_mean_vel[65:128,3],spread_coords_mean_vel[129:192,1], spread_coords_mean_vel[129:192,2],spread_coords_mean_vel[129:192,3],spread_coords_mean_acc[1:64,1], spread_coords_mean_acc[1:64,2],spread_coords_mean_acc[1:64,3],spread_coords_mean_acc[65:128,1], spread_coords_mean_acc[65:128,2],spread_coords_mean_acc[65:128,3],spread_coords_mean_acc[129:192,1], spread_coords_mean_acc[129:192,2],spread_coords_mean_acc[129:192,3],spread_coords_se_ori[1:64,1], spread_coords_se_ori[1:64,2],spread_coords_se_ori[1:64,3],spread_coords_se_ori[65:128,1], spread_coords_se_ori[65:128,2],spread_coords_se_ori[65:128,3],spread_coords_se_ori[129:192,1], spread_coords_se_ori[129:192,2],spread_coords_se_ori[129:192,3],spread_coords_se_vel[1:64,1], spread_coords_se_vel[1:64,2],spread_coords_se_vel[1:64,3],spread_coords_se_vel[65:128,1], spread_coords_se_vel[65:128,2],spread_coords_se_vel[65:128,3],spread_coords_se_vel[129:192,1], spread_coords_se_vel[129:192,2],spread_coords_se_vel[129:192,3],spread_coords_se_acc[1:64,1], spread_coords_se_acc[1:64,2],spread_coords_se_acc[1:64,3],spread_coords_se_acc[65:128,1], spread_coords_se_acc[65:128,2],spread_coords_se_acc[65:128,3],spread_coords_se_acc[129:192,1], spread_coords_se_acc[129:192,2],spread_coords_se_acc[129:192,3], spread_coords_mean_ori_total[1:64,1],spread_coords_mean_ori_total[65:128,1],spread_coords_mean_ori_total[129:192,1],spread_coords_mean_vel_total[1:64,1],spread_coords_mean_vel_total[65:128,1],spread_coords_mean_vel_total[129:192,1],spread_coords_mean_acc_total[1:64,1],spread_coords_mean_acc_total[65:128,1],spread_coords_mean_acc_total[129:192,1],spread_coords_se_ori_total[1:64,1],spread_coords_se_ori_total[65:128,1],spread_coords_se_ori_total[129:192,1],spread_coords_se_vel_total[1:64,1],spread_coords_se_vel_total[65:128,1],spread_coords_se_vel_total[129:192,1],spread_coords_se_acc_total[1:64,1],spread_coords_se_acc_total[65:128,1],spread_coords_se_acc_total[129:192,1])

row.names(table_of_phenotypes) <- sort(Genotype_names)


colnames(table_of_phenotypes) <- c("Amp_oriSDLD","Phase_oriSDLD", "RSQR_oriSDLD","Amp_velSDLD","Phase_velSDLD", "RSQR_velSDLD","Amp_accSDLD","Phase_accSDLD", "RSQR_accSDLD","Amp_oriSDSD","Phase_oriSDSD", "RSQR_oriSDSD","Amp_velSDSD","Phase_velSDSD", "RSQR_velSDSD","Amp_accSDSD","Phase_accSDSD", "RSQR_accSDSD","Amp_oriLDLD","Phase_oriLDLD", "RSQR_oriLDLD","Amp_velLDLD","Phase_velLDLD", "RSQR_velLDLD","Amp_accLDLD","Phase_accLDLD", "RSQR_accLDLD","Amp_diff_oriSDLD_SDSD","Phase_diff_oriSDLD_SDSD", "RSQR_diff_oriSDLD_SDSD","Amp_diff_oriSDLD_LDLD","Phase_diff_oriSDLD_LDLD", "RSQR_diff_oriSDLD_LDLD","Amp_diff_velSDLD_SDSD","Phase_diff_velSDLD_SDSD", "RSQR_diff_velSDLD_SDSD","Amp_diff_velSDLD_LDLD","Phase_diff_velSDLD_LDLD", "RSQR_diff_velSDLD_LDLD","Amp_diff_accSDLD_SDSD","Phase_diff_accSDLD_SDSD", "RSQR_diff_accSDLD_SDSD","Amp_diff_accSDLD_LDLD","Phase_diff_accSDLD_LDLD", "RSQR_diff_accSDLD_LDLD","dtw_norm_dist_meanSDLD","logdtw_norm_dist_meanSDLD", "dtw_norm_dist_sdSDLD","dtw_norm_dist_meanSDSD","logdtw_norm_dist_meanSDSD", "dtw_norm_dist_sdSDSD","dtw_norm_dist_meanLDLD","logdtw_norm_dist_meanLDLD", "dtw_norm_dist_sdLDLD","gc_warp_PC1_oriSDLD","gc_warp_PC2_oriSDLD","gc_warp_PC3_oriSDLD","gc_warp_PC1_oriSDSD","gc_warp_PC2_oriSDSD","gc_warp_PC3_oriSDSD", "gc_warp_PC1_oriLDLD","gc_warp_PC2_oriLDLD","gc_warp_PC3_oriLDLD","gc_warp_PC1_velSDLD","gc_warp_PC2_velSDLD","gc_warp_PC3_velSDLD","gc_warp_PC1_velSDSD","gc_warp_PC2_velSDSD","gc_warp_PC3_velSDSD", "gc_warp_PC1_velLDLD","gc_warp_PC2_velLDLD","gc_warp_PC3_velLDLD","gc_warp_PC1_accSDLD","gc_warp_PC2_accSDLD","gc_warp_PC3_accSDLD","gc_warp_PC1_accSDSD","gc_warp_PC2_accSDSD","gc_warp_PC3_accSDSD", "gc_warp_PC1_accLDLD","gc_warp_PC2_accLDLD","gc_warp_PC3_accLDLD","meandist_warp_PC1_oriSDLD","meandist_warp_PC2_oriSDLD","meandist_warp_PC3_oriSDLD","meandist_warp_PC1_oriSDSD","meandist_warp_PC2_oriSDSD","meandist_warp_PC3_oriSDSD","meandist_warp_PC1_oriLDLD","meandist_warp_PC2_oriLDLD","meandist_warp_PC3_oriLDLD","meandist_warp_PC1_velSDLD","meandist_warp_PC2_velSDLD","meandist_warp_PC3_velSDLD","meandist_warp_PC1_velSDSD","meandist_warp_PC2_velSDSD","meandist_warp_PC3_velSDSD","meandist_warp_PC1_velLDLD","meandist_warp_PC2_velLDLD","meandist_warp_PC3_velLDLD","meandist_warp_PC1_accSDLD","meandist_warp_PC2_accSDLD","meandist_warp_PC3_accSDLD","meandist_warp_PC1_accSDSD","meandist_warp_PC2_accSDSD","meandist_warp_PC3_accSDSD","meandist_warp_PC1_accLDLD","meandist_warp_PC2_accLDLD","meandist_warp_PC3_accLDLD","sedist_warp_PC1_oriSDLD","sedist_warp_PC2_oriSDLD","sedist_warp_PC3_oriSDLD","sedist_warp_PC1_oriSDSD","sedist_warp_PC2_oriSDSD","sedist_warp_PC3_oriSDSD","sedist_warp_PC1_oriLDLD","sedist_warp_PC2_oriLDLD","sedist_warp_PC3_oriLDLD","sedist_warp_PC1_velSDLD","sedist_warp_PC2_velSDLD","sedist_warp_PC3_velSDLD","sedist_warp_PC1_velSDSD","sedist_warp_PC2_velSDSD","sedist_warp_PC3_velSDSD","sedist_warp_PC1_velLDLD","sedist_warp_PC2_velLDLD","sedist_warp_PC3_velLDLD","sedist_warp_PC1_accSDLD","sedist_warp_PC2_accSDLD","sedist_warp_PC3_accSDLD","sedist_warp_PC1_accSDSD","sedist_warp_PC2_accSDSD","sedist_warp_PC3_accSDSD","sedist_warp_PC1_accLDLD","sedist_warp_PC2_accLDLD","sedist_warp_PC3_accLDLD","mean_eulid_ori_total123_SDLD","mean_eulid_ori_total123_SDSD","mean_eulid_ori_total123_LDLD","mean_eulid_vel_total123_SDLD","mean_eulid_vel_total123_SDSD","mean_eulid_vel_total123_LDLD","mean_eulid_acc_total123_SDLD","mean_eulid_acc_total123_SDSD","mean_eulid_acc_total123_LDLD","se_eulid_ori_total123_SDLD","se_eulid_ori_total123_SDSD","se_eulid_ori_total123_LDLD","se_eulid_vel_total123_SDLD","se_eulid_vel_total123_SDSD","se_eulid_vel_total123_LDLD","se_eulid_acc_total123_SDLD","se_eulid_acc_total123_SDSD","se_eulid_acc_total123_LDLD")


##
#save file
write.csv(table_of_phenotypes, "../TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/table_of_phenotypes_warping_full.csv")


###########################################################################


###   with varimax rotation!!! 




#ori
fpca_ori <- pca.fd(smooth_warp_all_ori$fd,3,harmfdPar = fdobj_outlier,centerfns = FALSE)
scores <- fpca_ori$scores
fpca_ori <- fda::varmx.pca.fd(fpca_ori)
#rownames(scores) <- colnames(Total_curves_SD_LD_rmout)
#this plots each of the principle compments
plot(fpca_ori$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_ori$varprop*100
sum(fpca_ori$varprop*100)



#vel
fpca_vel <- pca.fd(smooth_warp_all_vel$fd,3,harmfdPar = fdobj_outlier,centerfns = FALSE)
fpca_vel <- fda::varmx.pca.fd(fpca_vel)
scores <- fpca_vel$scores
#rownames(scores) <- colnames(Total_curves_SD_LD_rmout)
#this plots each of the principle compments
plot(fpca_vel$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD velgional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_vel$varprop*100
sum(fpca_vel$varprop*100)


#acc
fpca_acc <- pca.fd(smooth_warp_all_acc$fd,3,harmfdPar = fdobj_outlier,centerfns = FALSE)
fpca_acc <- fda::varmx.pca.fd(fpca_acc)
scores <- fpca_acc$scores
#rownames(scores) <- colnames(Total_curves_SD_LD_rmout)
#this plots each of the principle compments
plot(fpca_acc$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD accgional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_acc$varprop*100
sum(fpca_acc$varprop*100)


#make the my vetor twice as there is the SD and LD portion
my_vector_SDLD <- vector("numeric")
my_vector_SDSD <- vector("numeric")
my_vector_LDLD <- vector("numeric")


for (i in Genotype_names_sorted) {
  my_vector_SDLD[i] <- length(Total_curves_SD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_SDSD[i] <- length(Total_curves_SD_SD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_LDLD[i] <- length(Total_curves_LD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
}


my_vector_SDLD <- my_vector_SDLD[order(names(my_vector_SDLD))]
my_vector_SDSD <- my_vector_SDSD[order(names(my_vector_SDSD))]
my_vector_LDLD <- my_vector_LDLD[order(names(my_vector_LDLD))]

my_vector <- c(my_vector_SDLD,my_vector_SDSD,my_vector_LDLD)
cum <- cumsum(my_vector)
genotype_list <- list()


center_coords_ori <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_ori[i,] <- c(mean(fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}
center_coords_vel <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_vel[i,] <- c(mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}

center_coords_acc <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_acc[i,] <- c(mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}

########
#calculate the spread looking at distance in pc1, pc2 and pc3 separatly its just the distance calulated by absolute

spread_coords_mean_ori <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_ori[i,] <- c((sum(abs(center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_vel <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_vel[i,] <- c((sum(abs(center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_acc <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_acc[i,] <- c((sum(abs(center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_mean_ori_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_ori_total[i,1] <- sum(sqrt(((center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_vel_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_vel_total[i,1] <- sum(sqrt(((center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_acc_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_acc_total[i,1] <- sum(sqrt(((center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}


####################### now the standard error of the spread 


spread_coords_se_ori <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_ori[i,] <- c((sd(abs(center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_vel <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_vel[i,] <- c((sd(abs(center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_acc <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_acc[i,] <- c((sd(abs(center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_se_ori_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_ori_total[i,1] <- sd(sqrt(((center_coords_ori[i,1] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_ori[i,2] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_ori[i,3] - fpca_ori$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_vel_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_vel_total[i,1] <- sd(sqrt(((center_coords_vel[i,1] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_vel[i,2] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_vel[i,3] - fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_acc_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_acc_total[i,1] <- sd(sqrt(((center_coords_acc[i,1] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_acc[i,2] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_acc[i,3] - fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}





table_of_phenotypes_varmax <- cbind(center_coords_ori[1:64,1],center_coords_ori[1:64,2],center_coords_ori[1:64,3],center_coords_ori[65:128,1],center_coords_ori[65:128,2],center_coords_ori[65:128,3],center_coords_ori[129:192,1],center_coords_ori[129:192,2],center_coords_ori[129:192,3],center_coords_vel[1:64,1],center_coords_vel[1:64,2],center_coords_vel[1:64,3],center_coords_vel[65:128,1],center_coords_vel[65:128,2],center_coords_vel[65:128,3],center_coords_vel[129:192,1],center_coords_vel[129:192,2],center_coords_vel[129:192,3],center_coords_acc[1:64,1],center_coords_acc[1:64,2],center_coords_acc[1:64,3],center_coords_acc[65:128,1],center_coords_acc[65:128,2],center_coords_acc[65:128,3],center_coords_acc[129:192,1],center_coords_acc[129:192,2],center_coords_acc[129:192,3], spread_coords_mean_ori[1:64,1], spread_coords_mean_ori[1:64,2],spread_coords_mean_ori[1:64,3],spread_coords_mean_ori[65:128,1], spread_coords_mean_ori[65:128,2],spread_coords_mean_ori[65:128,3],spread_coords_mean_ori[129:192,1], spread_coords_mean_ori[129:192,2],spread_coords_mean_ori[129:192,3],spread_coords_mean_vel[1:64,1], spread_coords_mean_vel[1:64,2],spread_coords_mean_vel[1:64,3],spread_coords_mean_vel[65:128,1], spread_coords_mean_vel[65:128,2],spread_coords_mean_vel[65:128,3],spread_coords_mean_vel[129:192,1], spread_coords_mean_vel[129:192,2],spread_coords_mean_vel[129:192,3],spread_coords_mean_acc[1:64,1], spread_coords_mean_acc[1:64,2],spread_coords_mean_acc[1:64,3],spread_coords_mean_acc[65:128,1], spread_coords_mean_acc[65:128,2],spread_coords_mean_acc[65:128,3],spread_coords_mean_acc[129:192,1], spread_coords_mean_acc[129:192,2],spread_coords_mean_acc[129:192,3],spread_coords_se_ori[1:64,1], spread_coords_se_ori[1:64,2],spread_coords_se_ori[1:64,3],spread_coords_se_ori[65:128,1], spread_coords_se_ori[65:128,2],spread_coords_se_ori[65:128,3],spread_coords_se_ori[129:192,1], spread_coords_se_ori[129:192,2],spread_coords_se_ori[129:192,3],spread_coords_se_vel[1:64,1], spread_coords_se_vel[1:64,2],spread_coords_se_vel[1:64,3],spread_coords_se_vel[65:128,1], spread_coords_se_vel[65:128,2],spread_coords_se_vel[65:128,3],spread_coords_se_vel[129:192,1], spread_coords_se_vel[129:192,2],spread_coords_se_vel[129:192,3],spread_coords_se_acc[1:64,1], spread_coords_se_acc[1:64,2],spread_coords_se_acc[1:64,3],spread_coords_se_acc[65:128,1], spread_coords_se_acc[65:128,2],spread_coords_se_acc[65:128,3],spread_coords_se_acc[129:192,1], spread_coords_se_acc[129:192,2],spread_coords_se_acc[129:192,3], spread_coords_mean_ori_total[1:64,1],spread_coords_mean_ori_total[65:128,1],spread_coords_mean_ori_total[129:192,1],spread_coords_mean_vel_total[1:64,1],spread_coords_mean_vel_total[65:128,1],spread_coords_mean_vel_total[129:192,1],spread_coords_mean_acc_total[1:64,1],spread_coords_mean_acc_total[65:128,1],spread_coords_mean_acc_total[129:192,1],spread_coords_se_ori_total[1:64,1],spread_coords_se_ori_total[65:128,1],spread_coords_se_ori_total[129:192,1],spread_coords_se_vel_total[1:64,1],spread_coords_se_vel_total[65:128,1],spread_coords_se_vel_total[129:192,1],spread_coords_se_acc_total[1:64,1],spread_coords_se_acc_total[65:128,1],spread_coords_se_acc_total[129:192,1])


row.names(table_of_phenotypes_varmax) <- sort(Genotype_names)

colnames(table_of_phenotypes_varmax) <- c("varmax_gc_warp_PC1_oriSDLD","varmx_gc_warp_PC2_oriSDLD","varmx_gc_warp_PC3_oriSDLD","varmx_gc_warp_PC1_oriSDSD","varmx_gc_warp_PC2_oriSDSD","varmx_gc_warp_PC3_oriSDSD","varmx_gc_warp_PC1_oriLDLD","varmx_gc_warp_PC2_oriLDLD","varmx_gc_warp_PC3_oriLDLD","varmx_gc_warp_PC1_velSDLD","varmx_gc_warp_PC2_velSDLD","varmx_gc_warp_PC3_velSDLD","varmx_gc_warp_PC1_velSDSD","varmx_gc_warp_PC2_velSDSD","varmx_gc_warp_PC3_velSDSD","varmx_gc_warp_PC1_velLDLD","varmx_gc_warp_PC2_velLDLD","varmx_gc_warp_PC3_velLDLD","varmx_gc_warp_PC1_accSDLD","varmx_gc_warp_PC2_accSDLD","varmx_gc_warp_PC3_accSDLD","varmx_gc_warp_PC1_accSDSD","varmx_gc_warp_PC2_accSDSD","varmx_gc_warp_PC3_accSDSD","varmx_gc_warp_PC1_accLDLD","varmx_gc_warp_PC2_accLDLD","varmx_gc_warp_PC3_accLDLD","varmx_meandist_warp_PC1_oriSDLD","varmx_meandist_warp_PC2_oriSDLD","varmx_meandist_warp_PC3_oriSDLD","varmx_meandist_warp_PC1_oriSDSD","varmx_meandist_warp_PC2_oriSDSD","varmx_meandist_warp_PC3_oriSDSD","varmx_meandist_warp_PC1_oriLDLD","varmx_meandist_warp_PC2_oriLDLD","varmx_meandist_warp_PC3_oriLDLD","varmx_meandist_warp_PC1_velSDLD","varmx_meandist_warp_PC2_velSDLD","varmx_meandist_warp_PC3_velSDLD","varmx_meandist_warp_PC1_velSDSD","varmx_meandist_warp_PC2_velSDSD","varmx_meandist_warp_PC3_velSDSD","varmx_meandist_warp_PC1_velLDLD","varmx_meandist_warp_PC2_velLDLD","varmx_meandist_warp_PC3_velLDLD","varmx_meandist_warp_PC1_accSDLD","varmx_meandist_warp_PC2_accSDLD","varmx_meandist_warp_PC3_accSDLD","varmx_meandist_warp_PC1_accSDSD","varmx_meandist_warp_PC2_accSDSD","varmx_meandist_warp_PC3_accSDSD","varmx_meandist_warp_PC1_accLDLD","varmx_meandist_warp_PC2_accLDLD","varmx_meandist_warp_PC3_accLDLD","varmx_sedist_warp_PC1_oriSDLD","varmx_sedist_warp_PC2_oriSDLD","varmx_sedist_warp_PC3_oriSDLD","varmx_sedist_warp_PC1_oriSDSD","varmx_sedist_warp_PC2_oriSDSD","varmx_sedist_warp_PC3_oriSDSD","varmx_sedist_warp_PC1_oriLDLD","varmx_sedist_warp_PC2_oriLDLD","varmx_sedist_warp_PC3_oriLDLD","varmx_sedist_warp_PC1_velSDLD","varmx_sedist_warp_PC2_velSDLD","varmx_sedist_warp_PC3_velSDLD","varmx_sedist_warp_PC1_velSDSD","varmx_sedist_warp_PC2_velSDSD","varmx_sedist_warp_PC3_velSDSD","varmx_sedist_warp_PC1_velLDLD","varmx_sedist_warp_PC2_velLDLD","varmx_sedist_warp_PC3_velLDLD","varmx_sedist_warp_PC1_accSDLD","varmx_sedist_warp_PC2_accSDLD","varmx_sedist_warp_PC3_accSDLD","varmx_sedist_warp_PC1_accSDSD","varmx_sedist_warp_PC2_accSDSD","varmx_sedist_warp_PC3_accSDSD","varmx_sedist_warp_PC1_accLDLD","varmx_sedist_warp_PC2_accLDLD","varmx_sedist_warp_PC3_accLDLD","varmx_mean_eulid_ori_total123_SDLD","varmx_mean_eulid_ori_total123_SDSD","varmx_mean_eulid_ori_total123_LDLD","varmx_mean_eulid_vel_total123_SDLD","varmx_mean_eulid_vel_total123_SDSD","varmx_mean_eulid_vel_total123_LDLD","varmx_mean_eulid_acc_total123_SDLD","varmx_mean_eulid_acc_total123_SDSD","varmx_mean_eulid_acc_total123_LDLD","varmx_se_eulid_ori_total123_SDLD","varmx_se_eulid_ori_total123_SDSD","varmx_se_eulid_ori_total123_LDLD","varmx_se_eulid_vel_total123_SDLD","varmx_se_eulid_vel_total123_SDSD","varmx_se_eulid_vel_total123_LDLD","varmx_se_eulid_acc_total123_SDLD","varmx_se_eulid_acc_total123_SDSD","varmx_se_eulid_acc_total123_LDLD")

###


write.csv(table_of_phenotypes_varmax, "../TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/table_of_phenotypes_warping_full_varmax.csv")














