# FPCA of the sd curves of each genotype across all conditions. 

#need to estimate the curves separatly for each genotype and produce a sd curves .
#read in files 

Total_curves_SD_LD <- read.csv(file = "Total_curves_SDLD.csv")
Total_curves_SD_SD <- read.csv(file = "Total_curves_SDSD.csv")
Total_curves_LD_LD <- read.csv(file = "Total_curves_LDLD.csv")

# all curves were evalutaed under the same time frame so xan just combine the dataframes togetther 

new_ZT <- Total_curves_SD_LD[,1]

Total_curves_SD_LD <- Total_curves_SD_LD[,-1]
Total_curves_SD_SD <- Total_curves_SD_SD[,-1]
Total_curves_LD_LD <- Total_curves_LD_LD[,-1]

Total_curves_ALL_rmout <- cbind(Total_curves_SD_LD, Total_curves_SD_SD, Total_curves_LD_LD)

Total_curves_ALL_rmout <- apply(Total_curves_ALL_rmout,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))
Total_curves_ALL_rmout <- data.matrix(Total_curves_ALL_rmout)


Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

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
vel_genotypeSDLD <- lapply(X = seq(1,length(genotype_smoothSDLD)), FUN = function(x){deriv.fd(genotype_smoothSDLD[[x]]$fd,1)})
accel_genotypeSDLD <- lapply(X = seq(1,length(genotype_smoothSDLD)), FUN = function(x){deriv.fd(genotype_smoothSDLD[[x]]$fd,2)})

#on origional curves......
sd_genotype_oriSDLD <- lapply(X = seq(1,length(vel_genotypeSDLD)), FUN = function(x){sd.fd(genotype_smoothSDLD[[x]]$fd)})
#on velocity curves
sd_genotype_velSDLD <- lapply(X = seq(1,length(vel_genotypeSDLD)), FUN = function(x){sd.fd(vel_genotypeSDLD[[x]])})
#on accelaeration curves
sd_genotype_accSDLD <- lapply(X = seq(1,length(vel_genotypeSDLD)), FUN = function(x){sd.fd(accel_genotypeSDLD[[x]])})


###SDSD
vel_genotypeSDSD <- lapply(X = seq(1,length(genotype_smoothSDSD)), FUN = function(x){deriv.fd(genotype_smoothSDSD[[x]]$fd,1)})
accel_genotypeSDSD <- lapply(X = seq(1,length(genotype_smoothSDSD)), FUN = function(x){deriv.fd(genotype_smoothSDSD[[x]]$fd,2)})

#on origional curves......
sd_genotype_oriSDSD <- lapply(X = seq(1,length(vel_genotypeSDSD)), FUN = function(x){sd.fd(genotype_smoothSDSD[[x]]$fd)})
#on velocity curves
sd_genotype_velSDSD <- lapply(X = seq(1,length(vel_genotypeSDSD)), FUN = function(x){sd.fd(vel_genotypeSDSD[[x]])})
#on accelaeration curves
sd_genotype_accSDSD <- lapply(X = seq(1,length(vel_genotypeSDSD)), FUN = function(x){sd.fd(accel_genotypeSDSD[[x]])})


###LDLD
vel_genotypeLDLD <- lapply(X = seq(1,length(genotype_smoothLDLD)), FUN = function(x){deriv.fd(genotype_smoothLDLD[[x]]$fd,1)})
accel_genotypeLDLD <- lapply(X = seq(1,length(genotype_smoothLDLD)), FUN = function(x){deriv.fd(genotype_smoothLDLD[[x]]$fd,2)})

#on origional curves......
sd_genotype_oriLDLD <- lapply(X = seq(1,length(vel_genotypeLDLD)), FUN = function(x){sd.fd(genotype_smoothLDLD[[x]]$fd)})
#on velocity curves
sd_genotype_velLDLD <- lapply(X = seq(1,length(vel_genotypeLDLD)), FUN = function(x){sd.fd(vel_genotypeLDLD[[x]])})
#on accelaeration curves
sd_genotype_accLDLD <- lapply(X = seq(1,length(vel_genotypeLDLD)), FUN = function(x){sd.fd(accel_genotypeLDLD[[x]])})

##### CONTROLS TOGETHER 

vel_genotypecontrols <- lapply(X = seq(1,length(genotype_smoothcontrols)), FUN = function(x){deriv.fd(genotype_smoothcontrols[[x]]$fd,1)})
accel_genotypecontrols <- lapply(X = seq(1,length(genotype_smoothcontrols)), FUN = function(x){deriv.fd(genotype_smoothcontrols[[x]]$fd,2)})

#on origional curves......
sd_genotype_oricontrols <- lapply(X = seq(1,length(vel_genotypecontrols)), FUN = function(x){sd.fd(genotype_smoothcontrols[[x]]$fd)})
#on velocity curves
sd_genotype_velcontrols <- lapply(X = seq(1,length(vel_genotypecontrols)), FUN = function(x){sd.fd(vel_genotypecontrols[[x]])})
#on accelaeration curves
sd_genotype_acccontrols <- lapply(X = seq(1,length(vel_genotypecontrols)), FUN = function(x){sd.fd(accel_genotypecontrols[[x]])})



# need to evalulate the sd curves for the pre and post shift 

#FIRST ------- evaluate the curves from 48 hours before the missed cue..

time_48_to_96 <- seq(56,104, length.out = 100 )

SD_48_to_96SDLD <-  sapply(X = seq(1,length(vel_genotypeSDLD)), FUN = function(x){eval.fd(time_48_to_96, sd_genotype_velSDLD[[x]])})
colnames(SD_48_to_96SDLD) <- Genotype_names

SD_48_to_96SDSD <-  sapply(X = seq(1,length(vel_genotypeSDSD)), FUN = function(x){eval.fd(time_48_to_96, sd_genotype_velSDSD[[x]])})
colnames(SD_48_to_96SDSD) <- Genotype_names

SD_48_to_96LDLD <-  sapply(X = seq(1,length(vel_genotypeLDLD)), FUN = function(x){eval.fd(time_48_to_96, sd_genotype_velLDLD[[x]])})
colnames(SD_48_to_96LDLD) <- Genotype_names

SD_48_to_96controls <-  sapply(X = seq(1,length(vel_genotypecontrols)), FUN = function(x){eval.fd(time_48_to_96, sd_genotype_velcontrols[[x]])})
colnames(SD_48_to_96controls) <- Genotype_names

#SECOND ------- evaluate the curves from 48 hours AFTER the missed cue..

time_96_to_144 <- seq(104,152, length.out = 100 )

LD_96_to_144SDLD <- sapply(X = seq(1,length(vel_genotypeSDLD)), FUN = function(x){eval.fd(evalarg = time_96_to_144, sd_genotype_velSDLD[[x]])})
colnames(LD_96_to_144SDLD) <- Genotype_names

LD_96_to_144SDSD <- sapply(X = seq(1,length(vel_genotypeSDSD)), FUN = function(x){eval.fd(evalarg = time_96_to_144, sd_genotype_velSDSD[[x]])})
colnames(LD_96_to_144SDSD) <- Genotype_names

LD_96_to_144LDLD <- sapply(X = seq(1,length(vel_genotypeLDLD)), FUN = function(x){eval.fd(evalarg = time_96_to_144, sd_genotype_velLDLD[[x]])})
colnames(LD_96_to_144LDLD) <- Genotype_names

LD_96_to_144controls <- sapply(X = seq(1,length(vel_genotypecontrols)), FUN = function(x){eval.fd(evalarg = time_96_to_144, sd_genotype_velcontrols[[x]])})
colnames(LD_96_to_144controls) <- Genotype_names


#combine all should be two per genotype....

curves_pre_post <- cbind(SD_48_to_96SDLD,LD_96_to_144SDLD, SD_48_to_96SDSD,LD_96_to_144SDSD, SD_48_to_96LDLD,LD_96_to_144LDLD)


#if using controls together 
#curves_pre_post <- cbind(SD_48_to_96SDLD,LD_96_to_144SDLD, SD_48_to_96controls,LD_96_to_144controls)

##### so LDLD has no WxT2 so as only one curve so remove if using both controls separatly! 

curves_pre_post[,268]
curves_pre_post[,332]

curves_pre_post <- curves_pre_post[,- c(268,332)]


# now estimate the combined curves together..

time_both <- seq(1,48, length.out = 100 )

BASIS_bspline <- create.bspline.basis(rangeval = c(1,48),norder = 4,nbasis = 102)
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 10)
smooth_both_ori <- smooth.basis(argvals = time_both, y = curves_pre_post, fdParobj = fdobj)
#plot(smooth_both_ori)

smooth_both_vel <- smooth.basis(argvals = time_both, y = curves_pre_post, fdParobj = fdobj)
smooth_both_acc <- smooth.basis(argvals = time_both, y = curves_pre_post, fdParobj = fdobj)

ori_obj <- deriv.fd(smooth_both_ori$fd,0)
vel_obj <- deriv.fd(smooth_both_vel$fd,1)
acc_obj <- deriv.fd(smooth_both_acc$fd,2)



#FPCA
#ori
fpca_ori <- pca.fd(ori_obj,3,harmfdPar = fdPar(ori_obj),centerfns = FALSE)
scores_ori <- fpca_ori$scores
rownames(scores_ori) <- colnames(curves_pre_post)
#this plots each of the principle compments
plot(fpca_ori$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_ori$varprop*100
sum(fpca_ori$varprop*100)

#vel
fpca_vel <- pca.fd(vel_obj,3,harmfdPar = fdPar(vel_obj),centerfns = FALSE)
scores_vel <- fpca_vel$scores
rownames(scores_vel) <- colnames(curves_pre_post)
#this plots each of the principle compments
plot(fpca_vel$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD velocity curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_vel$varprop*100
sum(fpca_vel$varprop*100)

#acc
fpca_acc <- pca.fd(acc_obj,2,harmfdPar = fdPar(acc_obj),centerfns = FALSE)
scores_acc <- fpca_acc$scores
rownames(scores_acc) <- colnames(curves_pre_post)
#this plots each of the principle compments
plot(fpca_acc$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_acc$varprop*100
sum(fpca_acc$varprop*100)

#look at where each genotype is on PC2

plot(fpca_ori$scores[,1], fpca_ori$scores[,2], main = "ori")
abline(h=0)
points(fpca_ori$scores[1:64,1],fpca_ori$scores[1:64,2], pch = 19, col = "slateblue2")
points(fpca_ori$scores[65:128,1],fpca_ori$scores[65:128,2], pch = 19, col = "darkblue")

points(fpca_ori$scores[129:192,1],fpca_ori$scores[129:192,2], pch = 19, col = "red2")
points(fpca_ori$scores[193:256,1],fpca_ori$scores[193:256,2], pch = 19, col = "darkred")

points(fpca_ori$scores[257:319,1],fpca_ori$scores[257:319,2], pch = 19, col = "green3")
points(fpca_ori$scores[320:382,1],fpca_ori$scores[320:382,2], pch = 19, col = "darkgreen")


legend("topright",  c("SD_LD 48-96 h","SD_LD 96-144 h","SD_SD 48-96 h","SD_SD 96-144 h","LD_LD 48-96 h","LD_LD 96-144 h"), col =  c("slateblue2","darkblue","red2","darkred","green3","darkgreen"), pch = 19)

plot(fpca_vel$scores[,1], fpca_vel$scores[,2], main = "vel")
abline(h=0)

points(fpca_vel$scores[1:64,1],fpca_vel$scores[1:64,2], pch = 19, col = "slateblue2")
points(fpca_vel$scores[65:128,1],fpca_vel$scores[65:128,2], pch = 19, col = "darkblue")

points(fpca_vel$scores[129:192,1],fpca_vel$scores[129:192,2], pch = 19, col = "red2")
points(fpca_vel$scores[193:256,1],fpca_vel$scores[193:256,2], pch = 19, col = "darkred")

points(fpca_vel$scores[257:319,1],fpca_vel$scores[257:319,2], pch = 19, col = "green3")
points(fpca_vel$scores[320:382,1],fpca_vel$scores[320:382,2], pch = 19, col = "darkgreen")



plot(fpca_acc$scores[,1], fpca_acc$scores[,2], main = "acc")
abline(h=0)
points(fpca_acc$scores[1:64,1],fpca_acc$scores[1:64,2], pch = 19, col = "slateblue2")
points(fpca_acc$scores[65:128,1],fpca_acc$scores[65:128,2], pch = 19, col = "darkblue")

points(fpca_acc$scores[129:192,1],fpca_acc$scores[129:192,2], pch = 19, col = "red2")
points(fpca_acc$scores[193:256,1],fpca_acc$scores[193:256,2], pch = 19, col = "darkred")

points(fpca_acc$scores[257:319,1],fpca_acc$scores[257:319,2], pch = 19, col = "green3")
points(fpca_acc$scores[320:382,1],fpca_acc$scores[320:382,2], pch = 19, col = "darkgreen")



#group the controls together... 
#1:64
#65:128
#129:192
#192:256
#257:320
#321:384

#look at the movement in pc1 and pc2 from pre to post 
PC1_movement_ori_sd_SDLD <- fpca_ori$scores[1:64,1] - fpca_ori$scores[65:128,1] 
PC2_movement_ori_sd_SDLD <- fpca_ori$scores[1:64,2] - fpca_ori$scores[65:128,2] 

PC1_movement_ori_sd_SDSD <- fpca_ori$scores[129:192,1] - fpca_ori$scores[193:256,1] 
PC2_movement_ori_sd_SDSD <- fpca_ori$scores[129:192,2] - fpca_ori$scores[193:256,2] 

PC1_movement_ori_sd_LDLD <- fpca_ori$scores[257:319,1] - fpca_ori$scores[320:382,1] 
PC2_movement_ori_sd_LDLD <- fpca_ori$scores[257:319,2] - fpca_ori$scores[320:382,2] 

#VEL
PC1_movement_vel_sd_SDLD <- fpca_vel$scores[1:64,1] - fpca_vel$scores[65:128,1] 
PC2_movement_vel_sd_SDLD <- fpca_vel$scores[1:64,2] - fpca_vel$scores[65:128,2] 

PC1_movement_vel_sd_SDSD <- fpca_vel$scores[129:192,1] - fpca_vel$scores[193:256,1] 
PC2_movement_vel_sd_SDSD <- fpca_vel$scores[129:192,2] - fpca_vel$scores[193:256,2] 

PC1_movement_vel_sd_LDLD <- fpca_vel$scores[257:319,1] - fpca_vel$scores[320:382,1] 
PC2_movement_vel_sd_LDLD <- fpca_vel$scores[257:319,2] - fpca_vel$scores[320:382,2] 

#ACC
PC1_movement_acc_sd_SDLD <- fpca_acc$scores[1:64,1] - fpca_acc$scores[65:128,1] 
PC2_movement_acc_sd_SDLD <- fpca_acc$scores[1:64,2] - fpca_acc$scores[65:128,2] 

PC1_movement_acc_sd_SDSD <- fpca_acc$scores[129:192,1] - fpca_acc$scores[193:256,1] 
PC2_movement_acc_sd_SDSD <- fpca_acc$scores[129:192,2] - fpca_acc$scores[193:256,2] 

PC1_movement_acc_sd_LDLD <- fpca_acc$scores[257:319,1] - fpca_acc$scores[320:382,1] 
PC2_movement_acc_sd_LDLD <- fpca_acc$scores[257:319,2] - fpca_acc$scores[320:382,2] 



#also look at the absolute values of this 
PC1_movement_ori_sd_absSDLD <- abs(PC1_movement_ori_sd_SDLD)
PC2_movement_ori_sd_absSDLD <- abs(PC2_movement_ori_sd_SDLD)

PC1_movement_ori_sd_absSDSD <- abs(PC1_movement_ori_sd_SDSD)
PC2_movement_ori_sd_absSDSD <- abs(PC2_movement_ori_sd_SDSD)

PC1_movement_ori_sd_absLDLD <- abs(PC1_movement_ori_sd_LDLD)
PC2_movement_ori_sd_absLDLD <- abs(PC2_movement_ori_sd_LDLD)

#vel
PC1_movement_vel_sd_absSDLD <- abs(PC1_movement_vel_sd_SDLD)
PC2_movement_vel_sd_absSDLD <- abs(PC2_movement_vel_sd_SDLD)

PC1_movement_vel_sd_absSDSD <- abs(PC1_movement_vel_sd_SDSD)
PC2_movement_vel_sd_absSDSD <- abs(PC2_movement_vel_sd_SDSD)

PC1_movement_vel_sd_absLDLD <- abs(PC1_movement_vel_sd_LDLD)
PC2_movement_vel_sd_absLDLD <- abs(PC2_movement_vel_sd_LDLD)

#acc

PC1_movement_acc_sd_absSDLD <- abs(PC1_movement_acc_sd_SDLD)
PC2_movement_acc_sd_absSDLD <- abs(PC2_movement_acc_sd_SDLD)

PC1_movement_acc_sd_absSDSD <- abs(PC1_movement_acc_sd_SDSD)
PC2_movement_acc_sd_absSDSD <- abs(PC2_movement_acc_sd_SDSD)

PC1_movement_acc_sd_absLDLD <- abs(PC1_movement_acc_sd_LDLD)
PC2_movement_acc_sd_absLDLD <- abs(PC2_movement_acc_sd_LDLD)



#euclid distance between pre and post 

#sdld
distances_ori_sdSDLD <- sqrt(rowSums((fpca_ori$scores[1:64, 1:2] - fpca_ori$scores[65:128, 1:2])^2))
distances_vel_sdSDLD <- sqrt(rowSums((fpca_vel$scores[1:64, 1:2] - fpca_vel$scores[65:128, 1:2])^2))
distances_acc_sdSDLD <- sqrt(rowSums((fpca_acc$scores[1:64, 1:2] - fpca_acc$scores[65:128, 1:2])^2))

#sdsd
distances_ori_sdSDSD <- sqrt(rowSums((fpca_ori$scores[129:192, 1:2] - fpca_ori$scores[193:256, 1:2])^2))
distances_vel_sdSDSD <- sqrt(rowSums((fpca_vel$scores[129:192, 1:2] - fpca_vel$scores[193:256, 1:2])^2))
distances_acc_sdSDSD <- sqrt(rowSums((fpca_acc$scores[129:192, 1:2] - fpca_acc$scores[193:256, 1:2])^2))

#ldld
distances_ori_sdLDLD <- sqrt(rowSums((fpca_ori$scores[257:319, 1:2] - fpca_ori$scores[320:382, 1:2])^2))
distances_vel_sdLDLD <- sqrt(rowSums((fpca_vel$scores[257:319, 1:2] - fpca_vel$scores[320:382, 1:2])^2))
distances_acc_sdLDLD <- sqrt(rowSums((fpca_acc$scores[257:319, 1:2] - fpca_acc$scores[320:382, 1:2])^2))


#make a table of each of these.... 




table_of_phenotypes <- cbind(fpca_ori$scores[1:64,1],fpca_ori$scores[65:128,1],fpca_ori$scores[129:192,1],fpca_ori$scores[193:256,1],fpca_ori$scores[1:64,2],fpca_ori$scores[65:128,2],fpca_ori$scores[129:192,2],fpca_ori$scores[193:256,2],fpca_vel$scores[1:64,1],fpca_vel$scores[65:128,1],fpca_vel$scores[129:192,1],fpca_vel$scores[193:256,1],fpca_vel$scores[1:64,2],fpca_vel$scores[65:128,2],fpca_vel$scores[129:192,2],fpca_vel$scores[193:256,2],fpca_acc$scores[1:64,1],fpca_acc$scores[65:128,1],fpca_acc$scores[129:192,1],fpca_acc$scores[193:256,1],fpca_acc$scores[1:64,2],fpca_acc$scores[65:128,2],fpca_acc$scores[129:192,2],fpca_acc$scores[193:256,2],PC1_movement_ori_sd_SDLD, PC2_movement_ori_sd_SDLD, PC1_movement_ori_sd_SDSD,PC2_movement_ori_sd_SDSD ,PC1_movement_vel_sd_SDLD,PC2_movement_vel_sd_SDLD ,PC1_movement_vel_sd_SDSD ,PC2_movement_vel_sd_SDSD,PC1_movement_acc_sd_SDLD ,PC2_movement_acc_sd_SDLD ,PC1_movement_acc_sd_SDSD,PC2_movement_acc_sd_SDSD,PC1_movement_ori_sd_absSDLD,PC2_movement_ori_sd_absSDLD,PC1_movement_ori_sd_absSDSD,PC2_movement_ori_sd_absSDSD,PC1_movement_vel_sd_absSDLD,PC2_movement_vel_sd_absSDLD,PC1_movement_vel_sd_absSDSD,PC2_movement_vel_sd_absSDSD,PC1_movement_acc_sd_absSDLD,PC2_movement_acc_sd_absSDLD,PC1_movement_acc_sd_absSDSD,PC2_movement_acc_sd_absSDSD,distances_ori_sdSDLD,distances_vel_sdSDLD,distances_acc_sdSDLD,distances_ori_sdSDSD,distances_vel_sdSDSD,distances_acc_sdSDSD)

#
row.names(table_of_phenotypes) <- Genotype_names

colnames(table_of_phenotypes) <- c("sd_PC1_SDLDpre","sd_PC1_SDLDpost","sd_PC1_SDSDpre", "sd_PC1_SDSDpost","sd_PC2_SDLDpre","sd_PC2_SDLDpost","sd_PC2_SDSDpre", "sd_PC2_SDSDpost", "sd_PC1_SDLDpre_vel","sd_PC1_SDLDpost_vel","sd_PC1_SDSDpre_vel", "sd_PC1_SDSDpost_vel","sd_PC2_SDLDpre_vel","sd_PC2_SDLDpost_vel","sd_PC2_SDSDpre_vel", "sd_PC2_SDSDpost_vel","sd_PC1_SDLDpre_acc","sd_PC1_SDLDpost_acc","sd_PC1_SDSDpre_acc", "sd_PC1_SDSDpost_acc","sd_PC2_SDLDpre_acc","sd_PC2_SDLDpost_acc","sd_PC2_SDSDpre_acc", "sd_PC2_SDSDpost_acc","PC1_movement_ori_sd_SDLD"," PC2_movement_ori_sd_SDLD"," PC1_movement_ori_sd_SDSD","PC2_movement_ori_sd_SDSD ","PC1_movement_vel_sd_SDLD","PC2_movement_vel_sd_SDLD ","PC1_movement_vel_sd_SDSD ","PC2_movement_vel_sd_SDSD","PC1_movement_acc_sd_SDLD ","PC2_movement_acc_sd_SDLD ","PC1_movement_acc_sd_SDSD","PC2_movement_acc_sd_SDSD","PC1_movement_ori_sd_absSDLD","PC2_movement_ori_sd_absSDLD","PC1_movement_ori_sd_absSDSD","PC2_movement_ori_sd_absSDSD","PC1_movement_vel_sd_absSDLD","PC2_movement_vel_sd_absSDLD","PC1_movement_vel_sd_absSDSD","PC2_movement_vel_sd_absSDSD","PC1_movement_acc_sd_absSDLD","PC2_movement_acc_sd_absSDLD","PC1_movement_acc_sd_absSDSD","PC2_movement_acc_sd_absSDSD","distances_ori_sdSDLD","distances_vel_sdSDLD","distances_acc_sdSDLD","distances_ori_sdSDSD","distances_vel_sdSDSD","distances_acc_sdSDSD")

#because LDLD had a genotype missing (WxT2 do a separate table) can add in excel and just have it as

#add extra row to the LDLD sections 
table_of_phenotypesLDLD <- cbind(fpca_ori$scores[257:319,1],fpca_ori$scores[320:382,1],fpca_ori$scores[257:319,2],fpca_ori$scores[320:382,2],fpca_vel$scores[257:319,1],fpca_vel$scores[320:382,1],fpca_vel$scores[257:319,2],fpca_vel$scores[320:382,2],fpca_acc$scores[257:319,1],fpca_acc$scores[320:382,1],fpca_acc$scores[257:319,2],fpca_acc$scores[320:382,2],PC1_movement_acc_sd_LDLD ,PC2_movement_acc_sd_LDLD ,PC1_movement_ori_sd_LDLD,PC2_movement_ori_sd_LDLD,PC1_movement_vel_sd_LDLD ,PC2_movement_vel_sd_LDLD,PC1_movement_ori_sd_absLDLD,PC2_movement_ori_sd_absLDLD,PC1_movement_vel_sd_absLDLD,PC2_movement_vel_sd_absLDLD,PC1_movement_acc_sd_absLDLD,PC2_movement_acc_sd_absLDLD,distances_ori_sdLDLD,distances_vel_sdLDLD,distances_acc_sdLDLD)


colnames(table_of_phenotypesLDLD) <- c("sd_PC1_LDLDpre","sd_PC1_LDLDpost","sd_PC2_LDLDpre","sd_PC2_LDLDpost","sd_PC1_LDLDpre_vel","sd_PC1_LDLDpost_vel","sd_PC2_LDLDpre_vel","sd_PC2_LDLDpost_vel","sd_PC1_LDLDpre_acc","sd_PC1_LDLDpost_acc","sd_PC2_LDLDpre_acc","sd_PC2_LDLDpost_acc","PC1_movement_acc_sd_LDLD","PC2_movement_acc_sd_LDLD","PC1_movement_ori_sd_LDLD","PC2_movement_ori_sd_LDLD","PC1_movement_vel_sd_LDLD","PC2_movement_vel_sd_LDLD","PC1_movement_ori_sd_absLDLD","PC2_movement_ori_sd_absLDLD","PC1_movement_vel_sd_absLDLD","PC2_movement_vel_sd_absLDLD","PC1_movement_acc_sd_absLDLD","PC2_movement_acc_sd_absLDLD","distances_ori_sdLDLD","distances_vel_sdLDLD","distances_acc_sdLDLD")


sub_genotypes <- Genotype_names[-12]  
row.names(table_of_phenotypesLDLD) <- sub_genotypes


##
#save file
write.csv(table_of_phenotypes, "table_of_phenotypes_pre_post_sdcurves.csv")

write.csv(table_of_phenotypesLDLD, "table_of_phenotypes_pre_post_sdcurvesLDLD.csv")



###########################################################################


# same but with varimx rotation 



#FPCA
#ori
fpca_ori <- pca.fd(ori_obj,3,harmfdPar = fdPar(ori_obj),centerfns = FALSE)
fpca_ori <- fda::varmx.pca.fd(fpca_ori)
scores_ori <- fpca_ori$scores
rownames(scores_ori) <- colnames(curves_pre_post)
#this plots each of the principle compments
plot(fpca_ori$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_ori$varprop*100
sum(fpca_ori$varprop*100)

#vel
fpca_vel <- pca.fd(vel_obj,3,harmfdPar = fdPar(vel_obj),centerfns = FALSE)
fpca_vel <- fda::varmx.pca.fd(fpca_vel)
scores_vel <- fpca_vel$scores
rownames(scores_vel) <- colnames(curves_pre_post)
#this plots each of the principle compments
plot(fpca_vel$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD velocity curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_vel$varprop*100
sum(fpca_vel$varprop*100)

#acc
fpca_acc <- pca.fd(acc_obj,3,harmfdPar = fdPar(acc_obj),centerfns = FALSE)
fpca_acc <- fda::varmx.pca.fd(fpca_acc)
scores_acc <- fpca_acc$scores
rownames(scores_acc) <- colnames(curves_pre_post)
#this plots each of the principle compments
plot(fpca_acc$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
fpca_acc$varprop*100
sum(fpca_acc$varprop*100)

#look at where each genotype is on PC2

plot(fpca_ori$scores[,1], fpca_ori$scores[,2], main = "ori")
abline(h=0)
points(fpca_ori$scores[1:64,1],fpca_ori$scores[1:64,2], pch = 19, col = "slateblue2")
points(fpca_ori$scores[65:128,1],fpca_ori$scores[65:128,2], pch = 19, col = "darkblue")

points(fpca_ori$scores[129:192,1],fpca_ori$scores[129:192,2], pch = 19, col = "red2")
points(fpca_ori$scores[193:256,1],fpca_ori$scores[193:256,2], pch = 19, col = "darkred")

points(fpca_ori$scores[257:319,1],fpca_ori$scores[257:319,2], pch = 19, col = "green3")
points(fpca_ori$scores[320:382,1],fpca_ori$scores[320:382,2], pch = 19, col = "darkgreen")


legend("topright",  c("SD_LD 48-96 h","SD_LD 96-144 h","SD_SD 48-96 h","SD_SD 96-144 h","LD_LD 48-96 h","LD_LD 96-144 h"), col =  c("slateblue2","darkblue","red2","darkred","green3","darkgreen"), pch = 19)



plot(fpca_vel$scores[,1], fpca_vel$scores[,2], main = "vel")
abline(h=0)

points(fpca_vel$scores[1:64,1],fpca_vel$scores[1:64,2], pch = 19, col = "slateblue2")
points(fpca_vel$scores[65:128,1],fpca_vel$scores[65:128,2], pch = 19, col = "darkblue")

points(fpca_vel$scores[129:192,1],fpca_vel$scores[129:192,2], pch = 19, col = "red2")
points(fpca_vel$scores[193:256,1],fpca_vel$scores[193:256,2], pch = 19, col = "darkred")

points(fpca_vel$scores[257:319,1],fpca_vel$scores[257:319,2], pch = 19, col = "green3")
points(fpca_vel$scores[320:382,1],fpca_vel$scores[320:382,2], pch = 19, col = "darkgreen")



plot(fpca_acc$scores[,1], fpca_acc$scores[,2], main = "acc")
abline(h=0)
points(fpca_acc$scores[1:64,1],fpca_acc$scores[1:64,2], pch = 19, col = "slateblue2")
points(fpca_acc$scores[65:128,1],fpca_acc$scores[65:128,2], pch = 19, col = "darkblue")

points(fpca_acc$scores[129:192,1],fpca_acc$scores[129:192,2], pch = 19, col = "red2")
points(fpca_acc$scores[193:256,1],fpca_acc$scores[193:256,2], pch = 19, col = "darkred")

points(fpca_acc$scores[257:319,1],fpca_acc$scores[257:319,2], pch = 19, col = "green3")
points(fpca_acc$scores[320:382,1],fpca_acc$scores[320:382,2], pch = 19, col = "darkgreen")


#could look at the movement in pc1 and pc2 from pre to post 
PC1_movement_ori_sd_SDLD <- fpca_ori$scores[1:64,1] - fpca_ori$scores[65:128,1] 
PC2_movement_ori_sd_SDLD <- fpca_ori$scores[1:64,2] - fpca_ori$scores[65:128,2] 

PC1_movement_ori_sd_SDSD <- fpca_ori$scores[129:192,1] - fpca_ori$scores[193:256,1] 
PC2_movement_ori_sd_SDSD <- fpca_ori$scores[129:192,2] - fpca_ori$scores[193:256,2] 

PC1_movement_ori_sd_LDLD <- fpca_ori$scores[257:319,1] - fpca_ori$scores[320:382,1] 
PC2_movement_ori_sd_LDLD <- fpca_ori$scores[257:319,2] - fpca_ori$scores[320:382,2] 

#VEL
PC1_movement_vel_sd_SDLD <- fpca_vel$scores[1:64,1] - fpca_vel$scores[65:128,1] 
PC2_movement_vel_sd_SDLD <- fpca_vel$scores[1:64,2] - fpca_vel$scores[65:128,2] 

PC1_movement_vel_sd_SDSD <- fpca_vel$scores[129:192,1] - fpca_vel$scores[193:256,1] 
PC2_movement_vel_sd_SDSD <- fpca_vel$scores[129:192,2] - fpca_vel$scores[193:256,2] 

PC1_movement_vel_sd_LDLD <- fpca_vel$scores[257:319,1] - fpca_vel$scores[320:382,1] 
PC2_movement_vel_sd_LDLD <- fpca_vel$scores[257:319,2] - fpca_vel$scores[320:382,2] 

#ACC
PC1_movement_acc_sd_SDLD <- fpca_acc$scores[1:64,1] - fpca_acc$scores[65:128,1] 
PC2_movement_acc_sd_SDLD <- fpca_acc$scores[1:64,2] - fpca_acc$scores[65:128,2] 

PC1_movement_acc_sd_SDSD <- fpca_acc$scores[129:192,1] - fpca_acc$scores[193:256,1] 
PC2_movement_acc_sd_SDSD <- fpca_acc$scores[129:192,2] - fpca_acc$scores[193:256,2] 

PC1_movement_acc_sd_LDLD <- fpca_acc$scores[257:319,1] - fpca_acc$scores[320:382,1] 
PC2_movement_acc_sd_LDLD <- fpca_acc$scores[257:319,2] - fpca_acc$scores[320:382,2] 



#also look at the absolute values of this 
PC1_movement_ori_sd_absSDLD <- abs(PC1_movement_ori_sd_SDLD)
PC2_movement_ori_sd_absSDLD <- abs(PC2_movement_ori_sd_SDLD)

PC1_movement_ori_sd_absSDSD <- abs(PC1_movement_ori_sd_SDSD)
PC2_movement_ori_sd_absSDSD <- abs(PC2_movement_ori_sd_SDSD)

PC1_movement_ori_sd_absLDLD <- abs(PC1_movement_ori_sd_LDLD)
PC2_movement_ori_sd_absLDLD <- abs(PC2_movement_ori_sd_LDLD)

#vel
PC1_movement_vel_sd_absSDLD <- abs(PC1_movement_vel_sd_SDLD)
PC2_movement_vel_sd_absSDLD <- abs(PC2_movement_vel_sd_SDLD)

PC1_movement_vel_sd_absSDSD <- abs(PC1_movement_vel_sd_SDSD)
PC2_movement_vel_sd_absSDSD <- abs(PC2_movement_vel_sd_SDSD)

PC1_movement_vel_sd_absLDLD <- abs(PC1_movement_vel_sd_LDLD)
PC2_movement_vel_sd_absLDLD <- abs(PC2_movement_vel_sd_LDLD)

#acc

PC1_movement_acc_sd_absSDLD <- abs(PC1_movement_acc_sd_SDLD)
PC2_movement_acc_sd_absSDLD <- abs(PC2_movement_acc_sd_SDLD)

PC1_movement_acc_sd_absSDSD <- abs(PC1_movement_acc_sd_SDSD)
PC2_movement_acc_sd_absSDSD <- abs(PC2_movement_acc_sd_SDSD)

PC1_movement_acc_sd_absLDLD <- abs(PC1_movement_acc_sd_LDLD)
PC2_movement_acc_sd_absLDLD <- abs(PC2_movement_acc_sd_LDLD)



#euclid distance between pre and post 

#sdld
distances_ori_sdSDLD <- sqrt(rowSums((fpca_ori$scores[1:64, 1:2] - fpca_ori$scores[65:128, 1:2])^2))
distances_vel_sdSDLD <- sqrt(rowSums((fpca_vel$scores[1:64, 1:2] - fpca_vel$scores[65:128, 1:2])^2))
distances_acc_sdSDLD <- sqrt(rowSums((fpca_acc$scores[1:64, 1:2] - fpca_acc$scores[65:128, 1:2])^2))

#sdsd
distances_ori_sdSDSD <- sqrt(rowSums((fpca_ori$scores[129:192, 1:2] - fpca_ori$scores[193:256, 1:2])^2))
distances_vel_sdSDSD <- sqrt(rowSums((fpca_vel$scores[129:192, 1:2] - fpca_vel$scores[193:256, 1:2])^2))
distances_acc_sdSDSD <- sqrt(rowSums((fpca_acc$scores[129:192, 1:2] - fpca_acc$scores[193:256, 1:2])^2))

#ldld
distances_ori_sdLDLD <- sqrt(rowSums((fpca_ori$scores[257:319, 1:2] - fpca_ori$scores[320:382, 1:2])^2))
distances_vel_sdLDLD <- sqrt(rowSums((fpca_vel$scores[257:319, 1:2] - fpca_vel$scores[320:382, 1:2])^2))
distances_acc_sdLDLD <- sqrt(rowSums((fpca_acc$scores[257:319, 1:2] - fpca_acc$scores[320:382, 1:2])^2))


#make a table of each of these.... 

table_of_phenotypes <- cbind(fpca_ori$scores[1:64,1],fpca_ori$scores[65:128,1],fpca_ori$scores[129:192,1],fpca_ori$scores[193:256,1],fpca_ori$scores[1:64,2],fpca_ori$scores[65:128,2],fpca_ori$scores[129:192,2],fpca_ori$scores[193:256,2],fpca_vel$scores[1:64,1],fpca_vel$scores[65:128,1],fpca_vel$scores[129:192,1],fpca_vel$scores[193:256,1],fpca_vel$scores[1:64,2],fpca_vel$scores[65:128,2],fpca_vel$scores[129:192,2],fpca_vel$scores[193:256,2],fpca_acc$scores[1:64,1],fpca_acc$scores[65:128,1],fpca_acc$scores[129:192,1],fpca_acc$scores[193:256,1],fpca_acc$scores[1:64,2],fpca_acc$scores[65:128,2],fpca_acc$scores[129:192,2],fpca_acc$scores[193:256,2],PC1_movement_ori_sd_SDLD, PC2_movement_ori_sd_SDLD, PC1_movement_ori_sd_SDSD,PC2_movement_ori_sd_SDSD ,PC1_movement_vel_sd_SDLD,PC2_movement_vel_sd_SDLD ,PC1_movement_vel_sd_SDSD ,PC2_movement_vel_sd_SDSD,PC1_movement_acc_sd_SDLD ,PC2_movement_acc_sd_SDLD ,PC1_movement_acc_sd_SDSD,PC2_movement_acc_sd_SDSD,PC1_movement_ori_sd_absSDLD,PC2_movement_ori_sd_absSDLD,PC1_movement_ori_sd_absSDSD,PC2_movement_ori_sd_absSDSD,PC1_movement_vel_sd_absSDLD,PC2_movement_vel_sd_absSDLD,PC1_movement_vel_sd_absSDSD,PC2_movement_vel_sd_absSDSD,PC1_movement_acc_sd_absSDLD,PC2_movement_acc_sd_absSDLD,PC1_movement_acc_sd_absSDSD,PC2_movement_acc_sd_absSDSD,distances_ori_sdSDLD,distances_vel_sdSDLD,distances_acc_sdSDLD,distances_ori_sdSDSD,distances_vel_sdSDSD,distances_acc_sdSDSD)

#
row.names(table_of_phenotypes) <- Genotype_names

colnames(table_of_phenotypes) <- c("varmx_sd_PC1_SDLDpre","varmx_sd_PC1_SDLDpost","varmx_sd_PC1_SDSDpre", "varmx_sd_PC1_SDSDpost","varmx_sd_PC2_SDLDpre","varmx_sd_PC2_SDLDpost","varmx_sd_PC2_SDSDpre", "varmx_sd_PC2_SDSDpost", "varmx_sd_PC1_SDLDpre_vel","varmx_sd_PC1_SDLDpost_vel","varmx_sd_PC1_SDSDpre_vel", "varmx_sd_PC1_SDSDpost_vel","varmx_sd_PC2_SDLDpre_vel","varmx_sd_PC2_SDLDpost_vel","varmx_sd_PC2_SDSDpre_vel", "varmx_sd_PC2_SDSDpost_vel","varmx_sd_PC1_SDLDpre_acc","varmx_sd_PC1_SDLDpost_acc","varmx_sd_PC1_SDSDpre_acc", "varmx_sd_PC1_SDSDpost_acc","varmx_sd_PC2_SDLDpre_acc","varmx_sd_PC2_SDLDpost_acc","varmx_sd_PC2_SDSDpre_acc", "varmx_sd_PC2_SDSDpost_acc","varmx_PC1_movement_ori_sd_SDLD","varmx_ PC2_movement_ori_sd_SDLD","varmx_ PC1_movement_ori_sd_SDSD","varmx_PC2_movement_ori_sd_SDSD ","varmx_PC1_movement_vel_sd_SDLD","varmx_PC2_movement_vel_sd_SDLD ","varmx_PC1_movement_vel_sd_SDSD ","varmx_PC2_movement_vel_sd_SDSD","varmx_PC1_movement_acc_sd_SDLD ","varmx_PC2_movement_acc_sd_SDLD ","varmx_PC1_movement_acc_sd_SDSD","varmx_PC2_movement_acc_sd_SDSD","varmx_PC1_movement_ori_sd_absSDLD","varmx_PC2_movement_ori_sd_absSDLD","varmx_PC1_movement_ori_sd_absSDSD","varmx_PC2_movement_ori_sd_absSDSD","varmx_PC1_movement_vel_sd_absSDLD","varmx_PC2_movement_vel_sd_absSDLD","varmx_PC1_movement_vel_sd_absSDSD","varmx_PC2_movement_vel_sd_absSDSD","varmx_PC1_movement_acc_sd_absSDLD","varmx_PC2_movement_acc_sd_absSDLD","varmx_PC1_movement_acc_sd_absSDSD","varmx_PC2_movement_acc_sd_absSDSD","varmx_distances_ori_sdSDLD","varmx_distances_vel_sdSDLD","varmx_distances_acc_sdSDLD","varmx_distances_ori_sdSDSD","varmx_distances_vel_sdSDSD","varmx_distances_acc_sdSDSD")

#because LDLD had a genotype missing (WxT2 do a separate table) can add in excel and just have it as

#add extra row to the LDLD sections 
table_of_phenotypesLDLD <- cbind(fpca_ori$scores[257:319,1],fpca_ori$scores[320:382,1],fpca_ori$scores[257:319,2],fpca_ori$scores[320:382,2],fpca_vel$scores[257:319,1],fpca_vel$scores[320:382,1],fpca_vel$scores[257:319,2],fpca_vel$scores[320:382,2],fpca_acc$scores[257:319,1],fpca_acc$scores[320:382,1],fpca_acc$scores[257:319,2],fpca_acc$scores[320:382,2],PC1_movement_acc_sd_LDLD ,PC2_movement_acc_sd_LDLD ,PC1_movement_ori_sd_LDLD,PC2_movement_ori_sd_LDLD,PC1_movement_vel_sd_LDLD ,PC2_movement_vel_sd_LDLD,PC1_movement_ori_sd_absLDLD,PC2_movement_ori_sd_absLDLD,PC1_movement_vel_sd_absLDLD,PC2_movement_vel_sd_absLDLD,PC1_movement_acc_sd_absLDLD,PC2_movement_acc_sd_absLDLD,distances_ori_sdLDLD,distances_vel_sdLDLD,distances_acc_sdLDLD)


colnames(table_of_phenotypesLDLD) <- c("varmx_sd_PC1_LDLDpre","varmx_sd_PC1_LDLDpost","varmx_sd_PC2_LDLDpre","varmx_sd_PC2_LDLDpost","varmx_sd_PC1_LDLDpre_vel","varmx_sd_PC1_LDLDpost_vel","varmx_sd_PC2_LDLDpre_vel","varmx_sd_PC2_LDLDpost_vel","varmx_sd_PC1_LDLDpre_acc","varmx_sd_PC1_LDLDpost_acc","varmx_sd_PC2_LDLDpre_acc","varmx_sd_PC2_LDLDpost_acc","varmx_PC1_movement_acc_sd_LDLD","varmx_PC2_movement_acc_sd_LDLD","varmx_PC1_movement_ori_sd_LDLD","varmx_PC2_movement_ori_sd_LDLD","varmx_PC1_movement_vel_sd_LDLD","varmx_PC2_movement_vel_sd_LDLD","varmx_PC1_movement_ori_sd_absLDLD","varmx_PC2_movement_ori_sd_absLDLD","varmx_PC1_movement_vel_sd_absLDLD","varmx_PC2_movement_vel_sd_absLDLD","varmx_PC1_movement_acc_sd_absLDLD","varmx_PC2_movement_acc_sd_absLDLD","varmx_distances_ori_sdLDLD","varmx_distances_vel_sdLDLD","varmx_distances_acc_sdLDLD")


sub_genotypes <- Genotype_names[-12]  
row.names(table_of_phenotypesLDLD) <- sub_genotypes


##
#save file
write.csv(table_of_phenotypes, "table_of_phenotypes_pre_post_sdcurves_varmx.csv")

write.csv(table_of_phenotypesLDLD, "table_of_phenotypes_pre_post_sdcurvesLDLD_varmx.csv")

