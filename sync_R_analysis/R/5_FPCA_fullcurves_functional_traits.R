# use smooth datatables to make functional object and perform FPCA across all experiments and produce a trait table.  ===============

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

#normalise before analysis 
Total_curves_ALL_rmout <- apply(Total_curves_ALL_rmout,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))

Total_curves_ALL_rmout <- data.matrix(Total_curves_ALL_rmout)

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 300)

#using high smothing to get overall shape of curves 
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
sp_totalsmooth <- smooth.basis(argvals = new_ZT, y = Total_curves_ALL_rmout, fdParobj = fdobj)
#plot(sp_totalsmooth, lty = "solid")

ori_obj <- deriv.fd(sp_totalsmooth$fd,0)
vel_obj <- deriv.fd(sp_totalsmooth$fd,1)
acc_obj <- deriv.fd(sp_totalsmooth$fd,2)


#FPCA

fpca_ori <- pca.fd(ori_obj,3,harmfdPar = fdPar(ori_obj),centerfns = FALSE)
#fpca_ori <- fda::varmx.pca.fd(fpca_ori)
scores <- fpca_ori$scores
rownames(scores) <- colnames(Total_curves_ALL_rmout)
#this plots each of the principle compments
plot(fpca_ori$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2", "FPC3"), col = c("black", "red", "green"), lty = 1)
fpca_ori$varprop*100
sum(fpca_ori$varprop*100)

###########################
fpca_vel <- pca.fd(vel_obj,3,harmfdPar = fdPar(vel_obj),centerfns = FALSE)
#fpca_vel <- fda::varmx.pca.fd(fpca_vel)

scores <- fpca_vel$scores

rownames(scores) <- colnames(Total_curves_ALL_rmout)

#this plots each of the principle compments
plot(fpca_vel$harmonics, lty = "solid", lwd =2 )
legend("bottomright", c("FPC1", "FPC2", "FPC3"), col = c("black", "red", "green"), lty = 1)
fpca_vel$varprop*100
sum(fpca_vel$varprop*100)


###########################
fpca_acc <- pca.fd(acc_obj,3,harmfdPar = fdPar(acc_obj),centerfns = FALSE)
#fpca_acc <- fda::varmx.pca.fd(fpca_acc)
scores <- fpca_acc$scores

rownames(scores) <- colnames(Total_curves_ALL_rmout)

#this plots each of the principle compments
plot(fpca_acc$harmonics, lty = "solid", lwd =2 )
legend("bottomright", c("FPC1", "FPC2", "FPC3"), col = c("black", "red", "green"), lty = 1)
fpca_acc$varprop*100
sum(fpca_acc$varprop*100)



##############                SCORES
#library(dplyr)
#Sum the x and y cooridants and average...

#make the my vetor twice as there is the SD and LD portion
my_vector_SDLD <- vector("numeric")
my_vector_SDSD <- vector("numeric")
my_vector_LDLD <- vector("numeric")


Total_curves_ALL_rmout <- data.frame(Total_curves_ALL_rmout)

Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names_sorted <- sort(Genotype_names)
library(dplyr)
for (i in Genotype_names_sorted) {
  my_vector_SDLD[i] <- length(Total_curves_SD_LD %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_SDSD[i] <- length(Total_curves_SD_SD %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_LDLD[i] <- length(Total_curves_LD_LD %>%
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

#SPREAD OF DATA workout the average euclinean distance to each center then divide by the number of reps in each genotype.

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

#look at the difference in the spread from the SD_LD and the SD_SD and LD_LD 

#origional
#if above 1 then control is more variable if less than 1 then SD_LD is more variable

PC1_diff_SDLD_SDSD_ori_spread <- spread_coords_mean_ori[65:128,1] / spread_coords_mean_ori[1:64,1]
PC2_diff_SDLD_SDSD_ori_spread <- spread_coords_mean_ori[65:128,2] / spread_coords_mean_ori[1:64,2]
PC3_diff_SDLD_SDSD_ori_spread <- spread_coords_mean_ori[65:128,3] / spread_coords_mean_ori[1:64,3]


PC1_diff_SDLD_LDLD_ori_spread <- spread_coords_mean_ori[129:192,1] / spread_coords_mean_ori[1:64,1]
PC2_diff_SDLD_LDLD_ori_spread <- spread_coords_mean_ori[129:192,2] / spread_coords_mean_ori[1:64,2]
PC3_diff_SDLD_LDLD_ori_spread <- spread_coords_mean_ori[129:192,3] / spread_coords_mean_ori[1:64,3]


#veelocity

PC1_diff_SDLD_SDSD_vel_spread <- spread_coords_mean_vel[65:128,1] / spread_coords_mean_vel[1:64,1]
PC2_diff_SDLD_SDSD_vel_spread <- spread_coords_mean_vel[65:128,2] / spread_coords_mean_vel[1:64,2]
PC3_diff_SDLD_SDSD_vel_spread <- spread_coords_mean_vel[65:128,3] / spread_coords_mean_vel[1:64,3]


PC1_diff_SDLD_LDLD_vel_spread <- spread_coords_mean_vel[129:192,1] / spread_coords_mean_vel[1:64,1]
PC2_diff_SDLD_LDLD_vel_spread <- spread_coords_mean_vel[129:192,2] / spread_coords_mean_vel[1:64,2]
PC3_diff_SDLD_LDLD_vel_spread <- spread_coords_mean_vel[129:192,3] / spread_coords_mean_vel[1:64,3]

#acceleration


PC1_diff_SDLD_SDSD_acc_spread <- spread_coords_mean_acc[65:128,1] / spread_coords_mean_acc[1:64,1]
PC2_diff_SDLD_SDSD_acc_spread <- spread_coords_mean_acc[65:128,2] / spread_coords_mean_acc[1:64,2]
PC3_diff_SDLD_SDSD_acc_spread <- spread_coords_mean_acc[65:128,3] / spread_coords_mean_acc[1:64,3]


PC1_diff_SDLD_LDLD_acc_spread <- spread_coords_mean_acc[129:192,1] / spread_coords_mean_acc[1:64,1]
PC2_diff_SDLD_LDLD_acc_spread <- spread_coords_mean_acc[129:192,2] / spread_coords_mean_acc[1:64,2]
PC3_diff_SDLD_LDLD_acc_spread <- spread_coords_mean_acc[129:192,3] / spread_coords_mean_acc[1:64,3]


#make the my vetor twice as there is the SD and LD portion
pre_my_vector_SDLD <- vector("numeric")
pre_my_vector_SDSD <- vector("numeric")
pre_my_vector_LDLD <- vector("numeric")


Total_curves_ALL <- data.frame(Total_curves_ALL_rmout)

Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names_sorted <- sort(Genotype_names)
library(dplyr)
for (i in Genotype_names_sorted) {
  pre_my_vector_SDLD[i] <- length(Total_curves_SD_LD %>%
                                    dplyr::select(dplyr::contains(i)))
  pre_my_vector_SDSD[i] <- length(Total_curves_SD_SD %>%
                                    dplyr::select(dplyr::contains(i)))
  pre_my_vector_LDLD[i] <- length(Total_curves_LD_LD %>%
                                    dplyr::select(dplyr::contains(i)))
}

pre_my_vector_SDLD <- pre_my_vector_SDLD[order(names(pre_my_vector_SDLD))]
pre_my_vector_SDSD <- pre_my_vector_SDSD[order(names(pre_my_vector_SDSD))]
pre_my_vector_LDLD <- pre_my_vector_LDLD[order(names(pre_my_vector_LDLD))]

#TRAIT of number of curves removed......... 

removed_outlier_only_SDLD <- pre_my_vector_SDLD - my_vector_SDLD
removed_outlier_only_SDSD <- pre_my_vector_SDSD - my_vector_SDSD
removed_outlier_only_LDLD <- pre_my_vector_LDLD - my_vector_LDLD

removed_fnt1_only_SDLD <-  48 - pre_my_vector_SDLD
removed_fnt1_only_SDSD <-  12 - pre_my_vector_SDSD
removed_fnt1_only_LDLD <-  12 - pre_my_vector_LDLD

#divide by the number of samples to begin with as they will be different for the controls ?
removed_total_SDLD <- 48 - my_vector_SDLD
removed_total_SDSD <- 12 - my_vector_SDSD
removed_total_LDLD <- 12 - my_vector_LDLD
#


################################################################

################################################################

#        TRAIT TABLE

################################################################

################################################################

#make a table with each of the results


table_of_phenotypes <- cbind(center_coords_ori[1:64,1],center_coords_ori[1:64,2],center_coords_ori[1:64,3],
                             
                             center_coords_ori[65:128,1],center_coords_ori[65:128,2],center_coords_ori[65:128,3],
                             center_coords_ori[129:192,1],center_coords_ori[129:192,2],center_coords_ori[129:192,3],
                             center_coords_vel[1:64,1],center_coords_vel[1:64,2],center_coords_vel[1:64,3],
                             center_coords_vel[65:128,1],center_coords_vel[65:128,2],center_coords_vel[65:128,3],
                             center_coords_vel[129:192,1],center_coords_vel[129:192,2],center_coords_vel[129:192,3],
                             center_coords_acc[1:64,1],center_coords_acc[1:64,2],center_coords_acc[1:64,3],
                             center_coords_acc[65:128,1],center_coords_acc[65:128,2],center_coords_acc[65:128,3],
                             center_coords_acc[129:192,1],center_coords_acc[129:192,2],center_coords_acc[129:192,3],
                             spread_coords_mean_ori[1:64,1],spread_coords_mean_ori[1:64,2],spread_coords_mean_ori[1:64,3],
                             spread_coords_mean_ori[65:128,1],spread_coords_mean_ori[65:128,2],spread_coords_mean_ori[65:128,3],
                             spread_coords_mean_ori[129:192,1],spread_coords_mean_ori[129:192,2],spread_coords_mean_ori[129:192,3],
                             spread_coords_mean_vel[1:64,1],spread_coords_mean_vel[1:64,2],spread_coords_mean_vel[1:64,3],
                             spread_coords_mean_vel[65:128,1],spread_coords_mean_vel[65:128,2],spread_coords_mean_vel[65:128,3],
                             spread_coords_mean_vel[129:192,1],spread_coords_mean_vel[129:192,2],spread_coords_mean_vel[129:192,3],
                             spread_coords_mean_acc[1:64,1],spread_coords_mean_acc[1:64,2],spread_coords_mean_acc[1:64,3],
                             spread_coords_mean_acc[65:128,1],spread_coords_mean_acc[65:128,2],spread_coords_mean_acc[65:128,3],
                             spread_coords_mean_acc[129:192,1],spread_coords_mean_acc[129:192,2],spread_coords_mean_acc[129:192,3],
                             spread_coords_mean_ori_total[1:64,1],
                             spread_coords_mean_ori_total[65:128,1],
                             spread_coords_mean_ori_total[129:192,1],
                             spread_coords_mean_vel_total[1:64,1],
                             spread_coords_mean_vel_total[65:128,1],
                             spread_coords_mean_vel_total[129:192,1],
                             spread_coords_mean_acc_total[1:64,1],
                             spread_coords_mean_acc_total[65:128,1],
                             spread_coords_mean_acc_total[129:192,1],
                             spread_coords_se_ori[1:64,1],spread_coords_se_ori[1:64,2],spread_coords_se_ori[1:64,3],
                             spread_coords_se_ori[65:128,1],spread_coords_se_ori[65:128,2],spread_coords_se_ori[65:128,3],
                             spread_coords_se_ori[129:192,1],spread_coords_se_ori[129:192,2],spread_coords_se_ori[129:192,3],
                             spread_coords_se_vel[1:64,1],spread_coords_se_vel[1:64,2],spread_coords_se_vel[1:64,3],
                             spread_coords_se_vel[65:128,1],spread_coords_se_vel[65:128,2],spread_coords_se_vel[65:128,3],
                             spread_coords_se_vel[129:192,1],spread_coords_se_vel[129:192,2],spread_coords_se_vel[129:192,3],
                             spread_coords_se_acc[1:64,1],spread_coords_se_acc[1:64,2],spread_coords_se_acc[1:64,3],
                             spread_coords_se_acc[65:128,1],spread_coords_se_acc[65:128,2],spread_coords_se_acc[65:128,3],
                             spread_coords_se_acc[129:192,1],spread_coords_se_acc[129:192,2],spread_coords_se_acc[129:192,3],
                             spread_coords_se_ori_total[1:64,1],
                             spread_coords_se_ori_total[65:128,1],
                             spread_coords_se_ori_total[129:192,1],
                             spread_coords_se_vel_total[1:64,1],
                             spread_coords_se_vel_total[65:128,1],
                             spread_coords_se_vel_total[129:192,1],
                             spread_coords_se_acc_total[1:64,1],
                             spread_coords_se_acc_total[65:128,1],
                             spread_coords_se_acc_total[129:192,1],
                             
                             
                             
                             PC1_diff_SDLD_SDSD_ori_spread,PC2_diff_SDLD_SDSD_ori_spread,PC3_diff_SDLD_SDSD_ori_spread,
                             PC1_diff_SDLD_LDLD_ori_spread,PC2_diff_SDLD_LDLD_ori_spread,PC3_diff_SDLD_LDLD_ori_spread,
                             PC1_diff_SDLD_SDSD_vel_spread,PC2_diff_SDLD_SDSD_vel_spread,PC3_diff_SDLD_SDSD_vel_spread,
                             PC1_diff_SDLD_LDLD_vel_spread,PC2_diff_SDLD_LDLD_vel_spread,PC3_diff_SDLD_LDLD_vel_spread,
                             PC1_diff_SDLD_SDSD_acc_spread,PC2_diff_SDLD_SDSD_acc_spread,PC3_diff_SDLD_SDSD_acc_spread,
                             PC1_diff_SDLD_LDLD_acc_spread,PC2_diff_SDLD_LDLD_acc_spread,PC3_diff_SDLD_LDLD_acc_spread,
                             
                             removed_outlier_only_SDLD,removed_outlier_only_SDSD,removed_outlier_only_LDLD,removed_fnt1_only_SDLD,removed_fnt1_only_SDSD ,removed_fnt1_only_LDLD,removed_total_SDLD ,removed_total_SDSD,removed_total_LDLD)


#
row.names(table_of_phenotypes) <- Genotype_names_sorted

names1 <- c("PC1_gc_SDLD","PC2_gc_SDLD","PC3_gc_SDLD","PC1_gc_SDSD","PC2_gc_SDSD","PC3_gc_SDSD","PC1_gc_LDLD","PC2_gc_LDLD","PC3_gc_LDLD","PC1_gc_SDLD_vel","PC2_gc_SDLD_vel","PC3_gc_SDLD_vel","PC1_gc_SDSD_vel","PC2_gc_SDSD_vel","PC3_gc_SDSD_vel","PC1_gc_LDLD_vel","PC2_gc_LDLD_vel","PC3_gc_LDLD_vel","PC1_gc_SDLD_acc","PC2_gc_SDLD_acc","PC3_gc_SDLD_acc","PC1_gc_SDSD_acc","PC2_gc_SDSD_acc","PC3_gc_SDSD_acc","PC1_gc_LDLD_acc","PC2_gc_LDLD_acc","PC3_gc_LDLD_acc","PC1_mean_dist_SDLD_ori","PC2_mean_dist_SDLD_ori","PC3_mean_dist_SDLD_ori","PC1_mean_dist_SDSD_ori","PC2_mean_dist_SDSD_ori","PC3_mean_dist_SDSD_ori","PC1_mean_dist_LDLD_ori","PC2_mean_dist_LDLD_ori","PC3_mean_dist_LDLD_ori","PC1_mean_dist_SDLD_vel","PC2_mean_dist_SDLD_vel","PC3_mean_dist_SDLD_vel","PC1_mean_dist_SDSD_vel","PC2_mean_dist_SDSD_vel","PC3_mean_dist_SDSD_vel","PC1_mean_dist_LDLD_vel","PC2_mean_dist_LDLD_vel","PC3_mean_dist_LDLD_vel","PC1_mean_dist_SDLD_acc","PC2_mean_dist_SDLD_acc","PC3_mean_dist_SDLD_acc","PC1_mean_dist_SDSD_acc","PC2_mean_dist_SDSD_acc","PC3_mean_dist_SDSD_acc","PC1_mean_dist_LDLD_acc","PC2_mean_dist_LDLD_acc","PC3_mean_dist_LDLD_acc","mean_total_euclid_SDLD_ori","mean_total_euclid_SDSD_ori","mean_total_euclid_LDLD_ori","mean_total_euclid_SDLD_vel","mean_total_euclid_SDSD_vel","mean_total_euclid_LDLD_vel","mean_total_euclid_SDLD_acc","mean_total_euclid_SDSD_acc","mean_total_euclid_LDLD_acc","PC1_se_dist_SDLD_ori","PC2_se_dist_SDLD_ori","PC3_se_dist_SDLD_ori","PC1_se_dist_SDSD_ori","PC2_se_dist_SDSD_ori","PC3_se_dist_SDSD_ori","PC1_se_dist_LDLD_ori","PC2_se_dist_LDLD_ori","PC3_se_dist_LDLD_ori","PC1_se_dist_SDLD_vel","PC2_se_dist_SDLD_vel","PC3_se_dist_SDLD_vel","PC1_se_dist_SDSD_vel","PC2_se_dist_SDSD_vel","PC3_se_dist_SDSD_vel","PC1_se_dist_LDLD_vel","PC2_se_dist_LDLD_vel","PC3_se_dist_LDLD_vel","PC1_se_dist_SDLD_acc","PC2_se_dist_SDLD_acc","PC3_se_dist_SDLD_acc","PC1_se_dist_SDSD_acc","PC2_se_dist_SDSD_acc","PC3_se_dist_SDSD_acc","PC1_se_dist_LDLD_acc","PC2_se_dist_LDLD_acc","PC3_se_dist_LDLD_acc","se_total_euclid_SDLD_ori","se_total_euclid_SDSD_ori","se_total_euclid_LDLD_ori","se_total_euclid_SDLD_vel","se_total_euclid_SDSD_vel","se_total_euclid_LDLD_vel","se_total_euclid_SDLD_acc","se_total_euclid_SDSD_acc","se_total_euclid_LDLD_acc","PC1_diff_SDLD_SDSD_ori_spread","PC2_diff_SDLD_SDSD_ori_spread","PC3_diff_SDLD_SDSD_ori_spread","PC1_diff_SDLD_LDLD_ori_spread","PC2_diff_SDLD_LDLD_ori_spread","PC3_diff_SDLD_LDLD_ori_spread","PC1_diff_SDLD_SDSD_vel_spread","PC2_diff_SDLD_SDSD_vel_spread","PC3_diff_SDLD_SDSD_vel_spread","PC1_diff_SDLD_LDLD_vel_spread","PC2_diff_SDLD_LDLD_vel_spread","PC3_diff_SDLD_LDLD_vel_spread","PC1_diff_SDLD_SDSD_acc_spread","PC2_diff_SDLD_SDSD_acc_spread","PC3_diff_SDLD_SDSD_acc_spread","PC1_diff_SDLD_LDLD_acc_spread","PC2_diff_SDLD_LDLD_acc_spread","PC3_diff_SDLD_LDLD_acc_spread","removed_outlier_only_SDLD","removed_outlier_only_SDSD","removed_outlier_only_LDLD","removed_fnt1_only_SDLD","removed_fnt1_only_SDSD" ,"removed_fnt1_only_LDLD","removed_total_SDLD" ,"removed_total_SDSD","removed_total_LDLD")
colnames(table_of_phenotypes) <- names1


################################################################
#     add on traits just looking at FPCA of SDLD 


#on SDLD only 

Total_curves_SD_LD_rmout <- apply(Total_curves_SD_LD,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))
Total_curves_SD_LD_rmout <- data.matrix(Total_curves_SD_LD_rmout)
BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 300)
#using high smothing to get overall shape of curves 
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
sp_totalsmooth_SDLD <- smooth.basis(argvals = new_ZT, y = Total_curves_SD_LD_rmout, fdParobj = fdobj)
#plot(sp_totalsmooth, lty = "solid")
ori_objSDLD <- deriv.fd(sp_totalsmooth_SDLD$fd,0)
vel_objSDLD <- deriv.fd(sp_totalsmooth_SDLD$fd,1)
acc_objSDLD <- deriv.fd(sp_totalsmooth_SDLD$fd,2)

#FPCA
fpca_oriSDLD <- pca.fd(ori_objSDLD,3,harmfdPar = fdPar(ori_objSDLD),centerfns = FALSE)

scoresSDLD <- fpca_oriSDLD$scores
rownames(scoresSDLD) <- colnames(Total_curves_SD_LD_rmout)

###########################
fpca_velSDLD <- pca.fd(vel_objSDLD,3,harmfdPar = fdPar(vel_objSDLD),centerfns = FALSE)

###########################
fpca_accSDLD <- pca.fd(acc_objSDLD,3,harmfdPar = fdPar(acc_objSDLD),centerfns = FALSE)

##############                SCORES
#library(dplyr)
#Sum the x and y cooridants and average...

#make the my vetor twice as there is the SD and LD portion
my_vector_SDLD <- vector("numeric")

Total_curves_SD_LD_rmout <- data.frame(Total_curves_SD_LD_rmout)

Genotype_names_sorted <- sort(Genotype_names)
library(dplyr)
for (i in Genotype_names_sorted) {
  my_vector_SDLD[i] <- length(Total_curves_SD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
}
my_vector_SDLD <- my_vector_SDLD[order(names(my_vector_SDLD))]

my_vector <- my_vector_SDLD
cum <- cumsum(my_vector)
genotype_list <- list()

center_coords_oriSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_oriSDLD[i,] <- c(mean(fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}
center_coords_velSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_velSDLD[i,] <- c(mean(fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}

center_coords_accSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_accSDLD[i,] <- c(mean(fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}



#SPREAD OF DATA workout the average euclinean distance to each center then divide by the number of reps in each genotype.

spread_coords_mean_oriSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_oriSDLD[i,] <- c((sum(abs(center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_velSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_velSDLD[i,] <- c((sum(abs(center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_accSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_accSDLD[i,] <- c((sum(abs(center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_mean_oriSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_oriSDLD_total[i,1] <- sum(sqrt(((center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_velSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_velSDLD_total[i,1] <- sum(sqrt(((center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_accSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_accSDLD_total[i,1] <- sum(sqrt(((center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}




####################### now the standard error of the spread 


spread_coords_se_oriSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_oriSDLD[i,] <- c((sd(abs(center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_velSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_velSDLD[i,] <- c((sd(abs(center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_accSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_accSDLD[i,] <- c((sd(abs(center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_se_oriSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_oriSDLD_total[i,1] <- sd(sqrt(((center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_velSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_velSDLD_total[i,1] <- sd(sqrt(((center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_accSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_accSDLD_total[i,1] <- sd(sqrt(((center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

################################################################



#make a table with each of the results


table_of_phenotypesSDLD <- cbind(center_coords_oriSDLD[1:64,1],center_coords_oriSDLD[1:64,2],center_coords_oriSDLD[1:64,3],
                                 center_coords_velSDLD[1:64,1],center_coords_velSDLD[1:64,2],center_coords_velSDLD[1:64,3],
                                 center_coords_accSDLD[1:64,1],center_coords_accSDLD[1:64,2],center_coords_accSDLD[1:64,3],
                                 spread_coords_mean_oriSDLD[1:64,1],spread_coords_mean_oriSDLD[1:64,2],spread_coords_mean_oriSDLD[1:64,3],
                                 spread_coords_mean_velSDLD[1:64,1],spread_coords_mean_velSDLD[1:64,2],spread_coords_mean_velSDLD[1:64,3],
                                 spread_coords_mean_accSDLD[1:64,1],spread_coords_mean_accSDLD[1:64,2],spread_coords_mean_accSDLD[1:64,3],
                                 spread_coords_mean_oriSDLD_total[1:64,1],spread_coords_mean_velSDLD_total[1:64,1],spread_coords_mean_accSDLD_total[1:64,1],spread_coords_se_oriSDLD[1:64,1],spread_coords_se_oriSDLD[1:64,2],spread_coords_se_oriSDLD[1:64,3],
                                 spread_coords_se_velSDLD[1:64,1],spread_coords_se_velSDLD[1:64,2],spread_coords_se_velSDLD[1:64,3],
                                 spread_coords_se_accSDLD[1:64,1],spread_coords_se_accSDLD[1:64,2],spread_coords_se_accSDLD[1:64,3],
                                 spread_coords_se_oriSDLD_total[1:64,1],spread_coords_se_velSDLD_total[1:64,1],spread_coords_se_accSDLD_total[1:64,1])


#
row.names(table_of_phenotypesSDLD) <- Genotype_names_sorted

colnames(table_of_phenotypesSDLD) <- c("SDLD_PC1_center","SDLD_PC2_center","SDLD_PC3_center", "SDLD_PC1_center_vel","SDLD_PC2_center_vel","SDLD_PC3_center_vel", "SDLD_PC1_center_acc","SDLD_PC2_center_acc","SDLD_PC3_center_acc", "SDLD_PC1_mean_dsit","SDLD_PC2_mean_dsit","SDLD_PC3_mean_dsit","SDLD_PC1_mean_dsit_vel","SDLD_PC2_mean_dsit_vel","SDLD_PC3_mean_dsit_vel","SDLD_PC1_mean_dsit_acc","SDLD_PC2_mean_dsit_acc","SDLD_PC3_mean_dsit_acc", "SDLD_total_euclid_mean_dist_ori","SDLD_total_euclid_mean_dist_vel","SDLD_total_euclid_mean_dist_acc","SDLD_PC1_sd_dsit","SDLD_PC2_sd_dsit","SDLD_PC3_sd_dsit","SDLD_PC1_sd_dsit_vel","SDLD_PC2_sd_dsit_vel","SDLD_PC3_sd_dsit_vel","SDLD_PC1_sd_dsit_acc","SDLD_PC2_sd_dsit_acc","SDLD_PC3_sd_dsit_acc", "SDLD_total_euclid_sd_dist_ori","SDLD_total_euclid_sd_dist_vel","SDLD_total_euclid_sd_dist_acc")





#combine the two 

table_of_phenotypes_all <- cbind(table_of_phenotypes,table_of_phenotypesSDLD)

##
#save file 
write.csv(table_of_phenotypes, "table_of_phenotypes_fullcurves.csv")


#
##############################################

#same ut on varmx rotation of fpca 



#normalise before analysis 
Total_curves_ALL_rmout <- apply(Total_curves_ALL_rmout,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))

Total_curves_ALL_rmout <- data.matrix(Total_curves_ALL_rmout)

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 300)

#using high smothing to get overall shape of curves 
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
sp_totalsmooth <- smooth.basis(argvals = new_ZT, y = Total_curves_ALL_rmout, fdParobj = fdobj)
#plot(sp_totalsmooth, lty = "solid")

ori_obj <- deriv.fd(sp_totalsmooth$fd,0)
vel_obj <- deriv.fd(sp_totalsmooth$fd,1)
acc_obj <- deriv.fd(sp_totalsmooth$fd,2)


#FPCA

fpca_ori <- pca.fd(ori_obj,3,harmfdPar = fdPar(ori_obj),centerfns = FALSE)
fpca_ori <- fda::varmx.pca.fd(fpca_ori)
scores <- fpca_ori$scores
rownames(scores) <- colnames(Total_curves_ALL_rmout)

###########################
fpca_vel <- pca.fd(vel_obj,3,harmfdPar = fdPar(vel_obj),centerfns = FALSE)
fpca_vel <- fda::varmx.pca.fd(fpca_vel)
scores <- fpca_vel$scores
rownames(scores) <- colnames(Total_curves_ALL_rmout)

###########################
fpca_acc <- pca.fd(acc_obj,3,harmfdPar = fdPar(acc_obj),centerfns = FALSE)
fpca_acc <- fda::varmx.pca.fd(fpca_acc)
scores <- fpca_acc$scores
rownames(scores) <- colnames(Total_curves_ALL_rmout)

##############                SCORES
#library(dplyr)
#Sum the x and y cooridants and average...

#make the my vetor twice as there is the SD and LD portion
my_vector_SDLD <- vector("numeric")
my_vector_SDSD <- vector("numeric")
my_vector_LDLD <- vector("numeric")


Total_curves_ALL_rmout <- data.frame(Total_curves_ALL_rmout)

Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names_sorted <- sort(Genotype_names)
library(dplyr)
for (i in Genotype_names_sorted) {
  my_vector_SDLD[i] <- length(Total_curves_SD_LD %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_SDSD[i] <- length(Total_curves_SD_SD %>%
                                dplyr::select(dplyr::contains(i)))
  my_vector_LDLD[i] <- length(Total_curves_LD_LD %>%
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

#SPREAD OF DATA workout the average euclinean distance to each center then divide by the number of reps in each genotype.

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



#look at the difference in the spread from the SD_LD and the SD_SD and LD_LD 

#origional
#if above 1 then control is more variable if less than 1 then SD_LD is more variable

PC1_diff_SDLD_SDSD_ori_spread <- spread_coords_mean_ori[65:128,1] / spread_coords_mean_ori[1:64,1]
PC2_diff_SDLD_SDSD_ori_spread <- spread_coords_mean_ori[65:128,2] / spread_coords_mean_ori[1:64,2]
PC3_diff_SDLD_SDSD_ori_spread <- spread_coords_mean_ori[65:128,3] / spread_coords_mean_ori[1:64,3]


PC1_diff_SDLD_LDLD_ori_spread <- spread_coords_mean_ori[129:192,1] / spread_coords_mean_ori[1:64,1]
PC2_diff_SDLD_LDLD_ori_spread <- spread_coords_mean_ori[129:192,2] / spread_coords_mean_ori[1:64,2]
PC3_diff_SDLD_LDLD_ori_spread <- spread_coords_mean_ori[129:192,3] / spread_coords_mean_ori[1:64,3]


#veelocity

PC1_diff_SDLD_SDSD_vel_spread <- spread_coords_mean_vel[65:128,1] / spread_coords_mean_vel[1:64,1]
PC2_diff_SDLD_SDSD_vel_spread <- spread_coords_mean_vel[65:128,2] / spread_coords_mean_vel[1:64,2]
PC3_diff_SDLD_SDSD_vel_spread <- spread_coords_mean_vel[65:128,3] / spread_coords_mean_vel[1:64,3]


PC1_diff_SDLD_LDLD_vel_spread <- spread_coords_mean_vel[129:192,1] / spread_coords_mean_vel[1:64,1]
PC2_diff_SDLD_LDLD_vel_spread <- spread_coords_mean_vel[129:192,2] / spread_coords_mean_vel[1:64,2]
PC3_diff_SDLD_LDLD_vel_spread <- spread_coords_mean_vel[129:192,3] / spread_coords_mean_vel[1:64,3]

#acceleration


PC1_diff_SDLD_SDSD_acc_spread <- spread_coords_mean_acc[65:128,1] / spread_coords_mean_acc[1:64,1]
PC2_diff_SDLD_SDSD_acc_spread <- spread_coords_mean_acc[65:128,2] / spread_coords_mean_acc[1:64,2]
PC3_diff_SDLD_SDSD_acc_spread <- spread_coords_mean_acc[65:128,3] / spread_coords_mean_acc[1:64,3]


PC1_diff_SDLD_LDLD_acc_spread <- spread_coords_mean_acc[129:192,1] / spread_coords_mean_acc[1:64,1]
PC2_diff_SDLD_LDLD_acc_spread <- spread_coords_mean_acc[129:192,2] / spread_coords_mean_acc[1:64,2]
PC3_diff_SDLD_LDLD_acc_spread <- spread_coords_mean_acc[129:192,3] / spread_coords_mean_acc[1:64,3]


################################################################

#        TRAIT TABLE

################################################################

################################################################

#make a table with each of the results


table_of_phenotypes_varmx <- cbind(center_coords_ori[1:64,1],center_coords_ori[1:64,2],center_coords_ori[1:64,3],
                                   center_coords_ori[65:128,1],center_coords_ori[65:128,2],center_coords_ori[65:128,3],
                                   center_coords_ori[129:192,1],center_coords_ori[129:192,2],center_coords_ori[129:192,3],
                                   center_coords_vel[1:64,1],center_coords_vel[1:64,2],center_coords_vel[1:64,3],
                                   center_coords_vel[65:128,1],center_coords_vel[65:128,2],center_coords_vel[65:128,3],
                                   center_coords_vel[129:192,1],center_coords_vel[129:192,2],center_coords_vel[129:192,3],
                                   center_coords_acc[1:64,1],center_coords_acc[1:64,2],center_coords_acc[1:64,3],
                                   center_coords_acc[65:128,1],center_coords_acc[65:128,2],center_coords_acc[65:128,3],
                                   center_coords_acc[129:192,1],center_coords_acc[129:192,2],center_coords_acc[129:192,3],
                                   spread_coords_mean_ori[1:64,1],spread_coords_mean_ori[1:64,2],spread_coords_mean_ori[1:64,3],
                                   spread_coords_mean_ori[65:128,1],spread_coords_mean_ori[65:128,2],spread_coords_mean_ori[65:128,3],
                                   spread_coords_mean_ori[129:192,1],spread_coords_mean_ori[129:192,2],spread_coords_mean_ori[129:192,3],
                                   spread_coords_mean_vel[1:64,1],spread_coords_mean_vel[1:64,2],spread_coords_mean_vel[1:64,3],
                                   spread_coords_mean_vel[65:128,1],spread_coords_mean_vel[65:128,2],spread_coords_mean_vel[65:128,3],
                                   spread_coords_mean_vel[129:192,1],spread_coords_mean_vel[129:192,2],spread_coords_mean_vel[129:192,3],
                                   spread_coords_mean_acc[1:64,1],spread_coords_mean_acc[1:64,2],spread_coords_mean_acc[1:64,3],
                                   spread_coords_mean_acc[65:128,1],spread_coords_mean_acc[65:128,2],spread_coords_mean_acc[65:128,3],
                                   spread_coords_mean_acc[129:192,1],spread_coords_mean_acc[129:192,2],spread_coords_mean_acc[129:192,3],
                                   spread_coords_mean_ori_total[1:64,1],
                                   spread_coords_mean_ori_total[65:128,1],
                                   spread_coords_mean_ori_total[129:192,1],
                                   spread_coords_mean_vel_total[1:64,1],
                                   spread_coords_mean_vel_total[65:128,1],
                                   spread_coords_mean_vel_total[129:192,1],
                                   spread_coords_mean_acc_total[1:64,1],
                                   spread_coords_mean_acc_total[65:128,1],
                                   spread_coords_mean_acc_total[129:192,1],
                                   spread_coords_se_ori[1:64,1],spread_coords_se_ori[1:64,2],spread_coords_se_ori[1:64,3],
                                   spread_coords_se_ori[65:128,1],spread_coords_se_ori[65:128,2],spread_coords_se_ori[65:128,3],
                                   spread_coords_se_ori[129:192,1],spread_coords_se_ori[129:192,2],spread_coords_se_ori[129:192,3],
                                   spread_coords_se_vel[1:64,1],spread_coords_se_vel[1:64,2],spread_coords_se_vel[1:64,3],
                                   spread_coords_se_vel[65:128,1],spread_coords_se_vel[65:128,2],spread_coords_se_vel[65:128,3],
                                   spread_coords_se_vel[129:192,1],spread_coords_se_vel[129:192,2],spread_coords_se_vel[129:192,3],
                                   spread_coords_se_acc[1:64,1],spread_coords_se_acc[1:64,2],spread_coords_se_acc[1:64,3],
                                   spread_coords_se_acc[65:128,1],spread_coords_se_acc[65:128,2],spread_coords_se_acc[65:128,3],
                                   spread_coords_se_acc[129:192,1],spread_coords_se_acc[129:192,2],spread_coords_se_acc[129:192,3],
                                   spread_coords_se_ori_total[1:64,1],
                                   spread_coords_se_ori_total[65:128,1],
                                   spread_coords_se_ori_total[129:192,1],
                                   spread_coords_se_vel_total[1:64,1],
                                   spread_coords_se_vel_total[65:128,1],
                                   spread_coords_se_vel_total[129:192,1],
                                   spread_coords_se_acc_total[1:64,1],
                                   spread_coords_se_acc_total[65:128,1],
                                   spread_coords_se_acc_total[129:192,1],
                                   
                                   
                                   
                                   PC1_diff_SDLD_SDSD_ori_spread,PC2_diff_SDLD_SDSD_ori_spread,PC3_diff_SDLD_SDSD_ori_spread,
                                   PC1_diff_SDLD_LDLD_ori_spread,PC2_diff_SDLD_LDLD_ori_spread,PC3_diff_SDLD_LDLD_ori_spread,
                                   PC1_diff_SDLD_SDSD_vel_spread,PC2_diff_SDLD_SDSD_vel_spread,PC3_diff_SDLD_SDSD_vel_spread,
                                   PC1_diff_SDLD_LDLD_vel_spread,PC2_diff_SDLD_LDLD_vel_spread,PC3_diff_SDLD_LDLD_vel_spread,
                                   PC1_diff_SDLD_SDSD_acc_spread,PC2_diff_SDLD_SDSD_acc_spread,PC3_diff_SDLD_SDSD_acc_spread,
                                   PC1_diff_SDLD_LDLD_acc_spread,PC2_diff_SDLD_LDLD_acc_spread,PC3_diff_SDLD_LDLD_acc_spread
)


#
row.names(table_of_phenotypes_varmx) <- Genotype_names_sorted

names1_varmx <- c("varmx_PC1_gc_SDLD","varmx_PC2_gc_SDLD","varmx_PC3_gc_SDLD","varmx_PC1_gc_SDSD","varmx_PC2_gc_SDSD","varmx_PC3_gc_SDSD","varmx_PC1_gc_LDLD","varmx_PC2_gc_LDLD","varmx_PC3_gc_LDLD","varmx_PC1_gc_SDLD_vel","varmx_PC2_gc_SDLD_vel","varmx_PC3_gc_SDLD_vel","varmx_PC1_gc_SDSD_vel","varmx_PC2_gc_SDSD_vel","varmx_PC3_gc_SDSD_vel","varmx_PC1_gc_LDLD_vel","varmx_PC2_gc_LDLD_vel","varmx_PC3_gc_LDLD_vel","varmx_PC1_gc_SDLD_acc","varmx_PC2_gc_SDLD_acc","varmx_PC3_gc_SDLD_acc","varmx_PC1_gc_SDSD_acc","varmx_PC2_gc_SDSD_acc","varmx_PC3_gc_SDSD_acc","varmx_PC1_gc_LDLD_acc","varmx_PC2_gc_LDLD_acc","varmx_PC3_gc_LDLD_acc","varmx_PC1_mean_dist_SDLD_ori","varmx_PC2_mean_dist_SDLD_ori","varmx_PC3_mean_dist_SDLD_ori","varmx_PC1_mean_dist_SDSD_ori","varmx_PC2_mean_dist_SDSD_ori","varmx_PC3_mean_dist_SDSD_ori","varmx_PC1_mean_dist_LDLD_ori","varmx_PC2_mean_dist_LDLD_ori","varmx_PC3_mean_dist_LDLD_ori","varmx_PC1_mean_dist_SDLD_vel","varmx_PC2_mean_dist_SDLD_vel","varmx_PC3_mean_dist_SDLD_vel","varmx_PC1_mean_dist_SDSD_vel","varmx_PC2_mean_dist_SDSD_vel","varmx_PC3_mean_dist_SDSD_vel","varmx_PC1_mean_dist_LDLD_vel","varmx_PC2_mean_dist_LDLD_vel","varmx_PC3_mean_dist_LDLD_vel","varmx_PC1_mean_dist_SDLD_acc","varmx_PC2_mean_dist_SDLD_acc","varmx_PC3_mean_dist_SDLD_acc","varmx_PC1_mean_dist_SDSD_acc","varmx_PC2_mean_dist_SDSD_acc","varmx_PC3_mean_dist_SDSD_acc","varmx_PC1_mean_dist_LDLD_acc","varmx_PC2_mean_dist_LDLD_acc","varmx_PC3_mean_dist_LDLD_acc","varmx_mean_total_euclid_SDLD_ori","varmx_mean_total_euclid_SDSD_ori","varmx_mean_total_euclid_LDLD_ori","varmx_mean_total_euclid_SDLD_vel","varmx_mean_total_euclid_SDSD_vel","varmx_mean_total_euclid_LDLD_vel","varmx_mean_total_euclid_SDLD_acc","varmx_mean_total_euclid_SDSD_acc","varmx_mean_total_euclid_LDLD_acc","varmx_PC1_se_dist_SDLD_ori","varmx_PC2_se_dist_SDLD_ori","varmx_PC3_se_dist_SDLD_ori","varmx_PC1_se_dist_SDSD_ori","varmx_PC2_se_dist_SDSD_ori","varmx_PC3_se_dist_SDSD_ori","varmx_PC1_se_dist_LDLD_ori","varmx_PC2_se_dist_LDLD_ori","varmx_PC3_se_dist_LDLD_ori","varmx_PC1_se_dist_SDLD_vel","varmx_PC2_se_dist_SDLD_vel","varmx_PC3_se_dist_SDLD_vel","varmx_PC1_se_dist_SDSD_vel","varmx_PC2_se_dist_SDSD_vel","varmx_PC3_se_dist_SDSD_vel","varmx_PC1_se_dist_LDLD_vel","varmx_PC2_se_dist_LDLD_vel","varmx_PC3_se_dist_LDLD_vel","varmx_PC1_se_dist_SDLD_acc","varmx_PC2_se_dist_SDLD_acc","varmx_PC3_se_dist_SDLD_acc","varmx_PC1_se_dist_SDSD_acc","varmx_PC2_se_dist_SDSD_acc","varmx_PC3_se_dist_SDSD_acc","varmx_PC1_se_dist_LDLD_acc","varmx_PC2_se_dist_LDLD_acc","varmx_PC3_se_dist_LDLD_acc","varmx_se_total_euclid_SDLD_ori","varmx_se_total_euclid_SDSD_ori","varmx_se_total_euclid_LDLD_ori","varmx_se_total_euclid_SDLD_vel","varmx_se_total_euclid_SDSD_vel","varmx_se_total_euclid_LDLD_vel","varmx_se_total_euclid_SDLD_acc","varmx_se_total_euclid_SDSD_acc","varmx_se_total_euclid_LDLD_acc","varmx_PC1_diff_SDLD_SDSD_ori_spread","varmx_PC2_diff_SDLD_SDSD_ori_spread","varmx_PC3_diff_SDLD_SDSD_ori_spread","varmx_PC1_diff_SDLD_LDLD_ori_spread","varmx_PC2_diff_SDLD_LDLD_ori_spread","varmx_PC3_diff_SDLD_LDLD_ori_spread","varmx_PC1_diff_SDLD_SDSD_vel_spread","varmx_PC2_diff_SDLD_SDSD_vel_spread","varmx_PC3_diff_SDLD_SDSD_vel_spread","varmx_PC1_diff_SDLD_LDLD_vel_spread","varmx_PC2_diff_SDLD_LDLD_vel_spread","varmx_PC3_diff_SDLD_LDLD_vel_spread","varmx_PC1_diff_SDLD_SDSD_acc_spread","varmx_PC2_diff_SDLD_SDSD_acc_spread","varmx_PC3_diff_SDLD_SDSD_acc_spread","varmx_PC1_diff_SDLD_LDLD_acc_spread","varmx_PC2_diff_SDLD_LDLD_acc_spread","varmx_PC3_diff_SDLD_LDLD_acc_spread")
colnames(table_of_phenotypes_varmx) <- names1_varmx


################################################################
#     add on traits just looking at FPCA of SDLD 


#on SDLD only 

Total_curves_SD_LD_rmout <- apply(Total_curves_SD_LD,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))
Total_curves_SD_LD_rmout <- data.matrix(Total_curves_SD_LD_rmout)
BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 300)
#using high smothing to get overall shape of curves 
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
sp_totalsmooth_SDLD <- smooth.basis(argvals = new_ZT, y = Total_curves_SD_LD_rmout, fdParobj = fdobj)
#plot(sp_totalsmooth, lty = "solid")
ori_objSDLD <- deriv.fd(sp_totalsmooth_SDLD$fd,0)
vel_objSDLD <- deriv.fd(sp_totalsmooth_SDLD$fd,1)
acc_objSDLD <- deriv.fd(sp_totalsmooth_SDLD$fd,2)

#FPCA
fpca_oriSDLD <- pca.fd(ori_objSDLD,3,harmfdPar = fdPar(ori_objSDLD),centerfns = FALSE)
fpca_oriSDLD <- varmx.pca.fd(fpca_oriSDLD)
scoresSDLD <- fpca_oriSDLD$scores
rownames(scoresSDLD) <- colnames(Total_curves_SD_LD_rmout)

###########################
fpca_velSDLD <- pca.fd(vel_objSDLD,3,harmfdPar = fdPar(vel_objSDLD),centerfns = FALSE)
fpca_velSDLD <- varmx.pca.fd(fpca_velSDLD)
###########################
fpca_accSDLD <- pca.fd(acc_objSDLD,3,harmfdPar = fdPar(acc_objSDLD),centerfns = FALSE)
fpca_accSDLD <- varmx.pca.fd(fpca_accSDLD)
##############                SCORES
#library(dplyr)
#Sum the x and y cooridants and average...

#make the my vetor twice as there is the SD and LD portion
my_vector_SDLD <- vector("numeric")

Total_curves_SD_LD_rmout <- data.frame(Total_curves_SD_LD_rmout)

Genotype_names_sorted <- sort(Genotype_names)
library(dplyr)
for (i in Genotype_names_sorted) {
  my_vector_SDLD[i] <- length(Total_curves_SD_LD_rmout %>%
                                dplyr::select(dplyr::contains(i)))
}
my_vector_SDLD <- my_vector_SDLD[order(names(my_vector_SDLD))]

my_vector <- my_vector_SDLD
cum <- cumsum(my_vector)
genotype_list <- list()

center_coords_oriSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_oriSDLD[i,] <- c(mean(fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}
center_coords_velSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_velSDLD[i,] <- c(mean(fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}

center_coords_accSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_accSDLD[i,] <- c(mean(fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_acc$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}


#SPREAD OF DATA workout the average euclinean distance to each center then divide by the number of reps in each genotype.

spread_coords_mean_oriSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_oriSDLD[i,] <- c((sum(abs(center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_velSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_velSDLD[i,] <- c((sum(abs(center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

spread_coords_mean_accSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_accSDLD[i,] <- c((sum(abs(center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/my_vector[i]),(sum(abs(center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/my_vector[i]),(sum(abs(center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/my_vector[i])) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_mean_oriSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_oriSDLD_total[i,1] <- sum(sqrt(((center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_velSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_velSDLD_total[i,1] <- sum(sqrt(((center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}

spread_coords_mean_accSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_mean_accSDLD_total[i,1] <- sum(sqrt(((center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/ my_vector[i]
}




####################### now the standard error of the spread 


spread_coords_se_oriSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_oriSDLD[i,] <- c((sd(abs(center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_velSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_velSDLD[i,] <- c((sd(abs(center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

spread_coords_se_accSDLD <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_accSDLD[i,] <- c((sd(abs(center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]))/sqrt(my_vector[i])),(sd(abs(center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]))/sqrt(my_vector[i])),(sd(abs(center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))/sqrt(my_vector[i]))) 
}

#####################################

# do the total euclidean distance betwen all three of the pcs 


spread_coords_se_oriSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_oriSDLD_total[i,1] <- sd(sqrt(((center_coords_oriSDLD[i,1] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_oriSDLD[i,2] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_oriSDLD[i,3] - fpca_oriSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_velSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_velSDLD_total[i,1] <- sd(sqrt(((center_coords_velSDLD[i,1] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_velSDLD[i,2] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_velSDLD[i,3] - fpca_velSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}

spread_coords_se_accSDLD_total <-  matrix(ncol = 1, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  spread_coords_se_accSDLD_total[i,1] <- sd(sqrt(((center_coords_accSDLD[i,1] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],1])^2)+ ((center_coords_accSDLD[i,2] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],2])^2)+ (center_coords_accSDLD[i,3] - fpca_accSDLD$scores[(1 + (cum[i] - my_vector[i])):cum[i],3])^2))/sqrt(my_vector[i])
}


################################################################

#make a table with each of the results


table_of_phenotypesSDLD_varmx <- cbind(center_coords_oriSDLD[1:64,1],center_coords_oriSDLD[1:64,2],center_coords_oriSDLD[1:64,3],
                                       center_coords_velSDLD[1:64,1],center_coords_velSDLD[1:64,2],center_coords_velSDLD[1:64,3],
                                       center_coords_accSDLD[1:64,1],center_coords_accSDLD[1:64,2],center_coords_accSDLD[1:64,3],
                                       spread_coords_mean_oriSDLD[1:64,1],spread_coords_mean_oriSDLD[1:64,2],spread_coords_mean_oriSDLD[1:64,3],
                                       spread_coords_mean_velSDLD[1:64,1],spread_coords_mean_velSDLD[1:64,2],spread_coords_mean_velSDLD[1:64,3],
                                       spread_coords_mean_accSDLD[1:64,1],spread_coords_mean_accSDLD[1:64,2],spread_coords_mean_accSDLD[1:64,3],
                                       spread_coords_mean_oriSDLD_total[1:64,1],spread_coords_mean_velSDLD_total[1:64,1],spread_coords_mean_accSDLD_total[1:64,1],spread_coords_se_oriSDLD[1:64,1],spread_coords_se_oriSDLD[1:64,2],spread_coords_se_oriSDLD[1:64,3],
                                       spread_coords_se_velSDLD[1:64,1],spread_coords_se_velSDLD[1:64,2],spread_coords_se_velSDLD[1:64,3],
                                       spread_coords_se_accSDLD[1:64,1],spread_coords_se_accSDLD[1:64,2],spread_coords_se_accSDLD[1:64,3],
                                       spread_coords_se_oriSDLD_total[1:64,1],spread_coords_se_velSDLD_total[1:64,1],spread_coords_se_accSDLD_total[1:64,1])


#
row.names(table_of_phenotypesSDLD_varmx) <- Genotype_names_sorted

colnames(table_of_phenotypesSDLD_varmx) <- c("varmx_SDLD_PC1_center","varmx_SDLD_PC2_center","varmx_SDLD_PC3_center","varmx_SDLD_PC1_center_vel","varmx_SDLD_PC2_center_vel","varmx_SDLD_PC3_center_vel","varmx_SDLD_PC1_center_acc","varmx_SDLD_PC2_center_acc","varmx_SDLD_PC3_center_acc","varmx_SDLD_PC1_mean_dsit","varmx_SDLD_PC2_mean_dsit","varmx_SDLD_PC3_mean_dsit","varmx_SDLD_PC1_mean_dsit_vel","varmx_SDLD_PC2_mean_dsit_vel","varmx_SDLD_PC3_mean_dsit_vel","varmx_SDLD_PC1_mean_dsit_acc","varmx_SDLD_PC2_mean_dsit_acc","varmx_SDLD_PC3_mean_dsit_acc","varmx_SDLD_total_euclid_mean_dist_ori","varmx_SDLD_total_euclid_mean_dist_vel","varmx_SDLD_total_euclid_mean_dist_acc","varmx_SDLD_PC1_sd_dsit","varmx_SDLD_PC2_sd_dsit","varmx_SDLD_PC3_sd_dsit","varmx_SDLD_PC1_sd_dsit_vel","varmx_SDLD_PC2_sd_dsit_vel","varmx_SDLD_PC3_sd_dsit_vel","varmx_SDLD_PC1_sd_dsit_acc","varmx_SDLD_PC2_sd_dsit_acc","varmx_SDLD_PC3_sd_dsit_acc","varmx_SDLD_total_euclid_sd_dist_ori","varmx_SDLD_total_euclid_sd_dist_vel","varmx_SDLD_total_euclid_sd_dist_acc")

#combine the two 

table_of_phenotypes_all_varmx <- cbind(table_of_phenotypes_varmx,table_of_phenotypesSDLD_varmx)

##
#save file 
write.csv(table_of_phenotypes_all_varmx, "table_of_phenotypes_fullcurves_varmx.csv")
