# estimating functions based of the curves of pre and post shift in all conditions , doing FPCA and identifiying functional traits. 

#read in files
curves_pre_post_SDLD <- read.csv(file = "curves_pre_post_SDLD.csv")
curves_pre_post_SDSD  <- read.csv(file = "curves_pre_post_SDSD.csv")
curves_pre_post_LDLD  <- read.csv(file = "curves_pre_post_LDLD.csv")

# all curves were evalutaed under the same time frame (48 hours) so xan just combine the dataframes togetther
new_ZT <- seq(1,48,length.out = 100)

curves_pre_post_SDLD <- curves_pre_post_SDLD[,-1]
curves_pre_post_SDSD <- curves_pre_post_SDSD[,-1]
curves_pre_post_LDLD <- curves_pre_post_LDLD[,-1]

curves_SDLD_pre <- curves_pre_post_SDLD[,1:2121]
curves_SDLD_post <- curves_pre_post_SDLD[,2122:4242]

curves_SDSD_pre <- curves_pre_post_SDSD[,1:587]
curves_SDSD_post <- curves_pre_post_SDSD[,588:1174]

curves_LDLD_pre <- curves_pre_post_LDLD[,1:493]
curves_LDLD_post <- curves_pre_post_LDLD[,494:986]


Total_curves_ALL_rmout <- cbind(curves_SDLD_pre,curves_SDSD_pre,curves_LDLD_pre,curves_SDLD_post,curves_SDSD_post,curves_LDLD_post)


Total_curves_ALL_rmout <- apply(Total_curves_ALL_rmout,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))

Total_curves_ALL_rmout <- data.matrix(Total_curves_ALL_rmout)

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 100)

fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
smooth_both <- smooth.basis(argvals = new_ZT, y = Total_curves_ALL_rmout, fdParobj = fdobj)


ori_obj <- deriv.fd(smooth_both$fd,0)
vel_obj <- deriv.fd(smooth_both$fd,1)
acc_obj <- deriv.fd(smooth_both$fd,2)


#FPCA

fpca_ori <- pca.fd(ori_obj,3,harmfdPar = fdPar(ori_obj),centerfns = FALSE)
scores <- fpca_ori$scores
rownames(scores) <- colnames(Total_curves_ALL_rmout)
#this plots each of the principle compments
plot(fpca_ori$harmonics, lty = "solid", lwd =2, main = "SD_LD, SD_SD, LD_LD origional curves" )
legend("bottomright", c("FPC1", "FPC2"), col = c("black", "red"), lty = 1)
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
Total_curves_SD_LD_rmout <- data.frame(curves_SDLD_pre)
Total_curves_SD_SD_rmout <- data.frame(curves_SDSD_pre)
Total_curves_LD_LD_rmout <- data.frame(curves_LDLD_pre)



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

#combine twice because its the pre section and then the post section

my_vector <- c(my_vector_SDLD,my_vector_SDSD,my_vector_LDLD,my_vector_SDLD,my_vector_SDSD,my_vector_LDLD)

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



#    VELPCOTY






#1:64 is SDLD pre
#65:128 is SDSD pre
#129:192 is LDLD pre
#193:256 is SDLD post
#257:320 is SDSD post
#321:384 is LDLD post

PC1_movement_ori_SDLD <- center_coords_ori[193:256,1] - center_coords_ori[1:64,1]  
PC2_movement_ori_SDLD <- center_coords_ori[193:256,2] - center_coords_ori[1:64,2] 
#vel
PC1_movement_vel_SDLD <- center_coords_vel[193:256,1] - center_coords_vel[1:64,1]  
PC2_movement_vel_SDLD <- center_coords_vel[193:256,2] - center_coords_vel[1:64,2] 
PC3_movement_vel_SDLD <- center_coords_vel[193:256,3] - center_coords_vel[1:64,3] 
#acc
PC1_movement_acc_SDLD <- center_coords_acc[193:256,1] - center_coords_acc[1:64,1]  
PC2_movement_acc_SDLD <- center_coords_acc[193:256,2] - center_coords_acc[1:64,2] 
PC3_movement_acc_SDLD <- center_coords_acc[193:256,3] - center_coords_acc[1:64,3] 

#SDSD
PC1_movement_ori_SDSD <- center_coords_ori[257:320,1] - center_coords_ori[65:128,1]  
PC2_movement_ori_SDSD <- center_coords_ori[257:320,2] - center_coords_ori[65:128,2] 
#vel
PC1_movement_vel_SDSD <- center_coords_vel[257:320,1] - center_coords_vel[65:128,1]  
PC2_movement_vel_SDSD <- center_coords_vel[257:320,2] - center_coords_vel[65:128,2] 
PC3_movement_vel_SDSD <- center_coords_vel[257:320,3] - center_coords_vel[65:128,3] 
#acc
PC1_movement_acc_SDSD <- center_coords_acc[257:320,1] - center_coords_acc[65:128,1]  
PC2_movement_acc_SDSD <- center_coords_acc[257:320,2] - center_coords_acc[65:128,2] 
PC3_movement_acc_SDSD <- center_coords_acc[257:320,3] - center_coords_acc[65:128,3] 


#LDLD
PC1_movement_ori_LDLD <- center_coords_ori[321:384,1] - center_coords_ori[129:192,1]  
PC2_movement_ori_LDLD <- center_coords_ori[321:384,2] - center_coords_ori[129:192,2] 
#vel
PC1_movement_vel_LDLD <- center_coords_vel[321:384,1] - center_coords_vel[129:192,1]  
PC2_movement_vel_LDLD <- center_coords_vel[321:384,2] - center_coords_vel[129:192,2] 
PC3_movement_vel_LDLD <- center_coords_vel[321:384,3] - center_coords_vel[129:192,3] 
#acc
PC1_movement_acc_LDLD <- center_coords_acc[321:384,1] - center_coords_acc[129:192,1]  
PC2_movement_acc_LDLD <- center_coords_acc[321:384,2] - center_coords_acc[129:192,2] 
PC3_movement_acc_LDLD <- center_coords_acc[321:384,3] - center_coords_acc[129:192,3] 

#first eculidean distance on total

#ori curves
euclid_dist_pre_post_ori_t <- matrix(ncol = 1, nrow = (length(my_vector)/2))
for (i in 1:(length(my_vector)/2)) {
  
  euclid_dist_pre_post_ori_t[i,] <- sqrt(((center_coords_ori[i,1] - center_coords_ori[i +192 ,1])^2)+((center_coords_ori[i,2] - center_coords_ori[i+192,2])^2)+((center_coords_ori[i,3] - center_coords_ori[i+192,3])^2))
  
}

#vel curves
euclid_dist_pre_post_vel_t <- matrix(ncol = 1, nrow = (length(my_vector)/2))
for (i in 1:(length(my_vector)/2)) {
  
  euclid_dist_pre_post_vel_t[i,] <- sqrt(((center_coords_vel[i,1] - center_coords_vel[i +192 ,1])^2)+((center_coords_vel[i,2] - center_coords_vel[i+192,2])^2)+((center_coords_vel[i,3] - center_coords_vel[i+192,3])^2))
  
}

#acc curves
euclid_dist_pre_post_acc_t <- matrix(ncol = 1, nrow = (length(my_vector)/2))
for (i in 1:(length(my_vector)/2)) {
  
  euclid_dist_pre_post_acc_t[i,] <- sqrt(((center_coords_acc[i,1] - center_coords_acc[i +192 ,1])^2)+((center_coords_acc[i,2] - center_coords_acc[i+192,2])^2)+((center_coords_acc[i,3] - center_coords_acc[i+192,3])^2))
  
}


######## now look at the the change separately
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


#####Eulcidiena distance total of going to the genotype centers
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

#################################################



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



#need to compare the before and after for the mean distance and spread and the se of distance and spread.... 

#do before take away the after 


change_in_mean_dist_ori <- spread_coords_mean_ori[1:192,1:3] - spread_coords_mean_ori[193:384,1:3]
change_in_mean_dist_vel <- spread_coords_mean_vel[1:192,1:3] - spread_coords_mean_vel[193:384,1:3]
change_in_mean_dist_acc <- spread_coords_mean_acc[1:192,1:3] - spread_coords_mean_acc[193:384,1:3]

change_in_mean_total_euclid_ori <- spread_coords_mean_ori_total[1:192,1] - spread_coords_mean_ori_total[193:384,1]
change_in_mean_total_euclid_vel <- spread_coords_mean_vel_total[1:192,1] - spread_coords_mean_vel_total[193:384,1]
change_in_mean_total_euclid_acc <- spread_coords_mean_acc_total[1:192,1] - spread_coords_mean_acc_total[193:384,1]

change_in_se_dist_ori <- spread_coords_se_ori[1:192,1:3] - spread_coords_se_ori[193:384,1:3]
change_in_se_dist_vel <- spread_coords_se_vel[1:192,1:3] - spread_coords_se_vel[193:384,1:3]
change_in_se_dist_acc <- spread_coords_se_acc[1:192,1:3] - spread_coords_se_acc[193:384,1:3]

change_in_se_total_euclid_ori <- spread_coords_se_ori_total[1:192,1] - spread_coords_se_ori_total[193:384,1]
change_in_se_total_euclid_vel <- spread_coords_se_vel_total[1:192,1] - spread_coords_se_vel_total[193:384,1]
change_in_se_total_euclid_acc <- spread_coords_se_acc_total[1:192,1] - spread_coords_se_acc_total[193:384,1]




#
#look at the difference in the spread from the pre and post shift (mostly in SD_LD)
#origional
#if posotive means after shift is LESS spread out


diff_SDLD_ori_centerdist <- spread_coords_mean_ori[1:64,1] - spread_coords_mean_ori[193:256,1]

diff_SDSD_ori_centerdist <- spread_coords_mean_ori[65:128,1] - spread_coords_mean_ori[257:320,1]

diff_LDLD_ori_centerdist <- spread_coords_mean_ori[129:192,1] - spread_coords_mean_ori[321:384,1]

ratio_SDLD_ori_centerdist <- spread_coords_mean_ori[1:64,1] / spread_coords_mean_ori[193:256,1]

ratio_SDSD_ori_centerdist <- spread_coords_mean_ori[65:128,1] / spread_coords_mean_ori[257:320,1]

ratio_LDLD_ori_centerdist <- spread_coords_mean_ori[129:192,1] / spread_coords_mean_ori[321:384,1]

#velocity

diff_SDLD_vel_centerdist <- spread_coords_mean_vel[1:64,1] - spread_coords_mean_vel[193:256,1]

diff_SDSD_vel_centerdist <- spread_coords_mean_vel[65:128,1] - spread_coords_mean_vel[257:320,1]

diff_LDLD_vel_centerdist <- spread_coords_mean_vel[129:192,1] - spread_coords_mean_vel[321:384,1]

#ratio of spread so do pre/post 

ratio_SDLD_vel_centerdist <- spread_coords_mean_vel[1:64,1] / spread_coords_mean_vel[193:256,1]

ratio_SDSD_vel_centerdist <- spread_coords_mean_vel[65:128,1] / spread_coords_mean_vel[257:320,1]

ratio_LDLD_vel_centerdist <- spread_coords_mean_vel[129:192,1] / spread_coords_mean_vel[321:384,1]

#acceleration

diff_SDLD_acc_centerdist <- spread_coords_mean_acc[1:64,1] - spread_coords_mean_acc[193:256,1]

diff_SDSD_acc_centerdist <- spread_coords_mean_acc[65:128,1] - spread_coords_mean_acc[257:320,1]

diff_LDLD_acc_centerdist <- spread_coords_mean_acc[129:192,1] - spread_coords_mean_acc[321:384,1]

#ratio of spread so do pre/post 

ratio_SDLD_acc_centerdist <- spread_coords_mean_acc[1:64,1] / spread_coords_mean_acc[193:256,1]

ratio_SDSD_acc_centerdist <- spread_coords_mean_acc[65:128,1] / spread_coords_mean_acc[257:320,1]

ratio_LDLD_acc_centerdist <- spread_coords_mean_acc[129:192,1] / spread_coords_mean_acc[321:384,1]
################################################################

#        TRAIT TABLE

################################################################

#make a table with each of the results

table_of_phenotypes <- cbind(center_coords_ori[1:64,1],center_coords_ori[1:64,2],center_coords_ori[1:64,3],
                             center_coords_ori[65:128,1],center_coords_ori[65:128,2],center_coords_ori[65:128,3],
                             center_coords_ori[129:192,1],center_coords_ori[129:192,2],center_coords_ori[129:192,3],
                             center_coords_ori[193:256,1],center_coords_ori[193:256,2],center_coords_ori[193:256,3],
                             center_coords_ori[257:320,1],center_coords_ori[257:320,2],center_coords_ori[257:320,3],
                             center_coords_ori[321:384,1],center_coords_ori[321:384,2],center_coords_ori[321:384,3],
                             center_coords_vel[1:64,1],center_coords_vel[1:64,2],center_coords_vel[1:64,3],
                             center_coords_vel[65:128,1],center_coords_vel[65:128,2],center_coords_vel[65:128,3],
                             center_coords_vel[129:192,1],center_coords_vel[129:192,2],center_coords_vel[129:192,3],
                             center_coords_vel[193:256,1],center_coords_vel[193:256,2],center_coords_vel[193:256,3],
                             center_coords_vel[257:320,1],center_coords_vel[257:320,2],center_coords_vel[257:320,3],
                             center_coords_vel[321:384,1],center_coords_vel[321:384,2],center_coords_vel[321:384,3],
                             center_coords_acc[1:64,1],center_coords_acc[1:64,2],center_coords_acc[1:64,3],
                             center_coords_acc[65:128,1],center_coords_acc[65:128,2],center_coords_acc[65:128,3],
                             center_coords_acc[129:192,1],center_coords_acc[129:192,2],center_coords_acc[129:192,3],
                             center_coords_acc[193:256,1],center_coords_acc[193:256,2],center_coords_acc[193:256,3],
                             center_coords_acc[257:320,1],center_coords_acc[257:320,2],center_coords_acc[257:320,3],
                             center_coords_acc[321:384,1],center_coords_acc[321:384,2],center_coords_acc[321:384,3],euclid_dist_pre_post_ori_t[1:64],euclid_dist_pre_post_ori_t[65:128],euclid_dist_pre_post_ori_t[129:192],euclid_dist_pre_post_vel_t[1:64],euclid_dist_pre_post_vel_t[65:128],euclid_dist_pre_post_vel_t[129:192],euclid_dist_pre_post_acc_t[1:64],euclid_dist_pre_post_acc_t[65:128],euclid_dist_pre_post_acc_t[129:192],
                             spread_coords_mean_ori[1:64,1],spread_coords_mean_ori[1:64,2],spread_coords_mean_ori[1:64,3],
                             spread_coords_mean_ori[65:128,1],spread_coords_mean_ori[65:128,2],spread_coords_mean_ori[65:128,3],
                             spread_coords_mean_ori[129:192,1],spread_coords_mean_ori[129:192,2],spread_coords_mean_ori[129:192,3],
                             spread_coords_mean_ori[193:256,1],spread_coords_mean_ori[193:256,2],spread_coords_mean_ori[193:256,3],
                             spread_coords_mean_ori[257:320,1],spread_coords_mean_ori[257:320,2],spread_coords_mean_ori[257:320,3],
                             spread_coords_mean_ori[321:384,1],spread_coords_mean_ori[321:384,2],spread_coords_mean_ori[321:384,3],
                             
                             spread_coords_mean_vel[1:64,1],spread_coords_mean_vel[1:64,2],spread_coords_mean_vel[1:64,3],
                             spread_coords_mean_vel[65:128,1],spread_coords_mean_vel[65:128,2],spread_coords_mean_vel[65:128,3],
                             spread_coords_mean_vel[129:192,1],spread_coords_mean_vel[129:192,2],spread_coords_mean_vel[129:192,3],
                             spread_coords_mean_vel[193:256,1],spread_coords_mean_vel[193:256,2],spread_coords_mean_vel[193:256,3],
                             spread_coords_mean_vel[257:320,1],spread_coords_mean_vel[257:320,2],spread_coords_mean_vel[257:320,3],
                             spread_coords_mean_vel[321:384,1],spread_coords_mean_vel[321:384,2],spread_coords_mean_vel[321:384,3],
                             
                             spread_coords_mean_acc[1:64,1],spread_coords_mean_acc[1:64,2],spread_coords_mean_acc[1:64,3],
                             spread_coords_mean_acc[65:128,1],spread_coords_mean_acc[65:128,2],spread_coords_mean_acc[65:128,3],
                             spread_coords_mean_acc[129:192,1],spread_coords_mean_acc[129:192,2],spread_coords_mean_acc[129:192,3],
                             spread_coords_mean_acc[193:256,1],spread_coords_mean_acc[193:256,2],spread_coords_mean_acc[193:256,3],
                             spread_coords_mean_acc[257:320,1],spread_coords_mean_acc[257:320,2],spread_coords_mean_acc[257:320,3],
                             spread_coords_mean_acc[321:384,1],spread_coords_mean_acc[321:384,2],spread_coords_mean_acc[321:384,3],
                             spread_coords_mean_ori_total[1:64,1],
                             spread_coords_mean_ori_total[65:128,1],
                             spread_coords_mean_ori_total[129:192,1],
                             spread_coords_mean_ori_total[193:256,1],
                             spread_coords_mean_ori_total[257:320,1],
                             spread_coords_mean_ori_total[321:384,1],
                             spread_coords_mean_vel_total[1:64,1],
                             spread_coords_mean_vel_total[65:128,1],
                             spread_coords_mean_vel_total[129:192,1],
                             spread_coords_mean_vel_total[193:256,1],
                             spread_coords_mean_vel_total[257:320,1],
                             spread_coords_mean_vel_total[321:384,1],
                             spread_coords_mean_acc_total[1:64,1],
                             spread_coords_mean_acc_total[65:128,1],
                             spread_coords_mean_acc_total[129:192,1],
                             spread_coords_mean_acc_total[193:256,1],
                             spread_coords_mean_acc_total[257:320,1],
                             spread_coords_mean_acc_total[321:384,1],
                             
                             
                             spread_coords_se_ori[1:64,1],spread_coords_se_ori[1:64,2],spread_coords_se_ori[1:64,3],
                             spread_coords_se_ori[65:128,1],spread_coords_se_ori[65:128,2],spread_coords_se_ori[65:128,3],
                             spread_coords_se_ori[129:192,1],spread_coords_se_ori[129:192,2],spread_coords_se_ori[129:192,3],
                             spread_coords_se_ori[193:256,1],spread_coords_se_ori[193:256,2],spread_coords_se_ori[193:256,3],
                             spread_coords_se_ori[257:320,1],spread_coords_se_ori[257:320,2],spread_coords_se_ori[257:320,3],
                             spread_coords_se_ori[321:384,1],spread_coords_se_ori[321:384,2],spread_coords_se_ori[321:384,3],
                             
                             spread_coords_se_vel[1:64,1],spread_coords_se_vel[1:64,2],spread_coords_se_vel[1:64,3],
                             spread_coords_se_vel[65:128,1],spread_coords_se_vel[65:128,2],spread_coords_se_vel[65:128,3],
                             spread_coords_se_vel[129:192,1],spread_coords_se_vel[129:192,2],spread_coords_se_vel[129:192,3],
                             spread_coords_se_vel[193:256,1],spread_coords_se_vel[193:256,2],spread_coords_se_vel[193:256,3],
                             spread_coords_se_vel[257:320,1],spread_coords_se_vel[257:320,2],spread_coords_se_vel[257:320,3],
                             spread_coords_se_vel[321:384,1],spread_coords_se_vel[321:384,2],spread_coords_se_vel[321:384,3],
                             
                             spread_coords_se_acc[1:64,1],spread_coords_se_acc[1:64,2],spread_coords_se_acc[1:64,3],
                             spread_coords_se_acc[65:128,1],spread_coords_se_acc[65:128,2],spread_coords_se_acc[65:128,3],
                             spread_coords_se_acc[129:192,1],spread_coords_se_acc[129:192,2],spread_coords_se_acc[129:192,3],
                             spread_coords_se_acc[193:256,1],spread_coords_se_acc[193:256,2],spread_coords_se_acc[193:256,3],
                             spread_coords_se_acc[257:320,1],spread_coords_se_acc[257:320,2],spread_coords_se_acc[257:320,3],
                             spread_coords_se_acc[321:384,1],spread_coords_se_acc[321:384,2],spread_coords_se_acc[321:384,3],
                             spread_coords_se_ori_total[1:64,1],
                             spread_coords_se_ori_total[65:128,1],
                             spread_coords_se_ori_total[129:192,1],
                             spread_coords_se_ori_total[193:256,1],
                             spread_coords_se_ori_total[257:320,1],
                             spread_coords_se_ori_total[321:384,1],
                             spread_coords_se_vel_total[1:64,1],
                             spread_coords_se_vel_total[65:128,1],
                             spread_coords_se_vel_total[129:192,1],
                             spread_coords_se_vel_total[193:256,1],
                             spread_coords_se_vel_total[257:320,1],
                             spread_coords_se_vel_total[321:384,1],
                             spread_coords_se_acc_total[1:64,1],
                             spread_coords_se_acc_total[65:128,1],
                             spread_coords_se_acc_total[129:192,1],
                             spread_coords_se_acc_total[193:256,1],
                             spread_coords_se_acc_total[257:320,1],
                             spread_coords_se_acc_total[321:384,1],
                             
                             PC1_movement_ori_SDLD,PC2_movement_ori_SDLD ,PC1_movement_vel_SDLD ,PC2_movement_vel_SDLD ,PC3_movement_vel_SDLD ,PC1_movement_acc_SDLD   
                             ,PC2_movement_acc_SDLD,PC3_movement_acc_SDLD ,PC1_movement_ori_SDSD,PC2_movement_ori_SDSD ,PC1_movement_vel_SDSD ,PC2_movement_vel_SDSD ,PC3_movement_vel_SDSD ,PC1_movement_acc_SDSD   ,PC2_movement_acc_SDSD ,PC3_movement_acc_SDSD,PC1_movement_ori_LDLD  ,PC2_movement_ori_LDLD  ,PC1_movement_vel_LDLD ,PC2_movement_vel_LDLD ,PC3_movement_vel_LDLD ,PC1_movement_acc_LDLD  ,PC2_movement_acc_LDLD  ,PC3_movement_acc_LDLD,                
                             diff_SDLD_ori_centerdist,diff_SDSD_ori_centerdist,diff_LDLD_ori_centerdist,diff_SDLD_vel_centerdist,diff_SDSD_vel_centerdist,diff_LDLD_vel_centerdist,diff_SDLD_acc_centerdist,diff_SDSD_acc_centerdist,diff_LDLD_acc_centerdist,ratio_SDLD_ori_centerdist,ratio_SDSD_ori_centerdist,ratio_LDLD_ori_centerdist,ratio_SDLD_vel_centerdist,ratio_SDSD_vel_centerdist,ratio_LDLD_vel_centerdist,ratio_SDLD_acc_centerdist,ratio_SDSD_acc_centerdist,ratio_LDLD_acc_centerdist,change_in_mean_dist_ori[1:64,1:3],change_in_mean_dist_ori[65:128,1:3],change_in_mean_dist_ori[129:192,1:3],
                             change_in_mean_dist_vel[1:64,1:3],change_in_mean_dist_vel[65:128,1:3],change_in_mean_dist_vel[129:192,1:3],
                             change_in_mean_dist_acc[1:64,1:3],change_in_mean_dist_acc[65:128,1:3],change_in_mean_dist_acc[129:192,1:3],
                             change_in_mean_total_euclid_ori[1:64],change_in_mean_total_euclid_ori[65:128],change_in_mean_total_euclid_ori[129:192],
                             change_in_se_dist_ori[1:64,1:3],change_in_se_dist_ori[65:128,1:3],change_in_se_dist_ori[129:192,1:3],
                             change_in_se_dist_vel[1:64,1:3],change_in_se_dist_vel[65:128,1:3],change_in_se_dist_vel[129:192,1:3],
                             change_in_se_dist_acc[1:64,1:3],change_in_se_dist_acc[65:128,1:3],change_in_se_dist_acc[129:192,1:3],
                             change_in_se_total_euclid_ori[1:64],change_in_se_total_euclid_ori[65:128],change_in_se_total_euclid_ori[129:192])

#
row.names(table_of_phenotypes) <- sort(Genotype_names)

colnames(table_of_phenotypes) <- c("PC1_gc_SDLD_pre","PC2_gc_SDLD_pre","PC3_gc_SDLD_pre",
                                   "PC1_gc_SDSD_pre","PC2_gc_SDSD_pre","PC3_gc_SDSD_pre",
                                   "PC1_gc_LDLD_pre","PC2_gc_LDLD_pre","PC3_gc_LDLD_pre",
                                   "PC1_gc_SDLD_post","PC2_gc_SDLD_post","PC3_gc_SDLD_post",
                                   "PC1_gc_SDSD_post","PC2_gc_SDSD_post","PC3_gc_SDSD_post",
                                   "PC1_gc_LDLD_post","PC2_gc_LDLD_post","PC3_gc_LDLD_post",
                                   "PC1_gc_vel_SDLD_pre","PC2_gc_vel_SDLD_pre","PC3_gc_vel_SDLD_pre",
                                   "PC1_gc_vel_SDSD_pre","PC2_gc_vel_SDSD_pre","PC3_gc_vel_SDSD_pre",
                                   "PC1_gc_vel_LDLD_pre","PC2_gc_vel_LDLD_pre","PC3_gc_vel_LDLD_pre",
                                   "PC1_gc_vel_SDLD_post","PC2_gc_vel_SDLD_post","PC3_gc_vel_SDLD_post",
                                   "PC1_gc_vel_SDSD_post","PC2_gc_vel_SDSD_post","PC3_gc_vel_SDSD_post",
                                   "PC1_gc_vel_LDLD_post","PC2_gc_vel_LDLD_post","PC3_gc_vel_LDLD_post",
                                   "PC1_gc_acc_SDLD_pre","PC2_gc_acc_SDLD_pre","PC3_gc_acc_SDLD_pre",
                                   "PC1_gc_acc_SDSD_pre","PC2_gc_acc_SDSD_pre","PC3_gc_acc_SDSD_pre",
                                   "PC1_gc_acc_LDLD_pre","PC2_gc_acc_LDLD_pre","PC3_gc_acc_LDLD_pre",
                                   "PC1_gc_acc_SDLD_post","PC2_gc_acc_SDLD_post","PC3_gc_acc_SDLD_post",
                                   "PC1_gc_acc_SDSD_post","PC2_gc_acc_SDSD_post","PC3_gc_acc_SDSD_post",
                                   "PC1_gc_acc_LDLD_post","PC2_gc_acc_LDLD_post","PC3_gc_acc_LDLD_post",
                                   "euclid_dist_pre_post_ori_tSDLD","euclid_dist_pre_post_ori_tSDSD","euclid_dist_pre_post_ori_tLDLD","euclid_dist_pre_post_vel_tSDLD","euclid_dist_pre_post_vel_tSDSD","euclid_dist_pre_post_vel_tLDLD","euclid_dist_pre_post_acc_tSDLD","euclid_dist_pre_post_acc_tSDSD","euclid_dist_pre_post_acc_tLDLD",
                                   
                                   "PC1_mean_dist_SDLD_ori_pre","PC2_mean_dist_SDLD_ori_pre","PC3_mean_dist_SDLD_ori_pre",
                                   "PC1_mean_dist_SDSD_ori_pre","PC2_mean_dist_SDSD_ori_pre","PC3_mean_dist_SDSD_ori_pre",
                                   "PC1_mean_dist_LDLD_ori_pre","PC2_mean_dist_LDLD_ori_pre","PC3_mean_dist_LDLD_ori_pre",
                                   "PC1_mean_dist_SDLD_ori_post","PC2_mean_dist_SDLD_ori_post","PC3_mean_dist_SDLD_ori_post",
                                   "PC1_mean_dist_SDSD_ori_post","PC2_mean_dist_SDSD_ori_post","PC3_mean_dist_SDSD_ori_post",
                                   "PC1_mean_dist_LDLD_ori_post","PC2_mean_dist_LDLD_ori_post","PC3_mean_dist_LDLD_ori_post",
                                   
                                   "PC1_mean_dist_SDLD_vel_pre","PC2_mean_dist_SDLD_vel_pre","PC3_mean_dist_SDLD_vel_pre",
                                   "PC1_mean_dist_SDSD_vel_pre","PC2_mean_dist_SDSD_vel_pre","PC3_mean_dist_SDSD_vel_pre",
                                   "PC1_mean_dist_LDLD_vel_pre","PC2_mean_dist_LDLD_vel_pre","PC3_mean_dist_LDLD_vel_pre",
                                   "PC1_mean_dist_SDLD_vel_post","PC2_mean_dist_SDLD_vel_post","PC3_mean_dist_SDLD_vel_post",
                                   "PC1_mean_dist_SDSD_vel_post","PC2_mean_dist_SDSD_vel_post","PC3_mean_dist_SDSD_vel_post",
                                   "PC1_mean_dist_LDLD_vel_post","PC2_mean_dist_LDLD_vel_post","PC3_mean_dist_LDLD_vel_post",
                                   
                                   "PC1_mean_dist_SDLD_acc_pre","PC2_mean_dist_SDLD_acc_pre","PC3_mean_dist_SDLD_acc_pre",
                                   "PC1_mean_dist_SDSD_acc_pre","PC2_mean_dist_SDSD_acc_pre","PC3_mean_dist_SDSD_acc_pre",
                                   "PC1_mean_dist_LDLD_acc_pre","PC2_mean_dist_LDLD_acc_pre","PC3_mean_dist_LDLD_acc_pre",
                                   "PC1_mean_dist_SDLD_acc_post","PC2_mean_dist_SDLD_acc_post","PC3_mean_dist_SDLD_acc_post",
                                   "PC1_mean_dist_SDSD_acc_post","PC2_mean_dist_SDSD_acc_post","PC3_mean_dist_SDSD_acc_post",
                                   "PC1_mean_dist_LDLD_acc_post","PC2_mean_dist_LDLD_acc_post","PC3_mean_dist_LDLD_acc_post",
                                   
                                   "mean_total_euclid_SDLD_ori_pre","mean_total_euclid_SDSD_ori_pre","mean_total_euclid_LDLD_ori_pre",
                                   "mean_total_euclid_SDLD_ori_post","mean_total_euclid_SDSD_ori_post","mean_total_euclid_LDLD_ori_post",
                                   
                                   "mean_total_euclid_SDLD_vel_pre","mean_total_euclid_SDSD_vel_pre","mean_total_euclid_LDLD_vel_pre",
                                   "mean_total_euclid_SDLD_vel_post","mean_total_euclid_SDSD_vel_post","mean_total_euclid_LDLD_vel_post",
                                   
                                   "mean_total_euclid_SDLD_acc_pre","mean_total_euclid_SDSD_acc_pre","mean_total_euclid_LDLD_acc_pre",
                                   "mean_total_euclid_SDLD_acc_post","mean_total_euclid_SDSD_acc_post","mean_total_euclid_LDLD_acc_post",
                                   
                                   "PC1_se_dist_SDLD_ori_pre","PC2_se_dist_SDLD_ori_pre","PC3_se_dist_SDLD_ori_pre",
                                   "PC1_se_dist_SDSD_ori_pre","PC2_se_dist_SDSD_ori_pre","PC3_se_dist_SDSD_ori_pre",
                                   "PC1_se_dist_LDLD_ori_pre","PC2_se_dist_LDLD_ori_pre","PC3_se_dist_LDLD_ori_pre",
                                   "PC1_se_dist_SDLD_ori_post","PC2_se_dist_SDLD_ori_post","PC3_se_dist_SDLD_ori_post",
                                   "PC1_se_dist_SDSD_ori_post","PC2_se_dist_SDSD_ori_post","PC3_se_dist_SDSD_ori_post",
                                   "PC1_se_dist_LDLD_ori_post","PC2_se_dist_LDLD_ori_post","PC3_se_dist_LDLD_ori_post",
                                   
                                   "PC1_se_dist_SDLD_vel_pre","PC2_se_dist_SDLD_vel_pre","PC3_se_dist_SDLD_vel_pre",
                                   "PC1_se_dist_SDSD_vel_pre","PC2_se_dist_SDSD_vel_pre","PC3_se_dist_SDSD_vel_pre",
                                   "PC1_se_dist_LDLD_vel_pre","PC2_se_dist_LDLD_vel_pre","PC3_se_dist_LDLD_vel_pre",
                                   "PC1_se_dist_SDLD_vel_post","PC2_se_dist_SDLD_vel_post","PC3_se_dist_SDLD_vel_post",
                                   "PC1_se_dist_SDSD_vel_post","PC2_se_dist_SDSD_vel_post","PC3_se_dist_SDSD_vel_post",
                                   "PC1_se_dist_LDLD_vel_post","PC2_se_dist_LDLD_vel_post","PC3_se_dist_LDLD_vel_post",
                                   
                                   "PC1_se_dist_SDLD_acc_pre","PC2_se_dist_SDLD_acc_pre","PC3_se_dist_SDLD_acc_pre",
                                   "PC1_se_dist_SDSD_acc_pre","PC2_se_dist_SDSD_acc_pre","PC3_se_dist_SDSD_acc_pre",
                                   "PC1_se_dist_LDLD_acc_pre","PC2_se_dist_LDLD_acc_pre","PC3_se_dist_LDLD_acc_pre",
                                   "PC1_se_dist_SDLD_acc_post","PC2_se_dist_SDLD_acc_post","PC3_se_dist_SDLD_acc_post",
                                   "PC1_se_dist_SDSD_acc_post","PC2_se_dist_SDSD_acc_post","PC3_se_dist_SDSD_acc_post",
                                   "PC1_se_dist_LDLD_acc_post","PC2_se_dist_LDLD_acc_post","PC3_se_dist_LDLD_acc_post",
                                   
                                   "se_total_euclid_SDLD_ori_pre","se_total_euclid_SDSD_ori_pre","se_total_euclid_LDLD_ori_pre",
                                   "se_total_euclid_SDLD_ori_post","se_total_euclid_SDSD_ori_post","se_total_euclid_LDLD_ori_post",
                                   
                                   "se_total_euclid_SDLD_vel_pre","se_total_euclid_SDSD_vel_pre","se_total_euclid_LDLD_vel_pre",
                                   "se_total_euclid_SDLD_vel_post","se_total_euclid_SDSD_vel_post","se_total_euclid_LDLD_vel_post",
                                   
                                   "se_total_euclid_SDLD_acc_pre","se_total_euclid_SDSD_acc_pre","se_total_euclid_LDLD_acc_pre",
                                   "se_total_euclid_SDLD_acc_post","se_total_euclid_SDSD_acc_post","se_total_euclid_LDLD_acc_post",
                                   
                                   
                                   "PC1_movement_ori_SDLD"," PC2_movement_ori_SDLD ","PC1_movement_vel_SDLD ","PC2_movement_vel_SDLD ","PC3_movement_vel_SDLD ","PC1_movement_acc_SDLD   
      ","PC2_movement_acc_SDLD","PC3_movement_acc_SDLD ","PC1_movement_ori_SDSD","PC2_movement_ori_SDSD ","PC1_movement_vel_SDSD","PC2_movement_vel_SDSD ","PC3_movement_vel_SDSD ","PC1_movement_acc_SDSD   ","PC2_movement_acc_SDSD","PC3_movement_acc_SDSD","PC1_movement_ori_LDLD  ","PC2_movement_ori_LDLD  ","PC1_movement_vel_LDLD","PC2_movement_vel_LDLD ","PC3_movement_vel_LDLD ","PC1_movement_acc_LDLD  ","PC2_movement_acc_LDLD","PC3_movement_acc_LDLD",  
                                   
                                   "diff_SDLD_ori_centerdist","diff_SDSD_ori_centerdist","diff_LDLD_ori_centerdist","diff_SDLD_vel_centerdist","diff_SDSD_vel_centerdist","diff_LDLD_vel_centerdist","diff_SDLD_acc_centerdist","diff_SDSD_acc_centerdist","diff_LDLD_acc_centerdist","ratio_SDLD_ori_centerdist","ratio_SDSD_ori_centerdist","ratio_LDLD_ori_centerdist","ratio_SDLD_vel_centerdist","ratio_SDSD_vel_centerdist","ratio_LDLD_vel_centerdist","ratio_SDLD_acc_centerdist","ratio_SDSD_acc_centerdist","ratio_LDLD_acc_centerdist",
                                   
                                   "change_in_mean_dist_prepost_SDLD_ori_PC1","change_in_mean_dist_prepost_SDLD_ori_PC2","change_in_mean_dist_prepost_SDLD_ori_PC3",
                                   "change_in_mean_dist_prepost_SDSD_ori_PC1","change_in_mean_dist_prepost_SDSD_ori_PC2","change_in_mean_dist_prepost_SDSD_ori_PC3",
                                   "change_in_mean_dist_prepost_LDLD_ori_PC1","change_in_mean_dist_prepost_LDLD_ori_PC2","change_in_mean_dist_prepost_LDLD_ori_PC3",
                                   "change_in_mean_dist_prepost_SDLD_vel_PC1","change_in_mean_dist_prepost_SDLD_vel_PC2","change_in_mean_dist_prepost_SDLD_vel_PC3",
                                   "change_in_mean_dist_prepost_SDSD_vel_PC1","change_in_mean_dist_prepost_SDSD_vel_PC2","change_in_mean_dist_prepost_SDSD_vel_PC3",
                                   "change_in_mean_dist_prepost_LDLD_vel_PC1","change_in_mean_dist_prepost_LDLD_vel_PC2","change_in_mean_dist_prepost_LDLD_vel_PC3",
                                   "change_in_mean_dist_prepost_SDLD_acc_PC1","change_in_mean_dist_prepost_SDLD_acc_PC2","change_in_mean_dist_prepost_SDLD_acc_PC3",
                                   "change_in_mean_dist_prepost_SDSD_acc_PC1","change_in_mean_dist_prepost_SDSD_acc_PC2","change_in_mean_dist_prepost_SDSD_acc_PC3",
                                   "change_in_mean_dist_prepost_LDLD_acc_PC1","change_in_mean_dist_prepost_LDLD_acc_PC2","change_in_mean_dist_prepost_LDLD_acc_PC3",
                                   "change_in_mean_total_euclid_ori_SDLD","change_in_mean_total_euclid_ori_SDSD","change_in_mean_total_euclid_ori_LDLD",
                                   "change_in_se_dist_prepost_SDLD_ori_PC1","change_in_se_dist_prepost_SDLD_ori_PC2","change_in_se_dist_prepost_SDLD_ori_PC3",
                                   "change_in_se_dist_prepost_SDSD_ori_PC1","change_in_se_dist_prepost_SDSD_ori_PC2","change_in_se_dist_prepost_SDSD_ori_PC3",
                                   "change_in_se_dist_prepost_LDLD_ori_PC1","change_in_se_dist_prepost_LDLD_ori_PC2","change_in_se_dist_prepost_LDLD_ori_PC3",
                                   "change_in_se_dist_prepost_SDLD_vel_PC1","change_in_se_dist_prepost_SDLD_vel_PC2","change_in_se_dist_prepost_SDLD_vel_PC3",
                                   "change_in_se_dist_prepost_SDSD_vel_PC1","change_in_se_dist_prepost_SDSD_vel_PC2","change_in_se_dist_prepost_SDSD_vel_PC3",
                                   "change_in_se_dist_prepost_LDLD_vel_PC1","change_in_se_dist_prepost_LDLD_vel_PC2","change_in_se_dist_prepost_LDLD_vel_PC3",
                                   "change_in_se_dist_prepost_SDLD_acc_PC1","change_in_se_dist_prepost_SDLD_acc_PC2","change_in_se_dist_prepost_SDLD_acc_PC3",
                                   "change_in_se_dist_prepost_SDSD_acc_PC1","change_in_se_dist_prepost_SDSD_acc_PC2","change_in_se_dist_prepost_SDSD_acc_PC3",
                                   "change_in_se_dist_prepost_LDLD_acc_PC1","change_in_se_dist_prepost_LDLD_acc_PC2","change_in_se_dist_prepost_LDLD_acc_PC3",
                                   "change_in_se_total_euclid_ori_SDLD","change_in_se_total_euclid_ori_SDSD","change_in_se_total_euclid_ori_LDLD")


##
#save file
write.csv(table_of_phenotypes, "table_of_phenotypes_pre_post.csv")


#

###########################################################################################################



##################### the same but with varimax rotation 

Total_curves_ALL_rmout <- apply(Total_curves_ALL_rmout,2, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))

Total_curves_ALL_rmout <- data.matrix(Total_curves_ALL_rmout)

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 100)

fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
smooth_both <- smooth.basis(argvals = new_ZT, y = Total_curves_ALL_rmout, fdParobj = fdobj)


ori_obj <- deriv.fd(smooth_both$fd,0)
vel_obj <- deriv.fd(smooth_both$fd,1)
acc_obj <- deriv.fd(smooth_both$fd,2)


#FPCA

fpca_ori <- pca.fd(ori_obj,3,harmfdPar = fdPar(ori_obj),centerfns = FALSE)
fpca_ori <- varmx.pca.fd(fpca_ori)
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
Total_curves_SD_LD_rmout <- data.frame(curves_SDLD_pre)
Total_curves_SD_SD_rmout <- data.frame(curves_SDSD_pre)
Total_curves_LD_LD_rmout <- data.frame(curves_LDLD_pre)



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

#my_vector_SDLD <- my_vector_SDLD[order(names(my_vector_SDLD))]
#my_vector_SDSD <- my_vector_SDSD[order(names(my_vector_SDSD))]
#my_vector_LDLD <- my_vector_LDLD[order(names(my_vector_LDLD))]

#combine twice because its the pre section and then the post section

my_vector <- c(my_vector_SDLD,my_vector_SDSD,my_vector_LDLD,my_vector_SDLD,my_vector_SDSD,my_vector_LDLD)

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

#could look at the movement in pc1 and pc2 from pre to post in genotype centers 
#1:64 is SDLD pre
#65:128 is SDSD pre
#129:192 is LDLD pre
#193:256 is SDLD post
#257:320 is SDSD post
#321:384 is LDLD post

PC1_movement_ori_SDLD <- center_coords_ori[193:256,1] - center_coords_ori[1:64,1]  
PC2_movement_ori_SDLD <- center_coords_ori[193:256,2] - center_coords_ori[1:64,2] 
#vel
PC1_movement_vel_SDLD <- center_coords_vel[193:256,1] - center_coords_vel[1:64,1]  
PC2_movement_vel_SDLD <- center_coords_vel[193:256,2] - center_coords_vel[1:64,2] 
PC3_movement_vel_SDLD <- center_coords_vel[193:256,3] - center_coords_vel[1:64,3] 
#acc
PC1_movement_acc_SDLD <- center_coords_acc[193:256,1] - center_coords_acc[1:64,1]  
PC2_movement_acc_SDLD <- center_coords_acc[193:256,2] - center_coords_acc[1:64,2] 
PC3_movement_acc_SDLD <- center_coords_acc[193:256,3] - center_coords_acc[1:64,3] 

#SDSD
PC1_movement_ori_SDSD <- center_coords_ori[257:320,1] - center_coords_ori[65:128,1]  
PC2_movement_ori_SDSD <- center_coords_ori[257:320,2] - center_coords_ori[65:128,2] 
#vel
PC1_movement_vel_SDSD <- center_coords_vel[257:320,1] - center_coords_vel[65:128,1]  
PC2_movement_vel_SDSD <- center_coords_vel[257:320,2] - center_coords_vel[65:128,2] 
PC3_movement_vel_SDSD <- center_coords_vel[257:320,3] - center_coords_vel[65:128,3] 
#acc
PC1_movement_acc_SDSD <- center_coords_acc[257:320,1] - center_coords_acc[65:128,1]  
PC2_movement_acc_SDSD <- center_coords_acc[257:320,2] - center_coords_acc[65:128,2] 
PC3_movement_acc_SDSD <- center_coords_acc[257:320,3] - center_coords_acc[65:128,3] 


#LDLD
PC1_movement_ori_LDLD <- center_coords_ori[321:384,1] - center_coords_ori[129:192,1]  
PC2_movement_ori_LDLD <- center_coords_ori[321:384,2] - center_coords_ori[129:192,2] 
#vel
PC1_movement_vel_LDLD <- center_coords_vel[321:384,1] - center_coords_vel[129:192,1]  
PC2_movement_vel_LDLD <- center_coords_vel[321:384,2] - center_coords_vel[129:192,2] 
PC3_movement_vel_LDLD <- center_coords_vel[321:384,3] - center_coords_vel[129:192,3] 
#acc
PC1_movement_acc_LDLD <- center_coords_acc[321:384,1] - center_coords_acc[129:192,1]  
PC2_movement_acc_LDLD <- center_coords_acc[321:384,2] - center_coords_acc[129:192,2] 
PC3_movement_acc_LDLD <- center_coords_acc[321:384,3] - center_coords_acc[129:192,3] 

#first eculidean distance on total

#ori curves
euclid_dist_pre_post_ori_t <- matrix(ncol = 1, nrow = (length(my_vector)/2))
for (i in 1:(length(my_vector)/2)) {
  
  euclid_dist_pre_post_ori_t[i,] <- sqrt(((center_coords_ori[i,1] - center_coords_ori[i +192 ,1])^2)+((center_coords_ori[i,2] - center_coords_ori[i+192,2])^2)+((center_coords_ori[i,3] - center_coords_ori[i+192,3])^2))
  
}

#vel curves
euclid_dist_pre_post_vel_t <- matrix(ncol = 1, nrow = (length(my_vector)/2))
for (i in 1:(length(my_vector)/2)) {
  
  euclid_dist_pre_post_vel_t[i,] <- sqrt(((center_coords_vel[i,1] - center_coords_vel[i +192 ,1])^2)+((center_coords_vel[i,2] - center_coords_vel[i+192,2])^2)+((center_coords_vel[i,3] - center_coords_vel[i+192,3])^2))
  
}

#acc curves
euclid_dist_pre_post_acc_t <- matrix(ncol = 1, nrow = (length(my_vector)/2))
for (i in 1:(length(my_vector)/2)) {
  
  euclid_dist_pre_post_acc_t[i,] <- sqrt(((center_coords_acc[i,1] - center_coords_acc[i +192 ,1])^2)+((center_coords_acc[i,2] - center_coords_acc[i+192,2])^2)+((center_coords_acc[i,3] - center_coords_acc[i+192,3])^2))
  
}



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


#####Eulcidiena distance total of going to the genotype centers
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

#################################################



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



#need to compare the before and after for the mean distance and spread and the se of distance and spread.... 

#do before take away the after 


change_in_mean_dist_ori <- spread_coords_mean_ori[1:192,1:3] - spread_coords_mean_ori[193:384,1:3]
change_in_mean_dist_vel <- spread_coords_mean_vel[1:192,1:3] - spread_coords_mean_vel[193:384,1:3]
change_in_mean_dist_acc <- spread_coords_mean_acc[1:192,1:3] - spread_coords_mean_acc[193:384,1:3]

change_in_mean_total_euclid_ori <- spread_coords_mean_ori_total[1:192,1] - spread_coords_mean_ori_total[193:384,1]
change_in_mean_total_euclid_vel <- spread_coords_mean_vel_total[1:192,1] - spread_coords_mean_vel_total[193:384,1]
change_in_mean_total_euclid_acc <- spread_coords_mean_acc_total[1:192,1] - spread_coords_mean_acc_total[193:384,1]

change_in_se_dist_ori <- spread_coords_se_ori[1:192,1:3] - spread_coords_se_ori[193:384,1:3]
change_in_se_dist_vel <- spread_coords_se_vel[1:192,1:3] - spread_coords_se_vel[193:384,1:3]
change_in_se_dist_acc <- spread_coords_se_acc[1:192,1:3] - spread_coords_se_acc[193:384,1:3]

change_in_se_total_euclid_ori <- spread_coords_se_ori_total[1:192,1] - spread_coords_se_ori_total[193:384,1]
change_in_se_total_euclid_vel <- spread_coords_se_vel_total[1:192,1] - spread_coords_se_vel_total[193:384,1]
change_in_se_total_euclid_acc <- spread_coords_se_acc_total[1:192,1] - spread_coords_se_acc_total[193:384,1]




#look at the difference in the spread from the pre and post shift (mostly in SD_LD)
#origional
#if posotive means after shift is LESS spread out


diff_SDLD_ori_centerdist <- spread_coords_mean_ori[1:64,1] - spread_coords_mean_ori[193:256,1]

diff_SDSD_ori_centerdist <- spread_coords_mean_ori[65:128,1] - spread_coords_mean_ori[257:320,1]

diff_LDLD_ori_centerdist <- spread_coords_mean_ori[129:192,1] - spread_coords_mean_ori[321:384,1]

ratio_SDLD_ori_centerdist <- spread_coords_mean_ori[1:64,1] / spread_coords_mean_ori[193:256,1]

ratio_SDSD_ori_centerdist <- spread_coords_mean_ori[65:128,1] / spread_coords_mean_ori[257:320,1]

ratio_LDLD_ori_centerdist <- spread_coords_mean_ori[129:192,1] / spread_coords_mean_ori[321:384,1]


#velocity

diff_SDLD_vel_centerdist <- spread_coords_mean_vel[1:64,1] - spread_coords_mean_vel[193:256,1]

diff_SDSD_vel_centerdist <- spread_coords_mean_vel[65:128,1] - spread_coords_mean_vel[257:320,1]

diff_LDLD_vel_centerdist <- spread_coords_mean_vel[129:192,1] - spread_coords_mean_vel[321:384,1]
#try the ratio of spread so do pre/post 

ratio_SDLD_vel_centerdist <- spread_coords_mean_vel[1:64,1] / spread_coords_mean_vel[193:256,1]

ratio_SDSD_vel_centerdist <- spread_coords_mean_vel[65:128,1] / spread_coords_mean_vel[257:320,1]

ratio_LDLD_vel_centerdist <- spread_coords_mean_vel[129:192,1] / spread_coords_mean_vel[321:384,1]



#acceleration

diff_SDLD_acc_centerdist <- spread_coords_mean_acc[1:64,1] - spread_coords_mean_acc[193:256,1]
diff_SDSD_acc_centerdist <- spread_coords_mean_acc[65:128,1] - spread_coords_mean_acc[257:320,1]

diff_LDLD_acc_centerdist <- spread_coords_mean_acc[129:192,1] - spread_coords_mean_acc[321:384,1]

#try the ratio of spread so do pre/post 

ratio_SDLD_acc_centerdist <- spread_coords_mean_acc[1:64,1] / spread_coords_mean_acc[193:256,1]

ratio_SDSD_acc_centerdist <- spread_coords_mean_acc[65:128,1] / spread_coords_mean_acc[257:320,1]

ratio_LDLD_acc_centerdist <- spread_coords_mean_acc[129:192,1] / spread_coords_mean_acc[321:384,1]


################################################################

################################################################

#        TRAIT TABLE

################################################################

################################################################

#make a table with each of the results


table_of_phenotypes_varmx <- cbind(center_coords_ori[1:64,1],center_coords_ori[1:64,2],center_coords_ori[1:64,3],
                                   center_coords_ori[65:128,1],center_coords_ori[65:128,2],center_coords_ori[65:128,3],
                                   center_coords_ori[129:192,1],center_coords_ori[129:192,2],center_coords_ori[129:192,3],
                                   center_coords_ori[193:256,1],center_coords_ori[193:256,2],center_coords_ori[193:256,3],
                                   center_coords_ori[257:320,1],center_coords_ori[257:320,2],center_coords_ori[257:320,3],
                                   center_coords_ori[321:384,1],center_coords_ori[321:384,2],center_coords_ori[321:384,3],
                                   center_coords_vel[1:64,1],center_coords_vel[1:64,2],center_coords_vel[1:64,3],
                                   center_coords_vel[65:128,1],center_coords_vel[65:128,2],center_coords_vel[65:128,3],
                                   center_coords_vel[129:192,1],center_coords_vel[129:192,2],center_coords_vel[129:192,3],
                                   center_coords_vel[193:256,1],center_coords_vel[193:256,2],center_coords_vel[193:256,3],
                                   center_coords_vel[257:320,1],center_coords_vel[257:320,2],center_coords_vel[257:320,3],
                                   center_coords_vel[321:384,1],center_coords_vel[321:384,2],center_coords_vel[321:384,3],
                                   center_coords_acc[1:64,1],center_coords_acc[1:64,2],center_coords_acc[1:64,3],
                                   center_coords_acc[65:128,1],center_coords_acc[65:128,2],center_coords_acc[65:128,3],
                                   center_coords_acc[129:192,1],center_coords_acc[129:192,2],center_coords_acc[129:192,3],
                                   center_coords_acc[193:256,1],center_coords_acc[193:256,2],center_coords_acc[193:256,3],
                                   center_coords_acc[257:320,1],center_coords_acc[257:320,2],center_coords_acc[257:320,3],
                                   center_coords_acc[321:384,1],center_coords_acc[321:384,2],center_coords_acc[321:384,3],euclid_dist_pre_post_ori_t[1:64],euclid_dist_pre_post_ori_t[65:128],euclid_dist_pre_post_ori_t[129:192],euclid_dist_pre_post_vel_t[1:64],euclid_dist_pre_post_vel_t[65:128],euclid_dist_pre_post_vel_t[129:192],euclid_dist_pre_post_acc_t[1:64],euclid_dist_pre_post_acc_t[65:128],euclid_dist_pre_post_acc_t[129:192],
                                   spread_coords_mean_ori[1:64,1],spread_coords_mean_ori[1:64,2],spread_coords_mean_ori[1:64,3],
                                   spread_coords_mean_ori[65:128,1],spread_coords_mean_ori[65:128,2],spread_coords_mean_ori[65:128,3],
                                   spread_coords_mean_ori[129:192,1],spread_coords_mean_ori[129:192,2],spread_coords_mean_ori[129:192,3],
                                   spread_coords_mean_ori[193:256,1],spread_coords_mean_ori[193:256,2],spread_coords_mean_ori[193:256,3],
                                   spread_coords_mean_ori[257:320,1],spread_coords_mean_ori[257:320,2],spread_coords_mean_ori[257:320,3],
                                   spread_coords_mean_ori[321:384,1],spread_coords_mean_ori[321:384,2],spread_coords_mean_ori[321:384,3],
                                   
                                   spread_coords_mean_vel[1:64,1],spread_coords_mean_vel[1:64,2],spread_coords_mean_vel[1:64,3],
                                   spread_coords_mean_vel[65:128,1],spread_coords_mean_vel[65:128,2],spread_coords_mean_vel[65:128,3],
                                   spread_coords_mean_vel[129:192,1],spread_coords_mean_vel[129:192,2],spread_coords_mean_vel[129:192,3],
                                   spread_coords_mean_vel[193:256,1],spread_coords_mean_vel[193:256,2],spread_coords_mean_vel[193:256,3],
                                   spread_coords_mean_vel[257:320,1],spread_coords_mean_vel[257:320,2],spread_coords_mean_vel[257:320,3],
                                   spread_coords_mean_vel[321:384,1],spread_coords_mean_vel[321:384,2],spread_coords_mean_vel[321:384,3],
                                   
                                   spread_coords_mean_acc[1:64,1],spread_coords_mean_acc[1:64,2],spread_coords_mean_acc[1:64,3],
                                   spread_coords_mean_acc[65:128,1],spread_coords_mean_acc[65:128,2],spread_coords_mean_acc[65:128,3],
                                   spread_coords_mean_acc[129:192,1],spread_coords_mean_acc[129:192,2],spread_coords_mean_acc[129:192,3],
                                   spread_coords_mean_acc[193:256,1],spread_coords_mean_acc[193:256,2],spread_coords_mean_acc[193:256,3],
                                   spread_coords_mean_acc[257:320,1],spread_coords_mean_acc[257:320,2],spread_coords_mean_acc[257:320,3],
                                   spread_coords_mean_acc[321:384,1],spread_coords_mean_acc[321:384,2],spread_coords_mean_acc[321:384,3],
                                   spread_coords_mean_ori_total[1:64,1],
                                   spread_coords_mean_ori_total[65:128,1],
                                   spread_coords_mean_ori_total[129:192,1],
                                   spread_coords_mean_ori_total[193:256,1],
                                   spread_coords_mean_ori_total[257:320,1],
                                   spread_coords_mean_ori_total[321:384,1],
                                   spread_coords_mean_vel_total[1:64,1],
                                   spread_coords_mean_vel_total[65:128,1],
                                   spread_coords_mean_vel_total[129:192,1],
                                   spread_coords_mean_vel_total[193:256,1],
                                   spread_coords_mean_vel_total[257:320,1],
                                   spread_coords_mean_vel_total[321:384,1],
                                   spread_coords_mean_acc_total[1:64,1],
                                   spread_coords_mean_acc_total[65:128,1],
                                   spread_coords_mean_acc_total[129:192,1],
                                   spread_coords_mean_acc_total[193:256,1],
                                   spread_coords_mean_acc_total[257:320,1],
                                   spread_coords_mean_acc_total[321:384,1],
                                   
                                   
                                   spread_coords_se_ori[1:64,1],spread_coords_se_ori[1:64,2],spread_coords_se_ori[1:64,3],
                                   spread_coords_se_ori[65:128,1],spread_coords_se_ori[65:128,2],spread_coords_se_ori[65:128,3],
                                   spread_coords_se_ori[129:192,1],spread_coords_se_ori[129:192,2],spread_coords_se_ori[129:192,3],
                                   spread_coords_se_ori[193:256,1],spread_coords_se_ori[193:256,2],spread_coords_se_ori[193:256,3],
                                   spread_coords_se_ori[257:320,1],spread_coords_se_ori[257:320,2],spread_coords_se_ori[257:320,3],
                                   spread_coords_se_ori[321:384,1],spread_coords_se_ori[321:384,2],spread_coords_se_ori[321:384,3],
                                   
                                   spread_coords_se_vel[1:64,1],spread_coords_se_vel[1:64,2],spread_coords_se_vel[1:64,3],
                                   spread_coords_se_vel[65:128,1],spread_coords_se_vel[65:128,2],spread_coords_se_vel[65:128,3],
                                   spread_coords_se_vel[129:192,1],spread_coords_se_vel[129:192,2],spread_coords_se_vel[129:192,3],
                                   spread_coords_se_vel[193:256,1],spread_coords_se_vel[193:256,2],spread_coords_se_vel[193:256,3],
                                   spread_coords_se_vel[257:320,1],spread_coords_se_vel[257:320,2],spread_coords_se_vel[257:320,3],
                                   spread_coords_se_vel[321:384,1],spread_coords_se_vel[321:384,2],spread_coords_se_vel[321:384,3],
                                   
                                   spread_coords_se_acc[1:64,1],spread_coords_se_acc[1:64,2],spread_coords_se_acc[1:64,3],
                                   spread_coords_se_acc[65:128,1],spread_coords_se_acc[65:128,2],spread_coords_se_acc[65:128,3],
                                   spread_coords_se_acc[129:192,1],spread_coords_se_acc[129:192,2],spread_coords_se_acc[129:192,3],
                                   spread_coords_se_acc[193:256,1],spread_coords_se_acc[193:256,2],spread_coords_se_acc[193:256,3],
                                   spread_coords_se_acc[257:320,1],spread_coords_se_acc[257:320,2],spread_coords_se_acc[257:320,3],
                                   spread_coords_se_acc[321:384,1],spread_coords_se_acc[321:384,2],spread_coords_se_acc[321:384,3],
                                   spread_coords_se_ori_total[1:64,1],
                                   spread_coords_se_ori_total[65:128,1],
                                   spread_coords_se_ori_total[129:192,1],
                                   spread_coords_se_ori_total[193:256,1],
                                   spread_coords_se_ori_total[257:320,1],
                                   spread_coords_se_ori_total[321:384,1],
                                   spread_coords_se_vel_total[1:64,1],
                                   spread_coords_se_vel_total[65:128,1],
                                   spread_coords_se_vel_total[129:192,1],
                                   spread_coords_se_vel_total[193:256,1],
                                   spread_coords_se_vel_total[257:320,1],
                                   spread_coords_se_vel_total[321:384,1],
                                   spread_coords_se_acc_total[1:64,1],
                                   spread_coords_se_acc_total[65:128,1],
                                   spread_coords_se_acc_total[129:192,1],
                                   spread_coords_se_acc_total[193:256,1],
                                   spread_coords_se_acc_total[257:320,1],
                                   spread_coords_se_acc_total[321:384,1],
                                   
                                   PC1_movement_ori_SDLD,PC2_movement_ori_SDLD ,PC1_movement_vel_SDLD ,PC2_movement_vel_SDLD ,PC3_movement_vel_SDLD ,PC1_movement_acc_SDLD   
                                   ,PC2_movement_acc_SDLD,PC3_movement_acc_SDLD ,PC1_movement_ori_SDSD,PC2_movement_ori_SDSD ,PC1_movement_vel_SDSD ,PC2_movement_vel_SDSD ,PC3_movement_vel_SDSD ,PC1_movement_acc_SDSD   ,PC2_movement_acc_SDSD ,PC3_movement_acc_SDSD,PC1_movement_ori_LDLD  ,PC2_movement_ori_LDLD  ,PC1_movement_vel_LDLD ,PC2_movement_vel_LDLD ,PC3_movement_vel_LDLD ,PC1_movement_acc_LDLD  ,PC2_movement_acc_LDLD  ,PC3_movement_acc_LDLD,                
                                   diff_SDLD_ori_centerdist,diff_SDSD_ori_centerdist,diff_LDLD_ori_centerdist,diff_SDLD_vel_centerdist,diff_SDSD_vel_centerdist,diff_LDLD_vel_centerdist,diff_SDLD_acc_centerdist,diff_SDSD_acc_centerdist,diff_LDLD_acc_centerdist,ratio_SDLD_ori_centerdist,ratio_SDSD_ori_centerdist,ratio_LDLD_ori_centerdist,ratio_SDLD_vel_centerdist,ratio_SDSD_vel_centerdist,ratio_LDLD_vel_centerdist,ratio_SDLD_acc_centerdist,ratio_SDSD_acc_centerdist,ratio_LDLD_acc_centerdist,change_in_mean_dist_ori[1:64,1:3],change_in_mean_dist_ori[65:128,1:3],change_in_mean_dist_ori[129:192,1:3],
                                   change_in_mean_dist_vel[1:64,1:3],change_in_mean_dist_vel[65:128,1:3],change_in_mean_dist_vel[129:192,1:3],
                                   change_in_mean_dist_acc[1:64,1:3],change_in_mean_dist_acc[65:128,1:3],change_in_mean_dist_acc[129:192,1:3],
                                   change_in_mean_total_euclid_ori[1:64],change_in_mean_total_euclid_ori[65:128],change_in_mean_total_euclid_ori[129:192],
                                   change_in_se_dist_ori[1:64,1:3],change_in_se_dist_ori[65:128,1:3],change_in_se_dist_ori[129:192,1:3],
                                   change_in_se_dist_vel[1:64,1:3],change_in_se_dist_vel[65:128,1:3],change_in_se_dist_vel[129:192,1:3],
                                   change_in_se_dist_acc[1:64,1:3],change_in_se_dist_acc[65:128,1:3],change_in_se_dist_acc[129:192,1:3],
                                   change_in_se_total_euclid_ori[1:64],change_in_se_total_euclid_ori[65:128],change_in_se_total_euclid_ori[129:192])

#
row.names(table_of_phenotypes_varmx) <- sort(Genotype_names)

colnames(table_of_phenotypes_varmx)<-c("varmx_PC1_gc_SDLD_pre","varmx_PC2_gc_SDLD_pre","varmx_PC3_gc_SDLD_pre",
                                       "varmx_PC1_gc_SDSD_pre","varmx_PC2_gc_SDSD_pre","varmx_PC3_gc_SDSD_pre","varmx_PC1_gc_LDLD_pre","varmx_PC2_gc_LDLD_pre","varmx_PC3_gc_LDLD_pre","varmx_PC1_gc_SDLD_post","varmx_PC2_gc_SDLD_post","varmx_PC3_gc_SDLD_post","varmx_PC1_gc_SDSD_post","varmx_PC2_gc_SDSD_post","varmx_PC3_gc_SDSD_post","varmx_PC1_gc_LDLD_post","varmx_PC2_gc_LDLD_post","varmx_PC3_gc_LDLD_post","varmx_PC1_gc_vel_SDLD_pre","varmx_PC2_gc_vel_SDLD_pre","varmx_PC3_gc_vel_SDLD_pre","varmx_PC1_gc_vel_SDSD_pre","varmx_PC2_gc_vel_SDSD_pre","varmx_PC3_gc_vel_SDSD_pre","varmx_PC1_gc_vel_LDLD_pre","varmx_PC2_gc_vel_LDLD_pre","varmx_PC3_gc_vel_LDLD_pre","varmx_PC1_gc_vel_SDLD_post","varmx_PC2_gc_vel_SDLD_post","varmx_PC3_gc_vel_SDLD_post","varmx_PC1_gc_vel_SDSD_post","varmx_PC2_gc_vel_SDSD_post","varmx_PC3_gc_vel_SDSD_post","varmx_PC1_gc_vel_LDLD_post","varmx_PC2_gc_vel_LDLD_post","varmx_PC3_gc_vel_LDLD_post","varmx_PC1_gc_acc_SDLD_pre","varmx_PC2_gc_acc_SDLD_pre","varmx_PC3_gc_acc_SDLD_pre","varmx_PC1_gc_acc_SDSD_pre","varmx_PC2_gc_acc_SDSD_pre","varmx_PC3_gc_acc_SDSD_pre","varmx_PC1_gc_acc_LDLD_pre","varmx_PC2_gc_acc_LDLD_pre","varmx_PC3_gc_acc_LDLD_pre","varmx_PC1_gc_acc_SDLD_post","varmx_PC2_gc_acc_SDLD_post","varmx_PC3_gc_acc_SDLD_post","varmx_PC1_gc_acc_SDSD_post","varmx_PC2_gc_acc_SDSD_post","varmx_PC3_gc_acc_SDSD_post","varmx_PC1_gc_acc_LDLD_post","varmx_PC2_gc_acc_LDLD_post",
                                       
                                       "varmx_PC3_gc_acc_LDLD_post","varmx_euclid_dist_pre_post_ori_tSDLD","varmx_euclid_dist_pre_post_ori_tSDSD","varmx_euclid_dist_pre_post_ori_tLDLD","varmx_euclid_dist_pre_post_vel_tSDLD","varmx_euclid_dist_pre_post_vel_tSDSD","varmx_euclid_dist_pre_post_vel_tLDLD","varmx_euclid_dist_pre_post_acc_tSDLD","varmx_euclid_dist_pre_post_acc_tSDSD","varmx_euclid_dist_pre_post_acc_tLDLD","varmx_PC1_mean_dist_SDLD_ori_pre","varmx_PC2_mean_dist_SDLD_ori_pre","varmx_PC3_mean_dist_SDLD_ori_pre","varmx_PC1_mean_dist_SDSD_ori_pre","varmx_PC2_mean_dist_SDSD_ori_pre","varmx_PC3_mean_dist_SDSD_ori_pre","varmx_PC1_mean_dist_LDLD_ori_pre","varmx_PC2_mean_dist_LDLD_ori_pre","varmx_PC3_mean_dist_LDLD_ori_pre","varmx_PC1_mean_dist_SDLD_ori_post","varmx_PC2_mean_dist_SDLD_ori_post","varmx_PC3_mean_dist_SDLD_ori_post","varmx_PC1_mean_dist_SDSD_ori_post","varmx_PC2_mean_dist_SDSD_ori_post","varmx_PC3_mean_dist_SDSD_ori_post","varmx_PC1_mean_dist_LDLD_ori_post","varmx_PC2_mean_dist_LDLD_ori_post","varmx_PC3_mean_dist_LDLD_ori_post","varmx_PC1_mean_dist_SDLD_vel_pre","varmx_PC2_mean_dist_SDLD_vel_pre","varmx_PC3_mean_dist_SDLD_vel_pre","varmx_PC1_mean_dist_SDSD_vel_pre","varmx_PC2_mean_dist_SDSD_vel_pre",
                                       "varmx_PC3_mean_dist_SDSD_vel_pre","varmx_PC1_mean_dist_LDLD_vel_pre","varmx_PC2_mean_dist_LDLD_vel_pre","varmx_PC3_mean_dist_LDLD_vel_pre","varmx_PC1_mean_dist_SDLD_vel_post","varmx_PC2_mean_dist_SDLD_vel_post","varmx_PC3_mean_dist_SDLD_vel_post","varmx_PC1_mean_dist_SDSD_vel_post","varmx_PC2_mean_dist_SDSD_vel_post","varmx_PC3_mean_dist_SDSD_vel_post","varmx_PC1_mean_dist_LDLD_vel_post","varmx_PC2_mean_dist_LDLD_vel_post","varmx_PC3_mean_dist_LDLD_vel_post","varmx_PC1_mean_dist_SDLD_acc_pre","varmx_PC2_mean_dist_SDLD_acc_pre","varmx_PC3_mean_dist_SDLD_acc_pre","varmx_PC1_mean_dist_SDSD_acc_pre","varmx_PC2_mean_dist_SDSD_acc_pre","varmx_PC3_mean_dist_SDSD_acc_pre","varmx_PC1_mean_dist_LDLD_acc_pre","varmx_PC2_mean_dist_LDLD_acc_pre","varmx_PC3_mean_dist_LDLD_acc_pre","varmx_PC1_mean_dist_SDLD_acc_post","varmx_PC2_mean_dist_SDLD_acc_post","varmx_PC3_mean_dist_SDLD_acc_post","varmx_PC1_mean_dist_SDSD_acc_post","varmx_PC2_mean_dist_SDSD_acc_post","varmx_PC3_mean_dist_SDSD_acc_post","varmx_PC1_mean_dist_LDLD_acc_post","varmx_PC2_mean_dist_LDLD_acc_post","varmx_PC3_mean_dist_LDLD_acc_post",
                                       "varmx_mean_total_euclid_SDLD_ori_pre","varmx_mean_total_euclid_SDSD_ori_pre","varmx_mean_total_euclid_LDLD_ori_pre","varmx_mean_total_euclid_SDLD_ori_post","varmx_mean_total_euclid_SDSD_ori_post","varmx_mean_total_euclid_LDLD_ori_post","varmx_mean_total_euclid_SDLD_vel_pre","varmx_mean_total_euclid_SDSD_vel_pre","varmx_mean_total_euclid_LDLD_vel_pre","varmx_mean_total_euclid_SDLD_vel_post","varmx_mean_total_euclid_SDSD_vel_post","varmx_mean_total_euclid_LDLD_vel_post","varmx_mean_total_euclid_SDLD_acc_pre","varmx_mean_total_euclid_SDSD_acc_pre","varmx_mean_total_euclid_LDLD_acc_pre","varmx_mean_total_euclid_SDLD_acc_post","varmx_mean_total_euclid_SDSD_acc_post","varmx_mean_total_euclid_LDLD_acc_post","varmx_PC1_se_dist_SDLD_ori_pre","varmx_PC2_se_dist_SDLD_ori_pre","varmx_PC3_se_dist_SDLD_ori_pre","varmx_PC1_se_dist_SDSD_ori_pre","varmx_PC2_se_dist_SDSD_ori_pre","varmx_PC3_se_dist_SDSD_ori_pre","varmx_PC1_se_dist_LDLD_ori_pre","varmx_PC2_se_dist_LDLD_ori_pre","varmx_PC3_se_dist_LDLD_ori_pre","varmx_PC1_se_dist_SDLD_ori_post","varmx_PC2_se_dist_SDLD_ori_post","varmx_PC3_se_dist_SDLD_ori_post","varmx_PC1_se_dist_SDSD_ori_post","varmx_PC2_se_dist_SDSD_ori_post","varmx_PC3_se_dist_SDSD_ori_post","varmx_PC1_se_dist_LDLD_ori_post","varmx_PC2_se_dist_LDLD_ori_post","varmx_PC3_se_dist_LDLD_ori_post","varmx_PC1_se_dist_SDLD_vel_pre","varmx_PC2_se_dist_SDLD_vel_pre","varmx_PC3_se_dist_SDLD_vel_pre","varmx_PC1_se_dist_SDSD_vel_pre","varmx_PC2_se_dist_SDSD_vel_pre","varmx_PC3_se_dist_SDSD_vel_pre","varmx_PC1_se_dist_LDLD_vel_pre","varmx_PC2_se_dist_LDLD_vel_pre","varmx_PC3_se_dist_LDLD_vel_pre","varmx_PC1_se_dist_SDLD_vel_post",
                                       
                                       "varmx_PC2_se_dist_SDLD_vel_post","varmx_PC3_se_dist_SDLD_vel_post","varmx_PC1_se_dist_SDSD_vel_post","varmx_PC2_se_dist_SDSD_vel_post","varmx_PC3_se_dist_SDSD_vel_post","varmx_PC1_se_dist_LDLD_vel_post","varmx_PC2_se_dist_LDLD_vel_post","varmx_PC3_se_dist_LDLD_vel_post","varmx_PC1_se_dist_SDLD_acc_pre","varmx_PC2_se_dist_SDLD_acc_pre","varmx_PC3_se_dist_SDLD_acc_pre","varmx_PC1_se_dist_SDSD_acc_pre","varmx_PC2_se_dist_SDSD_acc_pre","varmx_PC3_se_dist_SDSD_acc_pre","varmx_PC1_se_dist_LDLD_acc_pre","varmx_PC2_se_dist_LDLD_acc_pre","varmx_PC3_se_dist_LDLD_acc_pre","varmx_PC1_se_dist_SDLD_acc_post","varmx_PC2_se_dist_SDLD_acc_post","varmx_PC3_se_dist_SDLD_acc_post","varmx_PC1_se_dist_SDSD_acc_post","varmx_PC2_se_dist_SDSD_acc_post","varmx_PC3_se_dist_SDSD_acc_post","varmx_PC1_se_dist_LDLD_acc_post","varmx_PC2_se_dist_LDLD_acc_post","varmx_PC3_se_dist_LDLD_acc_post","varmx_se_total_euclid_SDLD_ori_pre","varmx_se_total_euclid_SDSD_ori_pre","varmx_se_total_euclid_LDLD_ori_pre",
                                       
                                       "varmx_se_total_euclid_SDLD_ori_post","varmx_se_total_euclid_SDSD_ori_post","varmx_se_total_euclid_LDLD_ori_post","varmx_se_total_euclid_SDLD_vel_pre","varmx_se_total_euclid_SDSD_vel_pre","varmx_se_total_euclid_LDLD_vel_pre","varmx_se_total_euclid_SDLD_vel_post","varmx_se_total_euclid_SDSD_vel_post","varmx_se_total_euclid_LDLD_vel_post","varmx_se_total_euclid_SDLD_acc_pre","varmx_se_total_euclid_SDSD_acc_pre","varmx_se_total_euclid_LDLD_acc_pre","varmx_se_total_euclid_SDLD_acc_post","varmx_se_total_euclid_SDSD_acc_post","varmx_se_total_euclid_LDLD_acc_post","varmx_PC1_movement_ori_SDLD","varmx_PC2_movement_ori_SDLD","varmx_PC1_movement_vel_SDLD","varmx_PC2_movement_vel_SDLD","varmx_PC3_movement_vel_SDLD","varmx_PC1_movement_acc_SDLD","varmx_PC2_movement_acc_SDLD","varmx_PC3_movement_acc_SDLD","varmx_PC1_movement_ori_SDSD","varmx_PC2_movement_ori_SDSD","varmx_PC1_movement_vel_SDSD","varmx_PC2_movement_vel_SDSD","varmx_PC3_movement_vel_SDSD","varmx_PC1_movement_acc_SDSD","varmx_PC2_movement_acc_SDSD",
                                       
                                       "varmx_PC3_movement_acc_SDSD","varmx_PC1_movement_ori_LDLD","varmx_PC2_movement_ori_LDLD","varmx_PC1_movement_vel_LDLD","varmx_PC2_movement_vel_LDLD","varmx_PC3_movement_vel_LDLD","varmx_PC1_movement_acc_LDLD","varmx_PC2_movement_acc_LDLD","varmx_PC3_movement_acc_LDLD","varmx_diff_SDLD_ori_centerdist","varmx_diff_SDSD_ori_centerdist","varmx_diff_LDLD_ori_centerdist","varmx_diff_SDLD_vel_centerdist","varmx_diff_SDSD_vel_centerdist","varmx_diff_LDLD_vel_centerdist","varmx_diff_SDLD_acc_centerdist","varmx_diff_SDSD_acc_centerdist","varmx_diff_LDLD_acc_centerdist","varmx_ratio_SDLD_ori_centerdist","varmx_ratio_SDSD_ori_centerdist","varmx_ratio_LDLD_ori_centerdist","varmx_ratio_SDLD_vel_centerdist","varmx_ratio_SDSD_vel_centerdist","varmx_ratio_LDLD_vel_centerdist","varmx_ratio_SDLD_acc_centerdist","varmx_ratio_SDSD_acc_centerdist","varmx_ratio_LDLD_acc_centerdist","varmx_change_in_mean_dist_prepost_SDLD_ori_PC1","varmx_change_in_mean_dist_prepost_SDLD_ori_PC2","varmx_change_in_mean_dist_prepost_SDLD_ori_PC3","varmx_change_in_mean_dist_prepost_SDSD_ori_PC1","varmx_change_in_mean_dist_prepost_SDSD_ori_PC2","varmx_change_in_mean_dist_prepost_SDSD_ori_PC3","varmx_change_in_mean_dist_prepost_LDLD_ori_PC1",
                                       
                                       "varmx_change_in_mean_dist_prepost_LDLD_ori_PC2","varmx_change_in_mean_dist_prepost_LDLD_ori_PC3","varmx_change_in_mean_dist_prepost_SDLD_vel_PC1","varmx_change_in_mean_dist_prepost_SDLD_vel_PC2","varmx_change_in_mean_dist_prepost_SDLD_vel_PC3","varmx_change_in_mean_dist_prepost_SDSD_vel_PC1","varmx_change_in_mean_dist_prepost_SDSD_vel_PC2","varmx_change_in_mean_dist_prepost_SDSD_vel_PC3","varmx_change_in_mean_dist_prepost_LDLD_vel_PC1","varmx_change_in_mean_dist_prepost_LDLD_vel_PC2","varmx_change_in_mean_dist_prepost_LDLD_vel_PC3","varmx_change_in_mean_dist_prepost_SDLD_acc_PC1","varmx_change_in_mean_dist_prepost_SDLD_acc_PC2","varmx_change_in_mean_dist_prepost_SDLD_acc_PC3","varmx_change_in_mean_dist_prepost_SDSD_acc_PC1","varmx_change_in_mean_dist_prepost_SDSD_acc_PC2","varmx_change_in_mean_dist_prepost_SDSD_acc_PC3","varmx_change_in_mean_dist_prepost_LDLD_acc_PC1","varmx_change_in_mean_dist_prepost_LDLD_acc_PC2","varmx_change_in_mean_dist_prepost_LDLD_acc_PC3","varmx_change_in_mean_total_euclid_ori_SDLD","varmx_change_in_mean_total_euclid_ori_SDSD","varmx_change_in_mean_total_euclid_ori_LDLD","varmx_change_in_se_dist_prepost_SDLD_ori_PC1","varmx_change_in_se_dist_prepost_SDLD_ori_PC2","varmx_change_in_se_dist_prepost_SDLD_ori_PC3","varmx_change_in_se_dist_prepost_SDSD_ori_PC1","varmx_change_in_se_dist_prepost_SDSD_ori_PC2","varmx_change_in_se_dist_prepost_SDSD_ori_PC3","varmx_change_in_se_dist_prepost_LDLD_ori_PC1","varmx_change_in_se_dist_prepost_LDLD_ori_PC2","varmx_change_in_se_dist_prepost_LDLD_ori_PC3","varmx_change_in_se_dist_prepost_SDLD_vel_PC1",
                                       
                                       "varmx_change_in_se_dist_prepost_SDLD_vel_PC2","varmx_change_in_se_dist_prepost_SDLD_vel_PC3","varmx_change_in_se_dist_prepost_SDSD_vel_PC1","varmx_change_in_se_dist_prepost_SDSD_vel_PC2","varmx_change_in_se_dist_prepost_SDSD_vel_PC3","varmx_change_in_se_dist_prepost_LDLD_vel_PC1","varmx_change_in_se_dist_prepost_LDLD_vel_PC2","varmx_change_in_se_dist_prepost_LDLD_vel_PC3","varmx_change_in_se_dist_prepost_SDLD_acc_PC1","varmx_change_in_se_dist_prepost_SDLD_acc_PC2","varmx_change_in_se_dist_prepost_SDLD_acc_PC3","varmx_change_in_se_dist_prepost_SDSD_acc_PC1","varmx_change_in_se_dist_prepost_SDSD_acc_PC2","varmx_change_in_se_dist_prepost_SDSD_acc_PC3","varmx_change_in_se_dist_prepost_LDLD_acc_PC1","varmx_change_in_se_dist_prepost_LDLD_acc_PC2","varmx_change_in_se_dist_prepost_LDLD_acc_PC3","varmx_change_in_se_total_euclid_ori_SDLD","varmx_change_in_se_total_euclid_ori_SDSD","varmx_change_in_se_total_euclid_ori_LDLD")

##
#save file
write.csv(table_of_phenotypes_varmx, "table_of_phenotypes_pre_post_varmx.csv")


