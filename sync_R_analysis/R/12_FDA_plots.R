# Script for plots associated with FDA approach..... 


#setwd
setwd("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts")

# librarys used 
library(dplyr)
library(data.table)
library(tidyverse)
library(fda)

#this file looks at the QTL traits for pre and post shift of the experiment for SD_LD, SD_SD and LD_LD

#setwd
setwd("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/")

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
fpca_vel <- pca.fd(vel_obj,3,harmfdPar = fdPar(vel_obj),centerfns = FALSE)
#fpca_vel <- fda::varmx.pca.fd(fpca_vel)

scores <- fpca_vel$scores

rownames(scores) <- colnames(Total_curves_ALL_rmout)

#this plots each of the principle compments
svg("../../../My_Papers/photoperiod_shift/Figures/FPCA_vel_harmonics.svg", width = 7.2, height = 6)

plot(fpca_vel$harmonics, lty = "solid", lwd =2, xlab = "time (h)" )
legend("topright", c("FPC1", "FPC2", "FPC3"), col = c("black", "red", "green"), lty = 1,lwd =2, bty = "n")

dev.off()

fpca_vel$varprop*100
sum(fpca_vel$varprop*100)


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


center_coords_vel <- matrix(ncol = 3, nrow = length(my_vector))
for (i in 1:length(my_vector)) {
  center_coords_vel[i,] <- c(mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],1]),mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],2]),mean(fpca_vel$scores[(1 + (cum[i] - my_vector[i])):cum[i],3]))
}


########################################
#          plot for genotype centers...... 

svg("../../../My_Papers/photoperiod_shift/Figures/gc_plot.svg", width = 7.2, height = 6)

plot(fpca_vel$scores,xlab = 'PC Score 1',ylab = 'PC Score 2',col = "white",
     cex.lab = 1.5,cex.axis = 1.5,cex = 1)
#plot all 1st genotype for 
points(fpca_vel$scores[1:43,1],fpca_vel$scores[1:43,2], pch = 19, col = "grey44")
segments(fpca_vel$scores[1:43,1],fpca_vel$scores[1:43,2],center_coords_vel[1,1],center_coords_vel[1,2], col="grey44")
points(center_coords_vel[1,1],center_coords_vel[1,2], pch = 19, col ="slateblue2", lwd =4)

points(fpca_vel$scores[3492:3535,1],fpca_vel$scores[3492:3535,2], pch = 19, col = "grey44")
segments(fpca_vel$scores[3492:3535,1],fpca_vel$scores[3492:3535,2],center_coords_vel[193,1],center_coords_vel[193,2], col="grey44")
points(center_coords_vel[193,1],center_coords_vel[193,2], pch = 19, col ="darkblue", lwd =4)

dev.off()

# plot for each condidion together... 

par(mar = c(5, 5, 4, 2))
#FOR THE SD_LDPLANTS


svg("../../../My_Papers/photoperiod_shift/Figures/FPCA_SDLD.svg", width = 7.2, height = 6)

plot(center_coords_vel[1:64,1], center_coords_vel[1:64,2],xlab = 'PC Score 1',ylab = 'PC Score 2',
     cex.lab = 1.5,cex.axis = 1.5,cex = 1, col = "slateblue2", pch = 19, xlim = c(0,0.6), ylim = c(-0.4,0.55), main = "SD --> LD")
points(center_coords_vel[193:256,1],center_coords_vel[193:256,2],xlab = 'PC Score 1',ylab = 'PC Score 2',
       cex.lab = 1.5,cex.axis = 1.5,cex = 1, col = "darkblue", pch = 19)
legend("bottomleft",  c("SD_LD pre shift","SD_LD post shift"), col =  c("slateblue2","darkblue"), pch = 19,bty = "n", cex = 1.5)
dev.off()


svg("../../../My_Papers/photoperiod_shift/Figures/FPCA_SDSD.svg", width = 7.2, height = 6)
#FOR THE SD_SD PLANTS
plot(center_coords_vel[65:128,1], center_coords_vel[65:128,2],xlab = 'PC Score 1',ylab = 'PC Score 2',
     cex.lab = 1.5,cex.axis = 1.5,cex = 1, col = "red2", pch = 19, xlim = c(0,0.6), ylim = c(-0.4,0.55), main = "SD --> SD")
points(center_coords_vel[257:320,1],center_coords_vel[257:320,2],xlab = 'PC Score 1',ylab = 'PC Score 2',
       cex.lab = 1.5,cex.axis = 1.5,cex = 1, col = "darkred", pch = 19)
legend("bottomleft",  c("SD_SD pre shift","SD_SD post shift"), col =  c("red2","darkred"), pch = 19,bty="n", cex = 1.5)
dev.off()




#FOR THE LD_LD PLANTS

svg("../../../My_Papers/photoperiod_shift/Figures/FPCA_LDLD.svg", width = 7.2, height = 6)

plot(center_coords_vel[129:192,1], center_coords_vel[129:192,2],xlab = 'PC Score 1',ylab = 'PC Score 2',
     cex.lab = 1.5,cex.axis = 1.5,cex = 1, col = "green3", pch = 19, xlim = c(0,0.6), ylim = c(-0.4,0.55),main = "LD --> LD")
points(center_coords_vel[321:384,1],center_coords_vel[321:384,2],xlab = 'PC Score 1',ylab = 'PC Score 2',
       cex.lab = 1.5,cex.axis = 1.5,cex = 1, col = "darkgreen", pch = 19)
legend("bottomright",  c("LD_LD pre shift","LD_LD post shift"), col =  c("green3","darkgreen"), pch = 19,bty = "n", cex = 1.5)

dev.off()
##########################################



# plot of heatmap trait vs 










#load in catagories 

catagories <- read.csv("catagories_of_linked_traits_all.csv", header = FALSE)

Shape <- catagories[catagories$V2 == 'Shape',]

#then want to filter the functional traits file with these columns the genetic map starts at column 991 so crop 
trait_data <- read.csv("FULL_table_of_traits_map.csv")
#remove first two blank rows 
developmental <- trait_data[-c(1,2),985:990]
functional_data <- trait_data[-c(1,2),2:984]
genotypes <- trait_data[-c(1,2),1]
#then filter functional data for the traits in the shape cataorgie above 

# Extract column names of df1 that are in the first column of df2
filtered_columns <- colnames(functional_data) %in% Shape$V1

# Subset df1 to keep only these columns
filtered_functional_Shape <- functional_data[, filtered_columns]

filtered_functional_Shape <- filtered_functional_Shape[,-c(20,21,81,82,83,84,85,86,87,88)]
# want to find the max and min of each column 
# Find the maximum of each column
max_values <- apply(filtered_functional_Shape, 2, max, na.rm = TRUE)

# Find the minimum of each column
min_values <- apply(filtered_functional_Shape, 2, min, na.rm = TRUE)

both <- cbind(max_values,min_values)

# Calculate the range for each column by subtracting min from max
rangee <- abs(both[,1] - both[,2])

# Find the column with the largest range
which.max(rangee)


#read in mean curves of both the maximum and minamum genotypes of that trait..


#look at varmx_PC1_gc_LDLD genotypes. Add row names from triat table to see what genotypes 
rownames(filtered_functional_Shape) <- genotypes

#genotypes wxT 15 and WxT 58 


# load in the mean curves 

meancurves <- read.csv("genotype_means_LDLD.csv")


#wxT26 is 589:630 cols and Wxt34 is 843:876
wx26<-Total_curves_ALL_rmout[,589:630]

wx261 <- wx26+1
wx34<-Total_curves_ALL_rmout[,843:876]

Total_curves_ALL_rmout <- cbind(wx26,wx34)

BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_ZT),max(new_ZT)),norder = 4,nbasis = 300)

#using high smothing to get overall shape of curves 
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 5)
sp_totalsmooth <- smooth.basis(argvals = new_ZT, y = Total_curves_ALL_rmout, fdParobj = fdobj)
#plot(sp_totalsmooth, lty = "solid")

ori_obj <- deriv.fd(sp_totalsmooth$fd,0)
vel_obj <- deriv.fd(sp_totalsmooth$fd,1)
acc_obj <- deriv.fd(sp_totalsmooth$fd,2)

#oeval 
allwx1wx28<- eval.fd(time,vel_obj)

allwx1wx28[, 1:42] <- allwx1wx28[, 1:42] + 0.5



# Assuming `df` is your original dataframe and `time` is your time vector
allwx1wx28 <- data.frame(allwx1wx28)
# Reshape data into long format
df_long <- allwx1wx28 %>%
  mutate(time = time) %>%  # Add the time column
  pivot_longer(cols = starts_with("WxT_"), names_to = "variable", values_to = "value")

# Add a new column `group` to indicate whether the variable belongs to "WxT_1." or "WxT_28."
df_long <- df_long %>%
  mutate(group = ifelse(str_detect(variable, "WxT_34\\."), "Genotype 1", "Genotype 2"))


svg("../../../My_Papers/photoperiod_shift/Figures/Spread.svg", width = 7.2, height = 6)

ggplot(df_long, aes(x = time, y = value, group = variable, color = group)) +
  geom_line(alpha = 1) +
  geom_vline(xintercept = 104, linetype = "dashed", color = "black", linewidth = 1) + # Set transparency
  scale_color_manual(values = c(
    "Group 1" = "#F8766D",  # Red
    "Group 2" = "#619CFF"   # Blue
  )) +
  labs(title = "Spread",
       x = "Time",
       y = "") +
  theme_minimal()  +
  theme(axis.text.x = element_text(size = 14, color = "black"),  # Increase x-axis number size
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.text.y = element_blank(),  # Remove y-axis numbers
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Increase gap to plot edges
        legend.position = "none",  # Adjust legend position (optional)
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  ) + annotate("text", x = 104, y = 0.2, label = "Missed light cue", 
               hjust = 1.1, vjust = -1, size = 5, color = "black")  # Axis lines style

dev.off()
#

######################################################################################


#look at shift

shift <- catagories[catagories$V2 == 'Phase_shift',]

# Extract column names of df1 that are in the first column of df2
filtered_columns <- colnames(functional_data) %in% shift$V1

# Subset df1 to keep only these columns
filtered_functional_Shape <- functional_data[, filtered_columns]

filtered_functional_Shape <- filtered_functional_Shape[,-c(20,21,81,82,83,84,85,86,87,88)]
# want to find the max and min of each column 
# Find the maximum of each column
max_values <- apply(filtered_functional_Shape, 2, max, na.rm = TRUE)

# Find the minimum of each column
min_values <- apply(filtered_functional_Shape, 2, min, na.rm = TRUE)

both <- cbind(max_values,min_values)

# Calculate the range for each column by subtracting min from max
rangee <- abs(both[,1] - both[,2])

# Find the column with the largest range
which.max(rangee)

#look at euclid_dist_pre_post_vel_tSDLD. Add row names from triat table to see what genotypes 
rownames(filtered_functional_Shape) <- genotypes

#genotypes wxT 42 and WxT 50 


# want to look at the pre and post for these two genotypes separatly, want to stack the pre and post 

# load in the all curves

meancurves <- read.csv("genotype_means_SDLD.csv")

#cols 34, 42 
#meancurves <- read.csv("meancurves_SDLD_vel.csv")
meancurves_pre1 <- meancurves[1:150,58]
meancurves_post1 <- meancurves[151:300,58]
meancurves_pre2 <- meancurves[1:150,39] + 1
meancurves_post2 <- meancurves[151:300,39] + 1

time <- seq(1,48, length.out = 150)

meancurves <- data.matrix(cbind(meancurves_pre1,meancurves_post1,meancurves_pre2,meancurves_post2))

BASIS_bspline <- create.bspline.basis(rangeval = c(1,48),norder = 4,nbasis = 302)
fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 0.1)
smooth_both_ori <- smooth.basis(argvals = time, y = meancurves, fdParobj = fdobj)
#plot(smooth_both_ori)

ori_obj <- deriv.fd(smooth_both_ori$fd,0)
vel_obj <- deriv.fd(smooth_both_ori$fd,1)
acc_obj <- deriv.fd(smooth_both_ori$fd,2)
#oeval 
meancurves <- eval.fd(time,ori_obj)
meancurves <- data.frame(meancurves)


# Select only the columns you need (WxT_15. and WxT_58.) and add the time column
df_selected <- data.frame(time = time, Genotype_1_pre = meancurves$meancurves_pre1, Genotype_1_post = meancurves$meancurves_post1, Genotype_2_pre = meancurves$meancurves_pre2, Genotype_2_post = meancurves$meancurves_post2)
# Add an offset to one of the columns to stack the curves
df_selected <- df_selected %>%
  mutate(WxT_42_pre_offset = WxT_42_pre + 0.0)  # Adjust the offset value as needed

# Convert to long format for ggplot
df_long <- df_selected %>%
  pivot_longer(cols = c("Genotype_1_pre", "Genotype_1_post","Genotype_2_pre", "Genotype_2_post"), names_to = "variable", values_to = "value")


svg("../../../My_Papers/photoperiod_shift/Figures/Shift.svg", width = 7.2, height = 6)

# to plot 
ggplot(df_long, aes(x = time, y = value, color = variable)) +
  geom_line(aes(linetype = variable), linewidth = 1) +
  scale_linetype_manual(values = c(
    "Genotype_1_pre" = "solid", 
    "Genotype_1_post" = "dashed", 
    "Genotype_2_pre" = "solid", 
    "Genotype_2_post" = "dashed"
  )) +  # Custom linetypes
  scale_color_manual(values = c(
    "Genotype_1_pre" = "#F8766D", 
    "Genotype_1_post" = "#F8766D", 
    "Genotype_2_pre" = "#619CFF", 
    "Genotype_2_post" = "#619CFF"
  )) +  # Custom colors
  labs(title = "Shift", x = "Time (h)", y = "") +
  theme_minimal()  +
  theme(axis.text.x = element_text(size = 14, color = "black"),  # Increase x-axis number size
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.text.y = element_blank(),  # Remove y-axis numbers
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Increase gap to plot edges
        legend.position = "none",  # Adjust legend position (optional)
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )

dev.off()


#heat map plot 

##################################

#Heatmap showing Spearman correlations between functional traits  developmental traits at flowering. Functional traits are categorized into five groups: Shape, Spread of Shape, Spread of Shift, Shift, and Phase Shift,
# looking at the correlation between functional traits and developmental traits 

# ready in functional traits 
setwd("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/")

# read in developmental QTLs linked 

developmental <- read.csv("QTLS_developmental_traits.csv")
# take just the name of the traits 
developmental <- developmental[1:66,1]

#read in coralations
correlations <- read.csv("correlation_traits.csv")
rownames(correlations) <- correlations[,1]
correlations <- correlations[,-1]

#filter based on developmental traits
correlations_dev <- correlations %>% dplyr::select(any_of(developmental))
#read in catagories from post QTL table
catagories <- read.csv("catagories_post_QTL.csv", header = TRUE)

#filter catagories based on orrelation names developmental 
df_filtered <- catagories %>%
  filter(.[[1]] %in% colnames(correlations_dev))

#check they are the same 
identical(df_filtered[,1], colnames(correlations_dev))
catagories <- data.frame(df_filtered[,-1])

# Rename the category column (if it has extra symbols or spaces)
colnames(catagories) <- "catagories"  # Ensure it matches your desired label
rownames(catagories) <- colnames(correlations_dev)


# for values that are shape do absolute correlation 
shape <- rownames(catagories)[which(catagories[[1]] == "Shape")]

correlations_dev_abs <- correlations_dev %>%
  mutate(across(any_of(shape), abs))

# Define the color palette for the annotations
annotation_colors <- list(
  catagories = c(Shape = "lightblue", Spread_of_Shape = "darkgreen", Spread_of_Shift = "green", Shift = "pink", Phase_Shift = "pink4")
)

# Generate the heatmap with annotations
myheatmap <- pheatmap(correlations_dev_abs,
                      legend_breaks = c(-0.5, -0.4, -0.2, 0, 0.2, 0.4, 0.5),
                      clustering_method = "complete", 
                      color = colorRampPalette(c("blue", "white", "red"))(50), 
                      show_rownames = TRUE, 
                      show_colnames = FALSE, 
                      treeheight_row = 30, 
                      cluster_rows = TRUE, 
                      cluster_cols = TRUE,
                      annotation_col = catagories,        # Use cleaned category annotations
                      annotation_colors = annotation_colors  # Specify the colors for each category
)
###########



