
##       PLOTS OF QTLS ===============

#read in saved output from QTL analysis. 
load("QTL_maping_lod_filtered_table.RData")


#   For chromosome 3 
# want to have plots of FT se and the shared traits
# these are in SNP regions 597 and 619 

# take the outhk and keep just the 597 and 619 rows to find the shared traits


filtered <- summary_lod_bypheno[c(88,90),-c(1,2)]

#to see where the QTL starts not just the max . #then need to make just the columns I wat 
which(colnames(out_hk) == "FT_se") 
which(colnames(out_hk) == "sedist_warp_PC2_velSDLD") 
which(colnames(out_hk) == "varmx_se_total_euclid_SDLD_vel") 
look <-  out_hk[390:414,c(1,2,987,852,217)]

#remove columns with sum zeor to find the triats 

filtered <- filtered[,colSums(filtered, na.rm = TRUE) != 0]

#traits sharing 

triats_sharing <- colnames(filtered)

# do the plot 

#first do plot using all the information 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")
# looking at QTL plots 
par(mar = c(5, 5, 4, 2))

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr3_ftse.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 4), lodcolumn = which(colnames(out_hk) == "FT_se") - 2,  col = "grey",  ylab = "LOD score",  main = "Standard error of flowering time", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "FT_se") - 2], col = "grey",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr3_spread_of_shift.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 4), lodcolumn = which(colnames(out_hk) == "sedist_warp_PC2_velSDLD") - 2,  col = "#E37383",  ylab = "LOD score",  main = "Spread of Shift", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "sedist_warp_PC2_velSDLD") - 2], col = "#E37383",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr3_spread_of_shape.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 4), lodcolumn = which(colnames(out_hk) == "varmx_se_total_euclid_SDLD_vel") - 2,  col = "#2171b5",  ylab = "LOD score",  main = "Spread of Shape", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "varmx_se_total_euclid_SDLD_vel") - 2], col = "#2171b5",lwd=3)

dev.off()


############### only chr3 
#b,l,t,r
par(mar=c(4.5,4.5,2,10))

#save using the Export option using SVG and 650x 450!!!!! 
#You can do a dark grey for the physiological trait, light blue for spread of shift, pink for spread of shape and dark red for shift, or something like that


#svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr3_grouped.svg", width = 7.2, height = 6)
#to save 
#plot info here 
plot(out_hk, ylim = c(0,4), lodcolumn = c(which(colnames(out_hk) == "FT_se") - 2,which(colnames(out_hk) == "sedist_warp_PC2_velSDLD") - 2,which(colnames(out_hk) == "varmx_se_total_euclid_SDLD_vel") - 2), col = c("grey44","#E37383","#2171b5"), chr = 3, ylab = "LOD score", cex.lab = 1.2,  cex.axis = 1.2,lwd =2, xlab = "Chromosome 3 position (cM)" )
abline(h=threshold_5_HK[which(colnames(out_hk) == "FT_se") - 2], col="grey44",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "sedist_warp_PC2_velSDLD") - 2], col="#E37383",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "varmx_se_total_euclid_SDLD_vel") - 2], col="#2171b5",lwd =2)
# add line for centromere 
rect(xleft = 39.480, xright = 40.611, ybottom = 0, ytop = 6, col = rgb(.128, .128, .128, alpha = 0.2), border = NA)

#add line for KH17 star the middle of the gene is 38.34955
abline(v=34.15209502, col="black", lty="dashed", lwd =2 )
#pink for spread of shape #E37383
#33.15209502
#add ledgend 
legend(x=92, y = 3.5,
       inset = c(-0.4, 0),  # Move the legend outside the plot area
       legend = c("FT se", "Spread of Shift", "Spread of Shape", "Centromere","AtKH17"), 
       col = c("grey44", "#E37383", "#2171b5", rgb(0.128, 0.128, 0.128, alpha = 0.2),"black"), 
       lty = c(1, 1, 1, NA,2), 
       lwd = c(3, 3, 3, NA,3), 
       pch = c(NA, NA, NA, 15,NA), 
       pt.cex = 2, 
       bty = "n", 
       xpd = TRUE,  # Allow drawing outside the plot region
       y.intersp = 1.5)

#dev.off()

#box plots! 

# read in the genotype data  

genotype_data <- read.csv("QTL_plots_genotype_info.csv")

#make sure columns are as factors 
genotype_data$SNP597 <- as.factor(genotype_data$SNP597)
genotype_data$SNP619 <- as.factor(genotype_data$SNP619)

#make plots for both the QTL regions 
#for FT_se there will be two as it has a QTL in both the regions 

genotype_data_597 <- genotype_data[genotype_data[[8]] %in% c("W", "X"), ]
genotype_data_619 <- genotype_data[genotype_data[[9]] %in% c("W", "X"), ]
###
#to work out the significance levels of each!

#first do normaility test. If the p-value of the test is greater than α = .05, then the data is assumed to be normally distributed. Also check histagram and qqnorm plots

#FT_se 
#shapiro.test(genotype_data_597$FT_se)$p.value > 0.05
#hist(genotype_data_597$FT_se)
#qqnorm(genotype_data_597$FT_se)
#qqline(genotype_data_597$FT_se)

#shapiro.test(genotype_data_619$FT_se)$p.value > 0.05
#hist(genotype_data_619$FT_se)
#qqnorm(genotype_data_619$FT_se)
#qqline(genotype_data_619$FT_se)
#
#shapiro.test(genotype_data_597$sedist_warp_PC2_velSDLD)$p.value > 0.05
#hist(genotype_data_597$sedist_warp_PC2_velSDLD)
#qqnorm(genotype_data_597$sedist_warp_PC2_velSDLD)
#qqline(genotype_data_597$sedist_warp_PC2_velSDLD)
#
#shapiro.test(genotype_data_619$varmx_se_total_euclid_SDLD_vel)$p.value > 0.05
#hist(genotype_data_619$varmx_se_total_euclid_SDLD_vel)
#qqnorm(genotype_data_619$varmx_se_total_euclid_SDLD_vel)
#qqline(genotype_data_619$varmx_se_total_euclid_SDLD_vel)


#now do a t-test for each of the boxplots.... 

library(ggpubr)
t_test_result <- t.test(FT_se ~ SNP597, data = genotype_data_597)
print(t_test_result)


# If a p-value is less than 0.05, it is flagged with one star (*). If a p-value is less than 0.01, it is flagged with 2 stars (**). If a p-value is less than 0.001, it is flagged with three stars (***)

plot1 <-ggplot(genotype_data_597, aes(x = SNP597, y = FT_se)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 0.62, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP597", x = "Genotype", y = "Standard error of Flowering Time") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )


plot2 <- ggplot(genotype_data_619, aes(x = SNP619, y = FT_se)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 0.62, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP619", x = "Genotype", y = "Standard error of Flowering Time") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )

(plot1|plot2)

#save it 
svg("../../../My_Papers/photoperiod_shift/Figures/Boxplot_chr3.svg", width = 7.2, height = 6)
(plot1|plot2)

dev.off()


t_test_result <- t.test(FT_se ~ SNP619, data = genotype_data_619)
print(t_test_result)
#
#for this figure also want to look at the warping functions that are in the two groups in SNP 619 block 

#look at the split 
# Assuming the column with factor levels is called 'SNP619'
grouped_rows <- split(genotype_data_619$Genotype, genotype_data_619$SNP619)

# Print row names for each group
grouped_rows

#remove wxT 31 as not in final map removeed... 

grouped_rows$W <- grouped_rows$W[-21] 

# load in the warping functions 

#origional curves
load(file = "../../../Warping/56_152/full_length_curves/origional/reg_fulllength_to_means_oriSDLD.RData")

reg_oriSDLD <- regSDLD


#velocity curves
load(file = "56_152/full_length_curves/velocity/reg_fulllength_to_means_velSDLD.RData")
reg_velSDLD <- regSDLD


#evalulate the curves so can plot them?

#first create mean for each 
means<- lapply(X = seq(1,64), FUN = function(x){mean.fd(reg_oriSDLD[[x]]$warpfd)})

#eval 

new_time <- seq(56,152, length.out = 300)

evlameans <- sapply(X = seq(1,64), FUN = function(x){eval.fd(evalarg = new_time,means[[x]],0)})

Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names <- sort(Genotype_names)
colnames(evlameans) <- Genotype_names

#group for ws 

WS_group <- evlameans[,grouped_rows$W]
WS_group <- WS_group - new_time
WS_group <- data.frame(new_time,WS_group)

TNZ_group <- evlameans[,grouped_rows$X]
TNZ_group <- TNZ_group - new_time
#TNZ_group <- TNZ_group +3 # buffer
TNZ_group <- data.frame(new_time,TNZ_group)

# Add an identifier column and combine the dataframes
df1_long <- WS_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  mutate(Source = "DF1")

df2_long <- TNZ_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  mutate(Source = "DF2")

combined_df <- bind_rows(df1_long, df2_long)

# Plot with ggplot
ggplot(combined_df, aes(x = new_time, y = Value, group = interaction(Variable, Source), color = Source)) +
  geom_line() +
  scale_color_manual(values = c("DF1" = "blue", "DF2" = "red")) +  # Custom colors
  labs(title = "Combined Dataframes Plot",
       x = "Time",
       y = "Value",
       color = "Source") +
  theme_minimal()
#


# Calculate the mean and standard deviation for each group at each time point
WS_group_summary <- WS_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

TNZ_group_summary <- TNZ_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

# Combine summaries with a Source identifier
combined_summary <- bind_rows(
  WS_group_summary %>% mutate(Source = "Ws-2"),
  TNZ_group_summary %>% mutate(Source = "Tnz-1")
)

# Plot with ggplot
ggplot(combined_summary, aes(x = new_time, y = Mean, color = Source, fill = Source)) +
  geom_line(size = 1.2) +  # Plot mean lines
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2, color = NA) +  # Add shadow for noise
  geom_vline(xintercept = 104, linetype = "dashed", color = "black", size = 1) +  # Add dashed vertical line
  scale_color_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom line colors
  scale_fill_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom shadow colors
  labs(
    title = "",
    x = "Time (h)",
    y = "",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14,, color = "black"),  # Increase x-axis number size
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.text.y = element_blank(),  # Remove y-axis numbers
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Increase gap to plot edges
    legend.position = "right",  # Adjust legend position (optional)
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  ) + annotate("text", x = 104, y = 1.2, label = "Missed light cue", 
               hjust = 1.1, vjust = -1, size = 5, color = "black")

##################################################################################################################



#############################     CHR 5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




########################################################################################



##       PLOTS OF QTLS


#read in saved output 
load("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/QTL_maping_lod_filtered_table.RData")


#   For chromosome 5
# want to have plots of meanFT and the shared traits so make outhk 0 and 1s to find the SNP regions 

looking <- out_hk

for (i in 3:ncol(looking)) {
  looking[,i] <- ifelse(looking[,i] > threshold_5_HK[(i-2)], 1, 0)
}

triat <- looking[,-c(1,2)]
triat<- triat[,colSums(triat, na.rm = TRUE) != 0]
#and remove rows that contain zero too
triat <- triat[rowSums(triat, na.rm = TRUE) != 0 ,]

#look for rows where meanFT has a one 
traits_share_meanFT <- triat[triat[,273] != 0 ,]

#remove any rows or columns are zero
traits_share_meanFT<- traits_share_meanFT[,colSums(traits_share_meanFT, na.rm = TRUE) != 0]
#and remove rows that contain zero too
traits_share_meanFT <- traits_share_meanFT[rowSums(traits_share_meanFT, na.rm = TRUE) != 0 ,]

#look at which traits this is only with meanFT
traits_share_meanFT <- traits_share_meanFT[,1:7]

#remove any row that doesnt = 0 
#and remove rows that contain zero too
traits_share_meanFT <- traits_share_meanFT[rowSums(traits_share_meanFT[,1:6], na.rm = TRUE) != 0 ,]

# so there are 4 regions that qualifiy.... 262-264,1098-1099, 1100, 1104-1109.

# and 6 traits... 

] "PC2_gc_SDSD_acc"           "varmx_PC2_gc_SDSD_acc"     "PC2_gc_acc_SDSD_pre"       "PC2_gc_acc_SDSD_post"     
[5] "ratio_SDLD_acc_centerdist" "sedist_warp_PC1_velLDLD" 


#to see where the QTL starts not just the max . #then need to make just the columns I wat 
which(colnames(out_hk) == "FT_mean") 
which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") # this is shift 
which(colnames(out_hk) == "sedist_warp_PC1_velLDLD") #spread of shift 
look <-  out_hk[c(209,210,211,776,777,778),c(1,2,986,518,857)]

#first do plot using all the information 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")
# looking at QTL plots 
par(mar = c(5, 5, 4, 2))

#plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "PC2_gc_acc_SDSD_pre") - 2,  col = "red",  ylab = "LOD score",  main = "Standard error of flowering time", cex.lab = 1.2,  cex.axis = 1.2)
#abline(h = threshold_5_HK[which(colnames(out_hk) == "PC2_gc_acc_SDSD_pre") - 2], col = "red")


par(mar = c(5, 5, 4, 2))

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr5_meanft.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "FT_mean") - 2,  col = "grey44",  ylab = "LOD score",  main = "Mean flowering time", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "FT_mean") - 2], col = "grey44",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr5_shift.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2,  col = "darkred",  ylab = "LOD score",  main = "Shift", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2], col = "darkred",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr5_spread_of_shift.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "sedist_warp_PC1_velLDLD") - 2,  col = "#2171b5",  ylab = "LOD score",  main = "Spread of Shift", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "sedist_warp_PC1_velLDLD") - 2], col = "#2171b5",lwd=3)

dev.off()


#
#You can do a dark grey for the physiological trait, light blue for spread of shift, pink for spread of shape and dark red for shift, or something like that
############### only chr5 and 1 
#b,l,t,r

#green #2171b5 darkred

par(mar=c(4.5,4.5,2,10))

plot(out_hk, ylim = c(0,7.5), lodcolumn = c(which(colnames(out_hk) == "FT_mean") - 2,which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2,which(colnames(out_hk) == "sedist_warp_PC1_velLDLD") - 2), col = c("grey44","darkred","#E37383"), chr = c(1,5), ylab = "LOD score", cex.lab = 1.2,  cex.axis = 1.2,lwd =2 )
abline(h=threshold_5_HK[which(colnames(out_hk) == "FT_mean") - 2], col="grey44",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2], col="darkred",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "sedist_warp_PC1_velLDLD") - 2], col="#E37383",lwd =2)
#add maf2 line 
abline(v=(122.7437836+25+129.7021294), col="black", lty="dashed", lwd =2)

#for kh29 122.7670540134
abline(v=(122.7437836+25+122.7670540134), col="black", lty=3, lwd =2)

#add ledgend 
legend(x=300, y = 4.5,
       inset = c(-0.4, 0),  # Move the legend outside the plot area
       legend = c("FT mean", "Shift", "Spread of Shift", "MAF2","KH29"), 
       col = c("grey44", "darkred", "#E37383","black","black"), 
       lty = c(1, 1, 1,2,3), 
       lwd = c(3, 3, 3,3,3), 
       pt.cex = 2, 
       bty = "n", 
       xpd = TRUE,  # Allow drawing outside the plot region
       y.intersp = 1.5)


# try doing box plots! 

# read in the genotype data  

genotype_data <- read.csv("QTL_plots_genotype_info.csv")

#make sure columns are as factors 

genotype_data$SNP262 <- as.factor(genotype_data$SNP262)
genotype_data$SNP1104 <- as.factor(genotype_data$SNP1104)

#make plots for both the QTL regions 

#for FT_se there will be two as it has a QTL in both the regions 

genotype_data_262 <- genotype_data[genotype_data[[10]] %in% c("W", "X"), ]
genotype_data_1104 <- genotype_data[genotype_data[[11]] %in% c("W", "X"), ]
###
#to work out the significance levels of each!

#first do normaility test. If the p-value of the test is greater than α = .05, then the data is assumed to be normally distributed. Also check histagram and qqnorm plots

#FT_mean
shapiro.test(genotype_data_262$FT_mean)$p.value > 0.05
hist(genotype_data_262$FT_mean)
qqnorm(genotype_data_262$FT_mean)
qqline(genotype_data_262$FT_mean)

shapiro.test(genotype_data_1104$FT_mean)$p.value > 0.05
hist(genotype_data_1104$FT_mean)
qqnorm(genotype_data_1104$FT_mean)
qqline(genotype_data_1104$FT_mean)
#
shapiro.test(genotype_data_262$sedist_warp_PC1_velLDLD)$p.value > 0.05
hist(genotype_data_262$sedist_warp_PC1_velLDLD)
qqnorm(genotype_data_262$sedist_warp_PC1_velLDLD)
qqline(genotype_data_262$sedist_warp_PC1_velLDLD)
#
shapiro.test(genotype_data_1104$ratio_SDLD_acc_centerdist)$p.value > 0.05
hist(genotype_data_1104$ratio_SDLD_acc_centerdist)
qqnorm(genotype_data_1104$ratio_SDLD_acc_centerdist)
qqline(genotype_data_1104$ratio_SDLD_acc_centerdist)


#now do a t-test for each of the boxplots.... 

library(ggpubr)
t_test_result <- t.test(FT_mean ~ SNP262, data = genotype_data_262)
print(t_test_result)
t_test_result <- t.test(FT_mean ~ SNP1104, data = genotype_data_1104)
print(t_test_result)


# If a p-value is less than 0.05, it is flagged with one star (*). If a p-value is less than 0.01, it is flagged with 2 stars (**). If a p-value is less than 0.001, it is flagged with three stars (***)

plot1 <-ggplot(genotype_data_262, aes(x = SNP262, y = FT_mean)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 43, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP262", x = "Genotype", y = "Mean Flowering Time") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )


plot2 <- ggplot(genotype_data_1104, aes(x = SNP1104, y = FT_mean)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 43, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP1104", x = "Genotype", y = "Mean Flowering Time") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )

(plot1|plot2)
t_test_result <- t.test(FT_mean ~ SNP1104, data = genotype_data_1104)
print(t_test_result)

#plot
svg("../../../My_Papers/photoperiod_shift/Figures/Boxplot_chr5.svg", width = 7.2, height = 6)
(plot1|plot2)

dev.off()


#
#for this figure also want to look at the warping functions that are in the two groups in SNP 1104 block 

#look at the split 
# Assuming the column with factor levels is called 'SNP1104'
grouped_rows <- split(genotype_data_1104$Genotype, genotype_data_1104$SNP1104)

# Print row names for each group
grouped_rows

#remove wxT 31 as not in final map removeed... 

#grouped_rows$W <- grouped_rows$W[-21] 


# in this case looking at the shift between ws and tnz so want to plot the first half vs the second half for genotype s in WS group and genotypes in tnz group....



# could either look at the mean warping function between the ws and tnz 



meancurves <- read.csv("genotype_means_SDLD.csv")
##################



## do plot for spread in CHR5 ........................... 


# read in mean curves the trait was sedist_warp_PC1_velLDLD so use LDLD curves 
#origional curves
load(file = "../../../Warping/56_152/full_length_curves/origional/reg_fulllength_to_means_oriLDLD.RData")

reg_oriLDLD <- regLDLD


#first create mean for each 
means<- lapply(X = seq(1,64), FUN = function(x){mean.fd(reg_oriLDLD[[x]]$warpfd)})

#eval 
new_time <- seq(56,152, length.out = 300)

evlameans <- sapply(X = seq(1,64), FUN = function(x){eval.fd(evalarg = new_time,means[[x]],0)})

Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names <- sort(Genotype_names)
colnames(evlameans) <- Genotype_names

#group for ws 

WS_group <- evlameans[,grouped_rows$W]
WS_group <- WS_group - new_time
WS_group <- data.frame(new_time,WS_group)

TNZ_group <- evlameans[,grouped_rows$X]
TNZ_group <- TNZ_group - new_time
#TNZ_group <- TNZ_group +3 # buffer
TNZ_group <- data.frame(new_time,TNZ_group)

# Calculate the mean and standard deviation for each group at each time point
WS_group_summary <- WS_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

TNZ_group_summary <- TNZ_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

# Combine summaries with a Source identifier
combined_summary <- bind_rows(
  WS_group_summary %>% mutate(Source = "Ws-2"),
  TNZ_group_summary %>% mutate(Source = "Tnz-1")
)

# Plot with ggplot
ggplot(combined_summary, aes(x = new_time, y = Mean, color = Source, fill = Source)) +
  geom_line(size = 1.2) +  # Plot mean lines
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2, color = NA) +  # Add shadow for noise
  geom_vline(xintercept = 104, linetype = "dashed", color = "black", size = 1) +  # Add dashed vertical line
  scale_color_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom line colors
  scale_fill_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom shadow colors
  labs(
    title = "",
    x = "Time (h)",
    y = "",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14,, color = "black"),  # Increase x-axis number size
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.text.y = element_blank(),  # Remove y-axis numbers
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Increase gap to plot edges
    legend.position = "right",  # Adjust legend position (optional)
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  ) + annotate("text", x = 104, y = 1.2, label = "Missed light cue", 
               hjust = 1.1, vjust = -1, size = 5, color = "black")






##################################################################################################################
















#

meancurves <- read.csv("genotype_means_SDLD.csv")

#cols 34, 42 
#meancurves <- read.csv("meancurves_SDLD_vel.csv")
meancurves_pre <- meancurves[1:150,-1]
meancurves_post <- meancurves[151:300,-1]
time <- seq(1,48, length.out = 150)

WS_group_pre <- meancurves_pre[1:150,grouped_rows$W]
WS_group_post <- meancurves_post[1:150,grouped_rows$W]
WS_group_premean <- apply(WS_group_pre,1,mean)
WS_group_postmean <- apply(WS_group_post,1,mean)
WS_group_presd <- apply(WS_group_pre,1,sd)
WS_group_postsd <- apply(WS_group_post,1,sd)

TNZ_group_pre <- meancurves_pre[1:150,grouped_rows$X] 
TNZ_group_post <- meancurves_post[1:150,grouped_rows$X]
TNZ_group_premean <- apply(TNZ_group_pre,1,mean)
TNZ_group_postmean <- apply(TNZ_group_post,1,mean)
TNZ_group_presd <- apply(TNZ_group_pre,1,sd)
TNZ_group_postsd <- apply(TNZ_group_post,1,sd)

# nust do the 4 sets of curve anad calculate the mean and sd of these for the ribbion 

# Select only the columns you need (WxT_15. and WxT_58.) and add the time column
df_selected <- data.frame(time = time, WS_pre = WS_group_premean, WS_post =WS_group_postmean, TNZ_pre = TNZ_group_premean, TNZ_post = TNZ_group_postmean)
# Add an offset to one of the columns to stack the curves
# Convert to long format for ggplot
df_long <- df_selected %>%
  pivot_longer(cols = c("WS_pre", "WS_post","TNZ_pre", "TNZ_post"), names_to = "variable", values_to = "value")

# Combine data for mean and SD into long format before plotting
df_long <- df_selected %>%
  pivot_longer(
    cols = c("WS_pre", "WS_post", "TNZ_pre", "TNZ_post"),
    names_to = "variable", 
    values_to = "value"
  ) %>%
  pivot_longer(
    cols = c("WS_pre_sd", "WS_post_sd", "TNZ_pre_sd", "TNZ_post_sd"),
    names_to = "variable_sd", 
    values_to = "SD"
  ) %>%
  filter(
    (variable == "WS_pre" & variable_sd == "WS_pre_sd") |
      (variable == "WS_post" & variable_sd == "WS_post_sd") |
      (variable == "TNZ_pre" & variable_sd == "TNZ_pre_sd") |
      (variable == "TNZ_post" & variable_sd == "TNZ_post_sd")
  )  # Match SD values to the correct group

# Plot with ribbons
ggplot(df_long, aes(x = time, y = value, color = variable)) +
  geom_line(aes(linetype = variable), linewidth = 1) +
  geom_ribbon(
    aes(ymin = value - SD, ymax = value + SD, fill = variable), 
    alpha = 0.2, 
    color = NA
  ) +  # Add ribbons for mean ± SD
  scale_linetype_manual(values = c(
    "WS_pre" = "solid", 
    "WS_post" = "dashed", 
    "TNZ_pre" = "solid", 
    "TNZ_post" = "dashed"
  )) +  # Custom linetypes
  scale_color_manual(values = c(
    "WS_pre" = "#F8766D", 
    "WS_post" = "#F8766D", 
    "TNZ_pre" = "#619CFF", 
    "TNZ_post" = "#619CFF"
  )) +  # Custom colors
  scale_fill_manual(values = c(
    "WS_pre" = "#F8766D", 
    "WS_post" = "#F8766D", 
    "TNZ_pre" = "#619CFF", 
    "TNZ_post" = "#619CFF"
  )) +  # Custom fill colors
  labs(title = "Shift", x = "Time (h)", y = "", fill = "Group") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )







#############################   LEAVES!!!!!!!!!





########################################################################################



##       PLOTS OF QTLS


#read in saved output 
load("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/QTL_maping_lod_filtered_table.RData")


#   For chromosome 5
# want to have plots of meanLeaves and the shared traits so make outhk 0 and 1s to find the SNP regions 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")

looking <- out_hk

for (i in 3:ncol(looking)) {
  looking[,i] <- ifelse(looking[,i] > threshold_5_HK[(i-2)], 1, 0)
}

triat <- looking[,-c(1,2)]
triat<- triat[,colSums(triat, na.rm = TRUE) != 0]
#and remove rows that contain zero too
triat <- triat[rowSums(triat, na.rm = TRUE) != 0 ,]

#look for rows where mean leaves has a one 
traits_share_meanLeaves <- triat[triat[,275] != 0 ,]

#remove any rows or columns are zero
traits_share_meanLeaves<- traits_share_meanLeaves[,colSums(traits_share_meanLeaves, na.rm = TRUE) != 0]
#and remove rows that contain zero too
traits_share_meanLeaves <- traits_share_meanLeaves[rowSums(traits_share_meanLeaves, na.rm = TRUE) != 0 ,]

#remove any row that doesnt = 0 
#and remove rows that contain zero too wher trait is non-developmental
traits_share_meanLeaves <- traits_share_meanLeaves[rowSums(traits_share_meanLeaves[,1:8], na.rm = TRUE) != 0 ,]

# so there are 5 regions that qualifiy.... 262-264,281-289,1098-1099, 1100, 1104-1109.

# and 8 traits... 
[1] "PC2_gc_SDSD_acc"           "varmx_PC2_gc_SDSD_acc"     "PC2_gc_vel_SDSD_post"      "PC2_gc_acc_SDSD_pre"      
[5] "PC2_gc_acc_SDSD_post"      "ratio_SDLD_acc_centerdist" "RSQR_velSDSD"              "sedist_warp_PC1_velLDLD" 



#to see where the QTL starts not just the max . #then need to make just the columns I wat 
which(colnames(out_hk) == "leaves_mean") 
which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") # this is shift 
which(colnames(out_hk) == "PC2_gc_vel_SDSD_post") #this is shape 

#first do plot using all the information 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")
# looking at QTL plots 
par(mar = c(5, 5, 4, 2))

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr5_meanLeaves.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "leaves_mean") - 2,  col = "grey44",  ylab = "LOD score",  main = "Mean number of leaves", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "leaves_mean") - 2], col = "grey44",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_leaves_shift.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2,  col = "darkred",  ylab = "LOD score",  main = "Shift", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2], col = "darkred",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_leaves_shape.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "PC2_gc_vel_SDSD_post") - 2,  col = "#6baed6",  ylab = "LOD score",  main = "Shape", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "PC2_gc_vel_SDSD_post") - 2], col = "#6baed6",lwd=3)

dev.off()


#
#You can do a dark grey for the physiological trait, light blue for spread of shift, pink for spread of shape and dark red for shift, or something like that
############### only chr5 and 1 
#b,l,t,r

#green #2171b5 darkred

par(mar=c(4.5,4.5,2,10))

plot(out_hk, ylim = c(0,9), lodcolumn = c(which(colnames(out_hk) == "leaves_mean") - 2,which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2,which(colnames(out_hk) == "PC2_gc_vel_SDSD_post") - 2), col = c("grey44","darkred","#6baed6"), chr = c(1,5), ylab = "LOD score", cex.lab = 1.2,  cex.axis = 1.2,lwd =2 )
abline(h=threshold_5_HK[which(colnames(out_hk) == "leaves_mean") - 2], col="grey44",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2], col="darkred",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "PC2_gc_vel_SDSD_post") - 2], col="#6baed6",lwd =2)
#add maf2 line 
abline(v=(122.7437836+25+129.7021294), col="black", lty="dashed", lwd =2)

#for FT 122.7670540134
abline(v=(92.28854), col="black", lty=3, lwd =2)

#add ledgend 
legend(x=300, y = 4.5,
       inset = c(-0.4, 0),  # Move the legend outside the plot area
       legend = c("FT mean", "Shift", "Shape", "MAF2","FT"), 
       col = c("grey44", "darkred", "#6baed6","black","black"), 
       lty = c(1, 1, 1,2,3), 
       lwd = c(3, 3, 3,3,3), 
       pt.cex = 2, 
       bty = "n", 
       xpd = TRUE,  # Allow drawing outside the plot region
       y.intersp = 1.5)

#saved at "QTL_leaves_grouped"

# try doing box plots! 

# read in the genotype data  

genotype_data <- read.csv("QTL_plots_genotype_info.csv")

#make sure columns are as factors 

genotype_data$SNP281 <- as.factor(genotype_data$SNP281)
genotype_data$SNP1104 <- as.factor(genotype_data$SNP1104)

#make plots for both the QTL regions 

#for FT_se there will be two as it has a QTL in both the regions 

genotype_data_281 <- genotype_data[genotype_data[[12]] %in% c("W", "X"), ]
genotype_data_1104 <- genotype_data[genotype_data[[11]] %in% c("W", "X"), ]
###
#to work out the significance levels of each!

#first do normaility test. If the p-value of the test is greater than α = .05, then the data is assumed to be normally distributed. Also check histagram and qqnorm plots

#leaves_mean
shapiro.test(genotype_data_281$leaves_mean)$p.value > 0.05
hist(genotype_data_281$leaves_mean)
qqnorm(genotype_data_281$leaves_mean)
qqline(genotype_data_281$leaves_mean)

shapiro.test(genotype_data_1104$leaves_mean)$p.value > 0.05
hist(genotype_data_1104$leaves_mean)
qqnorm(genotype_data_1104$leaves_mean)
qqline(genotype_data_1104$leaves_mean)
#
shapiro.test(genotype_data_281$PC2_gc_vel_SDSD_post)$p.value > 0.05
hist(genotype_data_281$PC2_gc_vel_SDSD_post)
qqnorm(genotype_data_281$PC2_gc_vel_SDSD_post)
qqline(genotype_data_281$PC2_gc_vel_SDSD_post)
#
shapiro.test(genotype_data_1104$ratio_SDLD_acc_centerdist)$p.value > 0.05
hist(genotype_data_1104$ratio_SDLD_acc_centerdist)
qqnorm(genotype_data_1104$ratio_SDLD_acc_centerdist)
qqline(genotype_data_1104$ratio_SDLD_acc_centerdist)


#now do a t-test for each of the boxplots.... 

library(ggpubr)
t_test_result <- t.test(leaves_mean ~ SNP281, data = genotype_data_281)
print(t_test_result)
t_test_result <- t.test(leaves_mean ~ SNP1104, data = genotype_data_1104)
print(t_test_result)


# If a p-value is less than 0.05, it is flagged with one star (*). If a p-value is less than 0.01, it is flagged with 2 stars (**). If a p-value is less than 0.001, it is flagged with three stars (***)

plot1 <-ggplot(genotype_data_281, aes(x = SNP281, y = leaves_mean)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 32, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP281", x = "Genotype", y = "Mean Number of Leaves") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )


plot2 <- ggplot(genotype_data_1104, aes(x = SNP1104, y = leaves_mean)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 32, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP1104", x = "Genotype", y = "Mean Number of Leaves") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )

(plot1|plot2)
t_test_result <- t.test(FT_mean ~ SNP1104, data = genotype_data_1104)
print(t_test_result)

#plot
svg("../../../My_Papers/photoperiod_shift/Figures/Boxplot_mean_leaves.svg", width = 7.2, height = 6)
(plot1|plot2)

dev.off()


###########################################################################################




############################################################################################




########################################################################################



##       PLOTS OF QTLS


#read in saved output 
load("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/QTL_maping_lod_filtered_table.RData")


#   For chromosome 5
# want to have plots of seLeaves and the shared traits so make outhk 0 and 1s to find the SNP regions 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")

looking <- out_hk

for (i in 3:ncol(looking)) {
  looking[,i] <- ifelse(looking[,i] > threshold_5_HK[(i-2)], 1, 0)
}

triat <- looking[,-c(1,2)]
triat<- triat[,colSums(triat, na.rm = TRUE) != 0]
#and remove rows that contain zero too
triat <- triat[rowSums(triat, na.rm = TRUE) != 0 ,]

#look for rows where mean leaves has a one 
traits_share_seLeaves <- triat[triat[,276] != 0 ,]

#remove any rows or columns are zero
traits_share_seLeaves<- traits_share_seLeaves[,colSums(traits_share_seLeaves, na.rm = TRUE) != 0]
#and remove rows that contain zero too
traits_share_seLeaves <- traits_share_seLeaves[rowSums(traits_share_seLeaves, na.rm = TRUE) != 0 ,]

#remove any row that doesnt = 0 
#and remove rows that contain zero too wher trait is non-developmental
traits_share_seLeaves <- traits_share_seLeaves[rowSums(traits_share_seLeaves[,1:31], na.rm = TRUE) != 0 ,]

# so there are 4 regions that qualifiy.... 38, 39-42,736-739, 1104-1109.

# and 31 traits... 
[1] "PC3_gc_SDLD"                  "PC3_gc_LDLD"                  "PC2_gc_SDLD_vel"             
[4] "PC1_gc_LDLD_vel"              "PC2_gc_LDLD_vel"              "PC2_gc_SDLD_acc"             
[7] "PC2_gc_LDLD_vel"              "PC1_gc_LDLD_acc"              "varmx_PC2_gc_SDLD"           
[10] "varmx_PC2_gc_SDLD_vel"        "varmx_PC2_gc_LDLD_vel"        "varmx_PC3_gc_LDLD_vel"       
[13] "varmx_PC1_gc_SDLD_acc"        "varmx_PC2_gc_LDLD_vel"        "varmx_PC2_gc_LDLD_acc"       
[16] "varmx_PC3_gc_LDLD_acc"        "varmx_SDLD_PC2_center"        "varmx_SDLD_PC2_center_vel"   
[19] "varmx_PC3_gc_LDLD_acc" "PC3_gc_SDLD_pre"              "PC3_gc_SDLD_post"            
[22] "PC2_gc_vel_SDLD_pre"          "PC1_gc_vel_LDLD_pre"          "PC2_gc_vel_LDLD_pre"         
[25] "PC2_gc_vel_SDLD_post"         "PC2_gc_vel_LDLD_post"         "PC2_gc_acc_SDLD_pre"         
[28] "PC1_gc_acc_LDLD_pre"          "PC2_gc_acc_SDLD_post"         "PC1_gc_acc_LDLD_post"        
[31] "ratio_SDLD_acc_centerdist"    "FT_mean"                      "leaves_se"                 
[34] "leaves_se"                    "weights_mean"                 "weights_se"    


#   choose a trait that is on chr 1 and  chr 4 ideally one trait only on each so plot all 


#PC2_gc_LDLD_vel
#varmx_PC3_gc_LDLD_acc

#to see where the QTL starts not just the max . #then need to make just the columns I wat 
which(colnames(out_hk) == "leaves_se") 
which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") # this is shift 
which(colnames(out_hk) == "varmx_PC3_gc_LDLD_acc") #this is shape 

#first do plot using all the information 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")

ddsd <- out_hk[,colnames(traits_share_seLeaves)]
# looking at QTL plots 
par(mar = c(5, 5, 4, 2))

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr5_seLeaves.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "leaves_se") - 2,  col = "grey44",  ylab = "LOD score",  main = "Mean number of leaves", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "leaves_se") - 2], col = "grey44",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_leaves_shift.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2,  col = "darkred",  ylab = "LOD score",  main = "Shift", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2], col = "darkred",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_leaves_shape.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "varmx_PC2_gc_LDLD_vel") - 2,  col = "#6baed6",  ylab = "LOD score",  main = "Shape", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "varmx_PC2_gc_LDLD_vel") - 2], col = "#6baed6",lwd=3)

dev.off()


#
#You can do a dark grey for the physiological trait, light blue for spread of shift, pink for spread of shape and dark red for shift, or something like that
############### only chr5 and 1 
#b,l,t,r
ratio_SDLD_acc_centerdist
varmx_PC2_gc_LDLD_vel
#green #2171b5 darkred

par(mar=c(4.5,4.5,2,10))

plot(out_hk, ylim = c(0,9), lodcolumn = c(which(colnames(out_hk) == "leaves_se") - 2,which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2,which(colnames(out_hk) == "varmx_PC2_gc_LDLD_vel") - 2), col = c("grey44","darkred","#6baed6"), chr = c(1,5), ylab = "LOD score", cex.lab = 1.2,  cex.axis = 1.2,lwd =2 )
abline(h=threshold_5_HK[which(colnames(out_hk) == "leaves_se") - 2], col="grey44",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "ratio_SDLD_acc_centerdist") - 2], col="darkred",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "varmx_PC2_gc_LDLD_vel") - 2], col="#6baed6",lwd =2)
#add maf2 line 
abline(v=(122.7437836+25+129.7021294), col="black", lty="dashed", lwd =2)

#for FT 122.7670540134
#abline(v=(92.28854), col="black", lty=3, lwd =2)

#add ledgend 
legend(x=300, y = 4.5,
       inset = c(-0.4, 0),  # Move the legend outside the plot area
       legend = c("FT se", "Shift", "Shape", "MAF2"), 
       col = c("grey44", "darkred", "#6baed6","black"), 
       lty = c(1, 1, 1,2), 
       lwd = c(3, 3, 3,3), 
       pt.cex = 2, 
       bty = "n", 
       xpd = TRUE,  # Allow drawing outside the plot region
       y.intersp = 1.5)

#saved at "QTL_leaves_grouped"

# try doing box plots! 

# read in the genotype data  

genotype_data <- read.csv("QTL_plots_genotype_info.csv")

#make sure columns are as factors 

genotype_data$SNP42 <- as.factor(genotype_data$SNP42)
genotype_data$SNP1104 <- as.factor(genotype_data$SNP1104)

#make plots for both the QTL regions 

#for FT_se there will be two as it has a QTL in both the regions 

genotype_data_42 <- genotype_data[genotype_data[[13]] %in% c("W", "X"), ]
genotype_data_1104 <- genotype_data[genotype_data[[11]] %in% c("W", "X"), ]
###
#to work out the significance levels of each!

#first do normaility test. If the p-value of the test is greater than α = .05, then the data is assumed to be normally distributed. Also check histagram and qqnorm plots

#leaves_se
shapiro.test(genotype_data_42$leaves_se)$p.value > 0.05
hist(genotype_data_42$leaves_se)
qqnorm(genotype_data_42$leaves_se)
qqline(genotype_data_42$leaves_se)

shapiro.test(genotype_data_1104$leaves_se)$p.value > 0.05
hist(genotype_data_1104$leaves_se)
qqnorm(genotype_data_1104$leaves_se)
qqline(genotype_data_1104$leaves_se)
#



shapiro.test(genotype_data_42$ratio_SDLD_acc_centerdist)$p.value > 0.05
hist(genotype_data_42$ratio_SDLD_acc_centerdist)
qqnorm(genotype_data_42$ratio_SDLD_acc_centerdist)
qqline(genotype_data_42$ratio_SDLD_acc_centerdist)
#
shapiro.test(genotype_data_1104$varmx_PC2_gc_LDLD_vel)$p.value > 0.05
hist(genotype_data_1104$varmx_PC2_gc_LDLD_vel)
qqnorm(genotype_data_1104$varmx_PC2_gc_LDLD_vel)
qqline(genotype_data_1104$varmx_PC2_gc_LDLD_vel)


#now do a t-test for each of the boxplots.... 

library(ggpubr)
t_test_result <- t.test(leaves_se ~ SNP42, data = genotype_data_42)
print(t_test_result)
t_test_result <- t.test(leaves_se ~ SNP1104, data = genotype_data_1104)
print(t_test_result)


# If a p-value is less than 0.05, it is flagged with one star (*). If a p-value is less than 0.01, it is flagged with 2 stars (**). If a p-value is less than 0.001, it is flagged with three stars (***)

plot1 <-ggplot(genotype_data_42, aes(x = SNP42, y = leaves_se)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 1.1, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP42", x = "Genotype", y = "Standard error of Leaves") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )


plot2 <- ggplot(genotype_data_1104, aes(x = SNP1104, y = leaves_se)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 1.1, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP1104", x = "Genotype", y = "Standard error of Leaves") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )

(plot1|plot2)
t_test_result <- t.test(FT_mean ~ SNP1104, data = genotype_data_1104)
print(t_test_result)

#plot
svg("../../../My_Papers/photoperiod_shift/Figures/Boxplot_se_leaves.svg", width = 7.2, height = 6)
(plot1|plot2)

dev.off()


# as the best boxplit is looking at shappe want to plot the shape of the two groups of curves.......

#the trait is varmx_PC2_gc_LDLD_vel so load in the mean LDLD curves then group them by genotype 


#look at the split 
# Assuming the column with factor levels is called 'SNP1104'
grouped_rows <- split(genotype_data_42$Genotype, genotype_data_42$SNP42)

# Print row names for each group
grouped_rows


meancurves <- read.csv("genotype_means_LDLD.csv")
meancurves <- data.matrix(meancurves[,-1])
new_time <- seq(56,152,length.out=300)

#it was using the accelaeration curves do this....... 
BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_time),max(new_time)),norder = 4,nbasis = 300)

fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
smooth_both <- smooth.basis(argvals = new_time, y = meancurves, fdParobj = fdobj)
acc_obj <- deriv.fd(smooth_both$fd,2)


#evalulate 

eval_acc <- eval.fd(new_time,acc_obj,0)


##################
#group for ws 

WS_group <- eval_acc[,grouped_rows$W]
WS_group <- data.frame(new_time,WS_group)

TNZ_group <- eval_acc[,grouped_rows$X]
TNZ_group <- TNZ_group  # buffer
TNZ_group <- data.frame(new_time,TNZ_group)

# Calculate the mean and standard deviation for each group at each time point
WS_group_summary <- WS_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

TNZ_group_summary <- TNZ_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

# Combine summaries with a Source identifier
combined_summary <- bind_rows(
  WS_group_summary %>% mutate(Source = "Ws-2"),
  TNZ_group_summary %>% mutate(Source = "Tnz-1")
)


svg("../../../My_Papers/photoperiod_shift/Figures/shape_se_leaves.svg", width = 7.2, height = 6)

# Plot with ggplot
ggplot(combined_summary, aes(x = new_time, y = Mean, color = Source, fill = Source)) +
  geom_line(size = 1.2) +  # Plot mean lines
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2, color = NA) +  # Add shadow for noise
  geom_vline(xintercept = 104, linetype = "dashed", color = "black", size = 1) +  # Add dashed vertical line
  scale_color_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom line colors
  scale_fill_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom shadow colors
  labs(
    title = "",
    x = "Time (h)",
    y = "",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14,, color = "black"),  # Increase x-axis number size
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.text.y = element_blank(),  # Remove y-axis numbers
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Increase gap to plot edges
    legend.position = "right",  # Adjust legend position (optional)
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  ) + annotate("text", x = 104, y = 0.042, label = "Missed light cue", 
               hjust = 1.1, vjust = -1, size = 5, color = "black")

dev.off()


##############################################

#                          BIOMASS

###################################################################


########################################################################################



##       PLOTS OF QTLS


#read in saved output 
load("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/QTL_maping_lod_filtered_table.RData")


#   For chromosome 5
# want to have plots of meanBiomass and the shared traits so make outhk 0 and 1s to find the SNP regions 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")

looking <- out_hk

for (i in 3:ncol(looking)) {
  looking[,i] <- ifelse(looking[,i] > threshold_5_HK[(i-2)], 1, 0)
}

triat <- looking[,-c(1,2)]
triat<- triat[,colSums(triat, na.rm = TRUE) != 0]
#and remove rows that contain zero too
triat <- triat[rowSums(triat, na.rm = TRUE) != 0 ,]

#look for rows where mean leaves has a one 
traits_share_meanBiomass <- triat[triat[,277] != 0 ,]

#remove any rows or columns are zero
traits_share_meanBiomass<- traits_share_meanBiomass[,colSums(traits_share_meanBiomass, na.rm = TRUE) != 0]
#and remove rows that contain zero too
traits_share_meanBiomass <- traits_share_meanBiomass[rowSums(traits_share_meanBiomass, na.rm = TRUE) != 0 ,]

#remove any row that doesnt = 0 
#and remove rows that contain zero too wher trait is non-developmental
traits_share_meanBiomass <- traits_share_meanBiomass[rowSums(traits_share_meanBiomass[,1:26], na.rm = TRUE) != 0 ,]

# so there are 5 regions that qualifiy.... 201-218,295-299,1098,1100, 1104-1109.

# and 26 traits...
[1] "PC1_gc_SDLD_acc"                  "PC1_gc_SDSD_acc"                  "PC2_gc_SDSD_acc"                 
[4] "PC1_mean_dist_SDSD_ori"           "PC3_mean_dist_LDLD_ori"           "mean_total_euclid_SDSD_ori"      
[7] "varmx_PC2_gc_SDLD_acc"            "varmx_PC2_gc_SDSD_acc"            "varmx_PC3_gc_SDSD_acc"           
[10] "varmx_PC2_mean_dist_SDSD_ori"     "varmx_mean_total_euclid_SDSD_ori" "varmx_PC2_se_dist_SDSD_ori"      
[13] "varmx_SDLD_PC1_center_acc"        "PC3_gc_SDSD_post"                 "PC2_gc_vel_SDSD_post"            
[16] "PC1_gc_acc_SDSD_pre"              "PC2_gc_acc_SDSD_pre"              "PC2_gc_acc_SDSD_post"            
[19] "PC1_mean_dist_SDSD_ori_pre"       "PC3_mean_dist_LDLD_ori_pre"       "mean_total_euclid_LDLD_vel_pre"  
[22] "varmx_gc_warp_PC1_velSDLD"        "varmx_sd_PC1_LDLDpre"             "Amp_velLDLD"                     
[25] "gc_warp_PC1_velSDLD"              "varmx_gc_warp_PC1_velSDLD"        "FT_mean"                         
[28] "leaves_mean"                      "weights_mean"                        "weights_mean"                    
[31] "weights_se"



#   choose a trait that is on chr 2 and  chr 5 ideally one trait only on each so plot all 

chr5 <- traits_share_meanBiomass[20:24,]

#remove cols = 0 
#remove any rows or columns are zero

chr5<- chr5[,colSums(chr5, na.rm = TRUE) != 0]





varmx_gc_warp_PC1_velSDLD

#to see where the QTL starts not just the max . #then need to make just the columns I wat 
which(colnames(out_hk) == "weights_mean") 
which(colnames(out_hk) == "varmx_gc_warp_PC1_velSDLD") # this is shift 
which(colnames(out_hk) == "PC2_gc_acc_SDSD_pre") #this is shape 

#first do plot using all the information 
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")

ddsd <- out_hk[,colnames(traits_share_meanBiomass)]
# looking at QTL plots 
par(mar = c(5, 5, 4, 2))

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_chr5_meanBiomass.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "weights_mean") - 2,  col = "grey44",  ylab = "LOD score",  main = "Mean number of leaves", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "weights_mean") - 2], col = "grey44",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_biomass_phaseshift.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "varmx_gc_warp_PC1_velSDLD") - 2,  col = "red3",  ylab = "LOD score",  main = "Phase Shift", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "varmx_gc_warp_PC1_velSDLD") - 2], col = "red3",lwd=3)

dev.off()

svg("../../../My_Papers/photoperiod_shift/Figures/QTL_biomass_shape.svg", width = 7.2, height = 6)

plot(out_hk, ylim = c(0, 7.5), lodcolumn = which(colnames(out_hk) == "PC2_gc_acc_SDSD_pre") - 2,  col = "#6baed6",  ylab = "LOD score",  main = "Shape", cex.lab = 1.2,  cex.axis = 1.2,lwd=3)
abline(h = threshold_5_HK[which(colnames(out_hk) == "PC2_gc_acc_SDSD_pre") - 2], col = "#6baed6",lwd=3)

dev.off()


#
#You can do a dark grey for the physiological trait, light blue for spread of shift, pink for spread of shape and dark red for shift, or something like that
############### only chr5 and 1 
#b,l,t,r
varmx_gc_warp_PC1_velSDLD
PC2_gc_acc_SDSD_pre
#green #2171b5 darkred

par(mar=c(4.5,4.5,2,10))

plot(out_hk, ylim = c(0,9), lodcolumn = c(which(colnames(out_hk) == "weights_mean") - 2,which(colnames(out_hk) == "varmx_gc_warp_PC1_velSDLD") - 2,which(colnames(out_hk) == "PC2_gc_acc_SDSD_pre") - 2), col = c("grey44","red3","#6baed6"), chr = c(2,5), ylab = "LOD score", cex.lab = 1.2,  cex.axis = 1.2,lwd =2 )
abline(h=threshold_5_HK[which(colnames(out_hk) == "weights_mean") - 2], col="grey44",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "varmx_gc_warp_PC1_velSDLD") - 2], col="red3",lwd =2)
abline(h=threshold_5_HK[which(colnames(out_hk) == "PC2_gc_acc_SDSD_pre") - 2], col="#6baed6",lwd =2)
#add maf2 line 
abline(v=(66.16488966+25+129.7021294), col="black", lty="dashed", lwd =2)

#for FT 122.7670540134
#abline(v=(92.28854), col="black", lty=3, lwd =2)

#add ledgend 
legend(x=250, y = 4.5,
       inset = c(-0.4, 0),  # Move the legend outside the plot area
       legend = c("FT se", "Phase Shift", "Shape", "MAF2"), 
       col = c("grey44", "darkred", "#6baed6","black"), 
       lty = c(1, 1, 1,2), 
       lwd = c(3, 3, 3,3), 
       pt.cex = 2, 
       bty = "n", 
       xpd = TRUE,  # Allow drawing outside the plot region
       y.intersp = 1.5)

#saved at "QTL_leaves_grouped"

# try doing box plots! 

# read in the genotype data  

genotype_data <- read.csv("QTL_plots_genotype_info.csv")

#make sure columns are as factors 

genotype_data$SNP295 <- as.factor(genotype_data$SNP295)
genotype_data$SNP1098 <- as.factor(genotype_data$SNP1098)

#make plots for both the QTL regions 

#for FT_se there will be two as it has a QTL in both the regions 

genotype_data_295 <- genotype_data[genotype_data[[14]] %in% c("W", "X"), ]
genotype_data_1098 <- genotype_data[genotype_data[[15]] %in% c("W", "X"), ]
###
#to work out the significance levels of each!

#first do normaility test. If the p-value of the test is greater than α = .05, then the data is assumed to be normally distributed. Also check histagram and qqnorm plots

#weights_mean
shapiro.test(genotype_data_295$weights_mean)$p.value > 0.05
hist(genotype_data_295$weights_mean)
qqnorm(genotype_data_295$weights_mean)
qqline(genotype_data_295$weights_mean)

shapiro.test(genotype_data_1098$weights_mean)$p.value > 0.05
hist(genotype_data_1098$weights_mean)
qqnorm(genotype_data_1098$weights_mean)
qqline(genotype_data_1098$weights_mean)
#



shapiro.test(genotype_data_295$varmx_gc_warp_PC1_velSDLD)$p.value > 0.05
hist(genotype_data_295$varmx_gc_warp_PC1_velSDLD)
qqnorm(genotype_data_295$varmx_gc_warp_PC1_velSDLD)
qqline(genotype_data_295$varmx_gc_warp_PC1_velSDLD)
#
shapiro.test(genotype_data_1098$PC2_gc_acc_SDSD_pre)$p.value > 0.05
hist(genotype_data_1098$PC2_gc_acc_SDSD_pre)
qqnorm(genotype_data_1098$PC2_gc_acc_SDSD_pre)
qqline(genotype_data_1098$PC2_gc_acc_SDSD_pre)


#now do a t-test for each of the boxplots.... 

library(ggpubr)
t_test_result <- t.test(weights_mean ~ SNP295, data = genotype_data_295)
print(t_test_result)
t_test_result <- t.test(weights_mean ~ SNP1098, data = genotype_data_1098)
print(t_test_result)


# If a p-value is less than 0.05, it is flagged with one star (*). If a p-value is less than 0.01, it is flagged with 2 stars (**). If a p-value is less than 0.001, it is flagged with three stars (***)

plot1 <-ggplot(genotype_data_295, aes(x = SNP295, y = weights_mean)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 0.1, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP295", x = "Genotype", y = "Mean Biomass") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )


plot2 <- ggplot(genotype_data_1098, aes(x = SNP1098, y = weights_mean)) +
  geom_boxplot(width = 0.6, fill = "gray", color = "black", outlier.shape = NA) +  # Boxplot style
  #geom_jitter(width = 0.1, color = "black", size = 1.5) +  # Add individual data points
  stat_compare_means(method = "t.test", label.y = 0.1, label.x = 1) +
  stat_compare_means(method = "t.test", comparisons = list(c("W", "X")), label = "p.signif", vjust = 2, ref.group = 0.5) +  
  labs(title = "SNP1098", x = "Genotype", y = "Mean Biomass") +  # Labels
  scale_x_discrete(labels = c("Ws-2", "Tnz-1")) +  # Custom x-axis labels
  theme_classic() +  # Clean background
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(t = 15)),  # Increase space between x-axis title and plot
    axis.title.y = element_text(size = 12, face = "bold", colour = "black", 
                                margin = margin(r = 15)),  # Increase space between y-axis title and plot
    axis.text = element_text(size = 12, face = "bold.italic", colour = "black"),  # Axis tick font size
    axis.line = element_line(color = "black")  # Axis lines style
  )

(plot1|plot2)
t_test_result <- t.test(FT_mean ~ SNP1098, data = genotype_data_1098)
print(t_test_result)

#plot
svg("../../../My_Papers/photoperiod_shift/Figures/Boxplot_mean_biomass.svg", width = 7.2, height = 6)
(plot1|plot2)

dev.off()


# as the best boxplit is looking at shappe want to plot the shape of the two groups of curves.......

#the trait is PC2_gc_acc_SDSD_pre so load in the mean LDLD curves then group them by genotype 


#look at the split 
# Assuming the column with factor levels is called '1098'
grouped_rows <- split(genotype_data_1098$Genotype, genotype_data_1098$SNP1098)

# Print row names for each group
grouped_rows


meancurves <- read.csv("genotype_means_SDSD.csv")
meancurves <- data.matrix(meancurves[,-1])
new_time <- seq(56,152,length.out=300)

#it was using the accelaeration curves do this....... 
BASIS_bspline <- create.bspline.basis(rangeval = c(min(new_time),max(new_time)),norder = 4,nbasis = 300)

fdobj = fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
smooth_both <- smooth.basis(argvals = new_time, y = meancurves, fdParobj = fdobj)
acc_obj <- deriv.fd(smooth_both$fd,2)


#evalulate 

eval_acc <- eval.fd(new_time,acc_obj,0)


##################
#group for ws 

WS_group <- eval_acc[,grouped_rows$W]
WS_group <- data.frame(new_time,WS_group)

TNZ_group <- eval_acc[,grouped_rows$X]
TNZ_group <- TNZ_group  # buffer
TNZ_group <- data.frame(new_time,TNZ_group)

# Calculate the mean and standard deviation for each group at each time point
WS_group_summary <- WS_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

TNZ_group_summary <- TNZ_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

# Combine summaries with a Source identifier
combined_summary <- bind_rows(
  WS_group_summary %>% mutate(Source = "Ws-2"),
  TNZ_group_summary %>% mutate(Source = "Tnz-1")
)


svg("../../../My_Papers/photoperiod_shift/Figures/shape_mean_biomass.svg", width = 7.2, height = 6)

# Plot with ggplot
ggplot(combined_summary, aes(x = new_time, y = Mean, color = Source, fill = Source)) +
  geom_line(size = 1.2) +  # Plot mean lines
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2, color = NA) +  # Add shadow for noise
  geom_vline(xintercept = 104, linetype = "dashed", color = "black", size = 1) +  # Add dashed vertical line
  scale_color_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom line colors
  scale_fill_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom shadow colors
  labs(
    title = "",
    x = "Time (h)",
    y = "",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14,, color = "black"),  # Increase x-axis number size
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.text.y = element_blank(),  # Remove y-axis numbers
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Increase gap to plot edges
    legend.position = "right",  # Adjust legend position (optional)
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  ) + annotate("text", x = 104, y = 0.038, label = "Missed light cue", 
               hjust = 1.1, vjust = -1, size = 5, color = "black")

dev.off()




# for plot of SNP218 in terms of curves.... as is a shape parameter will look at shape of curves  but it is post shift so maybe is best to llok at the warping 


#########################################################################################################


#########################################################################################################



#look at the split 
# Assuming the column with factor levels is called 'SNP281'
grouped_rows <- split(genotype_data_281$Genotype, genotype_data_281$SNP281)

# Print row names for each group
grouped_rows

#remove wxT 31 as not in final map removeed... 

grouped_rows$W <- grouped_rows$W[-12] 

# load in the warping functions 

#origional curves
load(file = "../../../Warping/56_152/full_length_curves/origional/reg_fulllength_to_means_oriSDSD.RData")

reg_oriSDLD <- regSDSD

#56_152/full_length_curves/velocity/reg_fulllength_to_means_ve
#velocity curves
load(file = "../../../Warping/56_152/full_length_curves/velocity/reg_fulllength_to_means_velSDSD.RData")
reg_velSDLD <- regSDSD


#evalulate the curves so can plot them?

#first create mean for each 
means<- lapply(X = seq(1,64), FUN = function(x){mean.fd(reg_oriSDLD[[x]]$warpfd)})

#eval 

new_time <- seq(56,152, length.out = 300)

evlameans <- sapply(X = seq(1,64), FUN = function(x){eval.fd(evalarg = new_time,means[[x]],0)})

Genotype_names <- c("WxT_1.","WxT_2.","WxT_3.","WxT_4.","WxT_5.","WxT_6.","WxT_9.","WxT_10.","WxT_11.","WxT_12.","WxT_13.","WxT_14.","WxT_15.","WxT_16.","WxT_17.","WxT_18.","WxT_19.","WxT_21.","WxT_22.","WxT_23.","WxT_24.","WxT_25.","WxT_26.","WxT_27.","WxT_28.","WxT_29.","WxT_30.","WxT_31.","WxT_34.","WxT_35.","WxT_36.","WxT_38.","WxT_39.","WxT_40.","WxT_41.","WxT_42.","WxT_44.","WxT_45.","WxT_46.","WxT_47.","WxT_48.","WxT_49.","WxT_50.","WxT_51.","WxT_52.","WxT_55.","WxT_56.","WxT_57.","WxT_58.","WxT_59.","WxT_61.","WxT_62.","WxT_63.","WxT_64.","WxT_65.","WxT_66.","WxT_68.","WxT_69.","WxT_71.","WxT_73.","WxT_74.","WxT_76.","WxT_77.","WxT_78.")

Genotype_names <- sort(Genotype_names)
colnames(evlameans) <- Genotype_names

#group for ws 

WS_group <- evlameans[,grouped_rows$W]
WS_group <- WS_group - new_time
WS_group <- data.frame(new_time,WS_group)

TNZ_group <- evlameans[,grouped_rows$X]
TNZ_group <- TNZ_group - new_time
#TNZ_group <- TNZ_group +3 # buffer
TNZ_group <- data.frame(new_time,TNZ_group)

# Add an identifier column and combine the dataframes
df1_long <- WS_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  mutate(Source = "DF1")

df2_long <- TNZ_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  mutate(Source = "DF2")

combined_df <- bind_rows(df1_long, df2_long)

# Plot with ggplot
ggplot(combined_df, aes(x = new_time, y = Value, group = interaction(Variable, Source), color = Source)) +
  geom_line() +
  scale_color_manual(values = c("DF1" = "blue", "DF2" = "red")) +  # Custom colors
  labs(title = "Combined Dataframes Plot",
       x = "Time",
       y = "Value",
       color = "Source") +
  theme_minimal()
#


# Calculate the mean and standard deviation for each group at each time point
WS_group_summary <- WS_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

TNZ_group_summary <- TNZ_group %>%
  pivot_longer(-new_time, names_to = "Variable", values_to = "Value") %>%
  group_by(new_time) %>%
  summarise(Mean = mean(Value), SD = sd(Value))

# Combine summaries with a Source identifier
combined_summary <- bind_rows(
  WS_group_summary %>% mutate(Source = "Ws-2"),
  TNZ_group_summary %>% mutate(Source = "Tnz-1")
)

svg("../../../My_Papers/photoperiod_shift/Figures/SNP281_shape_mean_leaves.svg", width = 7.2, height = 6)

# Plot with ggplot
ggplot(combined_summary, aes(x = new_time, y = Mean, color = Source, fill = Source)) +
  geom_line(size = 1.2) +  # Plot mean lines
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2, color = NA) +  # Add shadow for noise
  geom_vline(xintercept = 104, linetype = "dashed", color = "black", size = 1) +  # Add dashed vertical line
  scale_color_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom line colors
  scale_fill_manual(values = c("Ws-2" = "blue", "Tnz-1" = "red")) +  # Custom shadow colors
  labs(
    title = "",
    x = "Time (h)",
    y = "",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14,, color = "black"),  # Increase x-axis number size
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.text.y = element_blank(),  # Remove y-axis numbers
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Increase gap to plot edges
    legend.position = "right",  # Adjust legend position (optional)
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  ) + annotate("text", x = 104, y = 1.25, label = "Missed light cue", 
               hjust = 1.1, vjust = -1, size = 5, color = "black")


dev.off()




