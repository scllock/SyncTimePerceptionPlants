#### transcriptome associated plots ==========

# for the expression data plots for paper



# looking at only the WS and TNZ from each sample.......


# file prep for elastic nets .....
# for the expression data 
library(dplyr)
library(data.table)
library(matrixStats)
library(ggplot2)
library(readxl)
library(tidyr)
library(pheatmap)
library(VennDiagram)


#read in the Zscore files for min TPM

setwd("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/Elastic_Net/Transcriptome_elasticnet/Results/")


#save 
Control_TPM10 <- read.csv(file = "minTPM/Treatment/Treatment/data/Zscore_Control_TPM10.csv")

Treatment_TPM10 <- read.csv( file = "minTPM/Treatment/Treatment/data/Zscore_Treatment_TPM10.csv")

Control_TPM10 <- read.csv("../../../Transcriptome_expression_data/Master_files/TPM10_gene_expression_C.csv")
#Control_TPM10 <- Control_TPM10[,-c(1,4,6)]
#rownames(Control_TPM10) <- Control_TPM10[,1]
#Control_TPM10 <- Control_TPM10[,-1]

#remove the parents from treatment 



Treatment_TPM10 <- read.csv("../../../Transcriptome_expression_data/Master_files/TPM10_gene_expression_T.csv")
#Treatment_TPM10 <- Treatment_TPM10[,-1]


##         Differentially expressed genes in WS vs TNz  

# wnt to look at the scatterplot and then take the euclidiean distance and take top 10% of genes and filter each of the groups for this gene list 

# need to calculate the the fold change between the control group and the treatment use the common gene set

Control_TPM10_WS <-  Control_TPM10$WS_2_C
Control_TPM10_TNZ <- Control_TPM10$TNZ_1_C

Treatment_TPM10_WS <- Treatment_TPM10$WS_2_T
Treatment_TPM10_TNZ <- Treatment_TPM10$TNZ_1_T

# Add a small constant to avoid log(0) typically RNA use 1
pseudocount <- 1

# Calculate log fold change
log10_fc_ws_tnz_C <- log10(Control_TPM10_WS + pseudocount) - log10(Control_TPM10_TNZ + pseudocount)

log10_fc_ws_tnz_T <- log10(Treatment_TPM10_WS + pseudocount) - log10(Treatment_TPM10_TNZ + pseudocount)

#then

log10_fc_ws_ws_CT <- log10(Control_TPM10_WS + pseudocount) - log10(Treatment_TPM10_WS + pseudocount)


log10_fc_tnz_tnz_CT <- log10(Control_TPM10_TNZ + pseudocount) - log10(Treatment_TPM10_TNZ + pseudocount)

#plot 

#A use control filtered for gene names 
logFC_C_T <- cbind(Control_TPM10[,1],log10_fc_ws_tnz_C,log10_fc_ws_tnz_T)

# Convert the matrix to a data frame

logFC_C_T <- as.data.frame(logFC_C_T)

# Ensure the X and Y columns are numeric
logFC_C_T$log10_fc_ws_tnz_C <- as.numeric(logFC_C_T$log10_fc_ws_tnz_C)
logFC_C_T$log10_fc_ws_tnz_T <- as.numeric(logFC_C_T$log10_fc_ws_tnz_T)

#calculate the eluclidean distance  for control vs Treatment

# (sqrt(x^2+y^2)>1)
euclid_dist_C_T <- sqrt((log10_fc_ws_tnz_C)^2 + (log10_fc_ws_tnz_T)^2)

# Find the threshold for the top 20% distances
threshold <- quantile(euclid_dist_C_T, 0.8)

# Keep only the points with the top 10% distances
C_T_top_20 <- logFC_C_T[euclid_dist_C_T > threshold, ]

C_T_top_20 <- as.data.frame(C_T_top_20)
C_T_top_20$log10_fc_ws_tnz_C <- as.numeric(C_T_top_20$log10_fc_ws_tnz_C)
C_T_top_20$log10_fc_ws_tnz_T <- as.numeric(C_T_top_20$log10_fc_ws_tnz_T)

#plot of all points and the ones to include in alanysis in blue  

ggplot(logFC_C_T, aes(x = log10_fc_ws_tnz_C, y = log10_fc_ws_tnz_T)) +
  geom_point(size = 0.8) + theme_classic() +
  labs(x = "Log10 Fold Change Ws-2 and Tnz Dark", y = "Log10 Fold Change Ws-2 and Tnz Light", title = "") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+geom_vline(xintercept=0, linetype="dashed", color = "red") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_point(data = C_T_top_20, aes(x = log10_fc_ws_tnz_C, y = log10_fc_ws_tnz_T), color = "blue", size = 0.8)


#     For WS vs TNZ ======================
logFC_Ws_Tnz <- cbind(Control_TPM10[,1],log10_fc_ws_ws_CT,log10_fc_tnz_tnz_CT )
# Convert the matrix to a data frame

logFC_Ws_Tnz <- as.data.frame(logFC_Ws_Tnz)

# Ensure the X and Y columns are numeric
logFC_Ws_Tnz$log10_fc_ws_ws_CT <- as.numeric(logFC_Ws_Tnz$log10_fc_ws_ws_CT)
logFC_Ws_Tnz$log10_fc_tnz_tnz_CT <- as.numeric(logFC_Ws_Tnz$log10_fc_tnz_tnz_CT)


euclid_dist_Ws_Tnz <- sqrt((log10_fc_ws_ws_CT)^2 + (log10_fc_tnz_tnz_CT)^2)

# Find the threshold for the top 10% distances
threshold_Ws_Tnz <- quantile(euclid_dist_Ws_Tnz, 0.9)

# Keep only the points with the top 10% distances
logFC_Ws_Tnz[,1] <- Control_TPM10$V1

Ws_Tnz_top_20 <- logFC_Ws_Tnz[euclid_dist_Ws_Tnz > threshold_Ws_Tnz, ]

Ws_Tnz_top_20 <- as.data.frame(Ws_Tnz_top_20)
Ws_Tnz_top_20$log10_fc_ws_ws_CT <- as.numeric(Ws_Tnz_top_20$log10_fc_ws_ws_CT)
Ws_Tnz_top_20$log10_fc_tnz_tnz_CT <- as.numeric(Ws_Tnz_top_20$log10_fc_tnz_tnz_CT)

#plot of all points and the ones to include in alanysis in blue  

ggplot(logFC_Ws_Tnz, aes(x = log10_fc_ws_ws_CT, y = log10_fc_tnz_tnz_CT)) +
  geom_point(size = 0.8) + theme_classic() +
  labs(x = "Log10 Fold Change Ws-2 Dark and Light", y = "Log10 Fold Change Tnz Dark and Light", title = "") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+geom_vline(xintercept=0, linetype="dashed", color = "red") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_point(data = Ws_Tnz_top_20, aes(x = log10_fc_ws_ws_CT, y = log10_fc_tnz_tnz_CT), color = "blue", size = 0.8)


# plots of TPM values ===============

# Tnz c vs tnz t 
TNZ_CT <- cbind(Control_TPM10$TNZ_1_C,Treatment_TPM10$TNZ_1_T )
TNZ_CT <- as.data.frame(TNZ_CT)
colnames(TNZ_CT) <- c("TNZ_1_C","TNZ_1_T")
rownames(TNZ_CT) <- Control_TPM10[,1]
TNZ_CT$TNZ_1_C <- log10(as.numeric(TNZ_CT$TNZ_1_C))
TNZ_CT$TNZ_1_T <- log10(as.numeric(TNZ_CT$TNZ_1_T))

ggplot(TNZ_CT, aes(x = TNZ_1_C, y = TNZ_1_T)) +
  geom_point(size = 0.8) + theme_classic() + xlim(0, 5) + ylim(0,5) +
  labs(x = "TPM Tnz Dark", y = "TPM Tnz Light", title = "") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+geom_vline(xintercept=0, linetype="dashed", color = "red") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_abline(slope=1, intercept = 0, col = "red", linewidth = 0.8)

# ws c vs ws t 

WS_CT <- cbind(Control_TPM10$WS_2_C,Treatment_TPM10$WS_2_T )
WS_CT <- as.data.frame(WS_CT)
colnames(WS_CT)<- c("WS_2_C","WS_2_T")
WS_CT$WS_2_C <- log10(as.numeric(WS_CT$WS_2_C))
WS_CT$WS_2_T <- log10(as.numeric(WS_CT$WS_2_T))

ggplot(WS_CT, aes(x = WS_2_C, y = WS_2_T)) +
  geom_point(size = 0.8) + theme_classic() + xlim(0, 5) + ylim(0,5) +
  labs(x = "TPM Ws-2 Dark", y = "TPM Ws-2 Light", title = "") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+geom_vline(xintercept=0, linetype="dashed", color = "red") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_abline(slope=1, intercept = 0, col = "red", linewidth = 0.8)

# ws c vs tnz c 

WS_TNZ_C <- cbind(Control_TPM10$WS_2_C,Control_TPM10$TNZ_1_C)
WS_TNZ_C <- as.data.frame(WS_TNZ_C)
colnames(WS_TNZ_C)<- c("WS_2_C","TNZ_1_C")
WS_TNZ_C$WS_2_C <- log10(as.numeric(WS_TNZ_C$WS_2_C))
WS_TNZ_C$TNZ_1_C <- log10(as.numeric(WS_TNZ_C$TNZ_1_C))


ggplot(WS_TNZ_C, aes(x = WS_2_C, y = TNZ_1_C)) +
  geom_point(size = 0.8) + theme_classic() +  xlim(0, 5) + ylim(0,5) +
  labs(x = "TPM Ws-2 Dark", y = "TPM Tnz Dark", title = "") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+geom_vline(xintercept=0, linetype="dashed", color = "red") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_abline(slope=1, intercept = 0, col = "red", linewidth = 0.8)


# ws T vs tnz T 

WS_TNZ_T <- cbind(Treatment_TPM10$WS_2_T,Treatment_TPM10$TNZ_1_T)
WS_TNZ_T <- as.data.frame(WS_TNZ_T)
colnames(WS_TNZ_T) <- c("WS_2_T","TNZ_1_T")
WS_TNZ_T$WS_2_T <- log10(as.numeric(WS_TNZ_T$WS_2_T))
WS_TNZ_T$TNZ_1_T <- log10(as.numeric(WS_TNZ_T$TNZ_1_T))



ggplot(WS_TNZ_T, aes(x = WS_2_T, y = TNZ_1_T)) +
  geom_point(size = 0.8) + theme_classic() + xlim(0, 5) + ylim(0,5) +
  labs(x = "TPM Ws-2 Light", y = "TPM Tnz Light", title = "") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+geom_vline(xintercept=0, linetype="dashed", color = "red") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + geom_abline(slope=1, intercept = 0, col = "red", linewidth = 0.8)


#   HEATMAPS =============

# read in dataframe
control_filtered <- Control_TPM10[,-1]
treatment_filtered <- Treatment_TPM10[,-1]

# Remove gene column
data_control <- control_filtered[, -1]
data_treatment <- treatment_filtered[, -1]

# Combine data
combined_data <- cbind(data_control, data_treatment)

# Add rownames
rownames(combined_data) <- control_filtered$gene_id

# Create annotation for the columns
annotation_col <- data.frame(
  Condition = rep(c("Dark", "Light"), each=ncol(data_control)),
  row.names = colnames(combined_data)
)

# Generate the heatmap
pheatmap(combined_data, annotation_col = annotation_col, 
         scale = "row", # scale genes
         legend_labels = c("-10", "-5", "0", "5","10", "title\n"),
         clustering_method = "complete", # clustering method
         color = colorRampPalette(c("blue", "white", "red"))(50), show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 30, cluster_rows = TRUE, cluster_cols = TRUE) # color scheme

#next want to produce heat maps of the genes that are differentially expressed between the control and treatment group and between WS and TNZ. 


#use gene list from the top 20% of eulcidiean distance filter combined data ased on the gene names that are in C_T_top_20 

combined_data_CT_top20 <- combined_data[euclid_dist_C_T > threshold, ]

# Create annotation for the columns
annotation_col <- data.frame(
  Condition = rep(c("Dark", "Light"), each=ncol(data_control)),
  row.names = colnames(combined_data_CT_top20)
)

# Generate the heatmap
pheatmap(combined_data_CT_top20, annotation_col = annotation_col, 
         scale = "row", # scale genes
         legend_labels = c("-10", "-5", "0", "5","10", "title\n"),
         color = colorRampPalette(c("blue", "white", "red"))(50), show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 30, cluster_rows = TRUE, cluster_cols = TRUE) # color scheme

#look at genes different between WS and TNZ 

#filter combined data ased on the gene names that are in C_T_top_20 

combined_data_WsTnz_top20 <- combined_data[euclid_dist_Ws_Tnz > threshold_Ws_Tnz, ]

# Create annotation for the columns
annotation_col <- data.frame(
  Condition = rep(c("Dark", "Light"), each=ncol(data_control)),
  row.names = colnames(combined_data_WsTnz_top20)
)

# Generate the heatmap
pheatmap(combined_data_WsTnz_top20, annotation_col = annotation_col, 
         scale = "row", # scale genes
         color = colorRampPalette(c("blue", "white", "red"))(50), show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 30, cluster_rows = TRUE, cluster_cols = TRUE) # color scheme









