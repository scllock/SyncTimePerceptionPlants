# Input data ready for Elastic net =====


# TPM tables read in raw tables

Control_data <- fread("tnz_ws2_TPM_C.csv")
Treatment_data <- fread("tnz_ws2_TPM_T.csv")

#edit** remove genotypes that did not pass qtl analysis  
Control_data <- Control_data[,c("WT_3_C","WT_4_C","WT_6_C","WT_31_C","WT_54_C"):=NULL]
Treatment_data <- Treatment_data[,c("WT_3_T","WT_4_T","WT_6_T","WT_31_T", "WT_54_T"):=NULL]


#filter for genes where at least 50% of genotypes have a TPM value >10 in each condidtion separately


# Define the threshold and percentage criteria
TPMthreshold <- 10
percentage <- 0.5

# Calculate the number of samples
num_samples <- ncol(Control_data[,-1])

# Create a logical vector indicating whether each gene meets the criteria
genes_to_keep_C <- apply(Control_data[,-1], 1, function(row) {
  sum(row <= TPMthreshold) <= (num_samples * percentage)
})

genes_to_keep_T <- apply(Treatment_data[,-1], 1, function(row) {
  sum(row <= TPMthreshold) <= (num_samples * percentage)
})

# Filter each data frame to keep only the genes that meet the criteria
filtered_df_C <- Control_data[genes_to_keep_C, ]
filtered_df_T <- Treatment_data[genes_to_keep_T, ]

#find the intersection of these genes so always using the same list of genes in future analysis  

genes_to_keep <- union(filtered_df_C$gene_id,filtered_df_T$gene_id)

setkey(Control_data, gene_id)
setkey(Treatment_data, gene_id)

filtered_df_C <- Control_data[genes_to_keep]
filtered_df_T <- Treatment_data[genes_to_keep]

#save output

#write.csv(filtered_df_C, file = "TPM10_gene_expression_C.csv")

#write.csv(filtered_df_T, file = "TPM10_gene_expression_T.csv")



#get zcores of each TPM table by normalising each gene 

Control_TPM10 <- filtered_df_C
Treatment_TPM10 <- filtered_df_T


genenames <- Control_TPM10$gene_id
Zscore_Control_TPM10 <- get_z_score_matrix(Control_TPM10[,-1])
Zscore_Control_TPM10 <- data.frame(Zscore_Control_TPM10)
rownames(Zscore_Control_TPM10) <- genenames

genenames <- Treatment_TPM10$gene_id
Zscore_Treatment_TPM10 <- get_z_score_matrix(Treatment_TPM10[,-1])
Zscore_Treatment_TPM10 <- data.frame(Zscore_Treatment_TPM10)
rownames(Zscore_Treatment_TPM10) <- genenames


#save output 
#write.csv(Zscore_Control_TPM10, file = "Zscore_Control_TPM10.csv")
#write.csv(Zscore_Treatment_TPM10, file = "Zscore_Treatment_TPM10.csv")


## Look for differentially expressed genes in Ws-2 vs Tnz-1 across experiments.

# need to calculate the the logfold change between the control group and the treatment  for parents and use the common gene set

Control_TPM10_WS <-  Control_TPM10$WS_2_C
Control_TPM10_TNZ <- Control_TPM10$TNZ_1_C

Treatment_TPM10_WS <- Treatment_TPM10$WS_2_T
Treatment_TPM10_TNZ <- Treatment_TPM10$TNZ_1_T

# Add a small constant to avoid log(0) typically RNA use 1
pseudocount <- 1

# Calculate log fold change

log10_fc_ws_ws_CT <- log10(Control_TPM10_WS + pseudocount) - log10(Treatment_TPM10_WS + pseudocount)
log10_fc_tnz_tnz_CT <- log10(Control_TPM10_TNZ + pseudocount) - log10(Treatment_TPM10_TNZ + pseudocount)

#combineinto matrix and convert to dataframe

logFC_Ws_Tnz <- cbind(Control_TPM10[,1],log10_fc_ws_ws_CT,log10_fc_tnz_tnz_CT )
logFC_Ws_Tnz <- as.data.frame(logFC_Ws_Tnz)

# Ensure the X and Y columns are numeric
logFC_Ws_Tnz$log10_fc_ws_ws_CT <- as.numeric(logFC_Ws_Tnz$log10_fc_ws_ws_CT)
logFC_Ws_Tnz$log10_fc_tnz_tnz_CT <- as.numeric(logFC_Ws_Tnz$log10_fc_tnz_tnz_CT)

#calulate the eluclidean distance for distane (sqrt(x^2+y^2)>1)
euclid_dist_Ws_Tnz <- sqrt((log10_fc_ws_ws_CT)^2 + (log10_fc_tnz_tnz_CT)^2)

# Find the threshold for the top 20% distances
threshold_Ws_Tnz <- quantile(euclid_dist_Ws_Tnz, 0.8)

# Keep only the points with the top 20% distances
Ws_Tnz_top_20 <- logFC_Ws_Tnz[euclid_dist_Ws_Tnz > threshold_Ws_Tnz, ]

Ws_Tnz_top_20 <- as.data.frame(Ws_Tnz_top_20)
Ws_Tnz_top_20$log10_fc_ws_ws_CT <- as.numeric(Ws_Tnz_top_20$log10_fc_ws_ws_CT)
Ws_Tnz_top_20$log10_fc_tnz_tnz_CT <- as.numeric(Ws_Tnz_top_20$log10_fc_tnz_tnz_CT)

#plot of all points and the ones to include in alanysis in blue  
ggplot(logFC_Ws_Tnz, aes(x = log10_fc_ws_ws_CT, y = log10_fc_tnz_tnz_CT)) +
  geom_point(size = 0.8) +
  theme_classic() +
  labs(
    x = expression(Log[10] ~ "Fold Change WS Control and Treatment"),
    y = expression(Log[10] ~ "Fold Change TNZ Control and Treatment"),
    title = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 14),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
    axis.text.x = element_text(hjust = 1, size = 12, colour = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  ) +
  geom_point(data = Ws_Tnz_top_20, aes(x = log10_fc_ws_ws_CT, y = log10_fc_tnz_tnz_CT), 
             color = "blue", size = 0.8)

# filter each origional dataset for the genes that pass the threshold. 

Zscore_DEG_Ws_Tnz_control <- Zscore_Control_TPM10[euclid_dist_Ws_Tnz > threshold_Ws_Tnz, ]
Zscore_DEG_Ws_Tnz_treatment <- Zscore_Treatment_TPM10[euclid_dist_Ws_Tnz > threshold_Ws_Tnz, ]

#save files 
write.csv(Zscore_DEG_Ws_Tnz_control, "Zscore_DEG_Ws_Tnz_control.csv")
write.csv(Zscore_DEG_Ws_Tnz_treatment, "Zscore_DEG_Ws_Tnz_treatment.csv")

#These files are used for the elastic net input. 