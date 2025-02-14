# load elastic net results ==================

l1_ratio_meanft_C <- read_excel("Results/DEG_Tnz_Ws_new//Control/analysis/hyperparam_lists/meanft_hyperparams.xlsx")

l1_ratio_meanleaves_C <- read_excel("Results/DEG_Tnz_Ws_new//Control/analysis/hyperparam_lists/leavesmean_hyperparams.xlsx")

l1_ratio_meanbiomass_C <- read_excel("Results/DEG_Tnz_Ws_new//Control/analysis/hyperparam_lists/biomassmean_hyperparams.xlsx")


l1_ratio_meanft_T <- read_excel("Results/DEG_Tnz_Ws_new//Treatment/analysis/hyperparam_lists/meanft_hyperparams.xlsx")

l1_ratio_meanleaves_T <- read_excel("Results/DEG_Tnz_Ws_new//Treatment/analysis/hyperparam_lists/leavesmean_hyperparams.xlsx")

l1_ratio_meanbiomass_T <- read_excel("Results/DEG_Tnz_Ws_new//Treatment/analysis/hyperparam_lists/biomassmean_hyperparams.xlsx")

#control
median(l1_ratio_meanft_C$l1_ratio)
median(l1_ratio_meanleaves_C$l1_ratio)
median(l1_ratio_meanbiomass_C$l1_ratio)

#treatment
median(l1_ratio_meanft_T$l1_ratio)
median(l1_ratio_meanleaves_T$l1_ratio)
median(l1_ratio_meanbiomass_T$l1_ratio)


#Violin plots ========== 

#load gene list from ensembl 

listEnsemblGenomes()
ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")
searchDatasets(ensembl_plants, pattern = "Arabidopsis")
ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart",
                                         dataset = "athaliana_eg_gene")


#TNZ vs WS 
#Here change the file read in for different traits (meanft, meanleaves ect) and diffenet conditions. 

#Control
EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Control/analysis/coefficients_lists/meanft_coefs.xlsx")
#there are 63 genotypes so want every column where the number of non zero is > 31
#Identify columns where the first row has values > 30
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)

#Subset the data frame to keep only those columns
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]

#remove the unwanted rows  
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz_filtered[-c(1,2,3,4,5),-1]

#mean of each column 
mean_coeff <- apply(EN_coef_ws_tnz_filtered, 2,mean)
mean_coeff <-data.frame(mean_coeff)

# Calculate the magnitude of change (e.g., max change in each row)
magchange <- apply(EN_coef_ws_tnz_filtered, 2, function(x) max(abs(x)))

# Find the top 10 highest values and their names
top_10 <- names(sort(magchange, decreasing = TRUE)[1:10])

#filter just the top 10 genes 
EN_coef_ws_tnz_toplot<- as.data.frame(t(EN_coef_ws_tnz_filtered[,top_10]))

colnames(EN_coef_ws_tnz_toplot) <- colnames(Control_TPM10[,-1])

#here use the sort hand gene names.......... 
#use row names from the top 10 genes 

gene_info <- getBM(attributes = c("ensembl_gene_id", "description", "external_gene_name", "chromosome_name"),
                   filters = "ensembl_gene_id",   # Search by gene names
                   values = rownames(EN_coef_ws_tnz_toplot),            # Your gene list
                   mart = ensembl_arabidopsis)

# Add rownames as a new column to keep gene IDs
EN_coef_ws_tnz_toplot$gene_id <- gene_info$external_gene_name
#rename using the short names.....
rownames(EN_coef_ws_tnz_toplot) <- gene_info$external_gene_name


# Reshape the data from wide to long format
df_long <- gather(EN_coef_ws_tnz_toplot, key = "condition", value = "expression", -gene_id)


#svg("../../My_Papers/photoperiod_shift/Figures/violin_meanft_C.svg", width = 7.2, height = 6)

# Create the violin plot ordring based on size 
ggplot(df_long, aes(x = reorder(gene_id, expression, FUN = median), y = expression)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") + 
  theme_classic() + 
  labs(x = "Gene ID", y = "Coefficients of gene expression") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.y = element_text(margin = margin(r = 20), size = 16),
        axis.title.x = element_text(margin = margin(t = 20), size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, color = "black"))

#save as high aquality.... 
#dev.off()

#Treatment 
EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Treatment/analysis/coefficients_lists/meanft_coefs.xlsx")
#there are 63 genotypes so want every column where the number of noon zero is > 31
# Step 1: Identify columns where the first row has values > 30
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)

# Step 2: Subset the data frame to keep only those columns
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]

#remove the bits at the top 

EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz_filtered[-c(1,2,3,4,5),-1]


# Calculate the magnitude of change (e.g., max change in each row)
magchange <- apply(EN_coef_ws_tnz_filtered, 2, function(x) max(abs(x)))

# Find the top 10 highest values and their names
top_10 <- names(sort(magchange, decreasing = TRUE)[1:10])


#filter just the top 10 genes 
EN_coef_ws_tnz_toplot<- as.data.frame(t(EN_coef_ws_tnz_filtered[,top_10]))

colnames(EN_coef_ws_tnz_toplot) <- colnames(Control_TPM10[,-1])

#use row names from the top 10 genes 

gene_info <- getBM(attributes = c("ensembl_gene_id", "description", "external_gene_name", "chromosome_name"),
                   filters = "ensembl_gene_id",   # Search by gene names
                   values = rownames(EN_coef_ws_tnz_toplot),            # Your gene list
                   mart = ensembl_arabidopsis)

# Add rownames as a new column to keep gene IDs
EN_coef_ws_tnz_toplot$gene_id <- gene_info$external_gene_name
#rename using the short names.....
rownames(EN_coef_ws_tnz_toplot) <- gene_info$external_gene_name


# Reshape the data from wide to long format
df_long <- gather(EN_coef_ws_tnz_toplot, key = "condition", value = "expression", -gene_id)


#svg("../../My_Papers/photoperiod_shift/Figures/violin_meanft_T.svg", width = 7.2, height = 6)

# Create the violin plot ordring based on size 
ggplot(df_long, aes(x = reorder(gene_id, expression, FUN = median), y = expression)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") + 
  theme_classic() + 
  labs(x = "Gene ID", y = "Coefficients of gene expression") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.y = element_text(margin = margin(r = 20), size = 16),
        axis.title.x = element_text(margin = margin(t = 20), size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, color = "black"))

#save as high aquality.... 
#dev.off()
















################################

############################################


# again but dont filter for top 10 just do all ....
setwd("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/Elastic_Net/Transcriptome_elasticnet")

#DEG WS TNZ
EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Control/analysis/coefficients_lists/meanft_coefs.xlsx")
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]
genes9 <- colnames(EN_coef_ws_tnz_filtered[,-1])

EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Treatment/analysis/coefficients_lists/meanft_coefs.xlsx")
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]
genes10 <- colnames(EN_coef_ws_tnz_filtered[,-1])

EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Control/analysis/coefficients_lists/leavesmean_coefs.xlsx")
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]
genes11 <- colnames(EN_coef_ws_tnz_filtered[,-1])

EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Treatment/analysis/coefficients_lists/leavesmean_coefs.xlsx")
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]
genes12 <- colnames(EN_coef_ws_tnz_filtered[,-1])

############################################
EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Control/analysis/coefficients_lists/biomassmean_coefs.xlsx")
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]
genes13 <- colnames(EN_coef_ws_tnz_filtered[,-1])

EN_coef_ws_tnz <- read_excel("Results/DEG_Tnz_Ws_fixed//Treatment/analysis/coefficients_lists/biomassmean_coefs.xlsx")
cols_to_keep <- which(EN_coef_ws_tnz[1, ] > 30)
EN_coef_ws_tnz_filtered <- EN_coef_ws_tnz[, cols_to_keep]
genes14 <- colnames(EN_coef_ws_tnz_filtered[,-1])

length(genes9) <- 1200 
length(genes10) <- 1200   
length(genes11) <- 1200   
length(genes12) <- 1200   
length(genes13) <- 1200   
length(genes14) <- 1200   

# Combine vectors into a data frame

df <- data.frame(genes9,genes10,genes11,genes12,genes13,genes14)

colnames(df) <- c("DEG_ws_tnz_Control_meanft","DEG_ws_tnz_Treatment_meanft","DEG_ws_tnz_Control_meanleaves","DEG_ws_tnz_Treatment_meanleaves","DEG_ws_tnz_Control_meanbiomass","DEG_ws_tnz_Treatment_meanbiomass")
#save gene of intrest 

#write.csv(df, "Results/genes_of_intrest_post_EN_DEG_ws_tnz.csv")


# #DEG ws-2 vs tnz

wstnz_biomass_Control <- read_excel("Results/DEG_Tnz_Ws_fixed///Control/analysis/predictions/biomassmean_predictions.xlsx")
wstnz_biomass_Treatment <- read_excel("Results/DEG_Tnz_Ws_fixed//Treatment/analysis/predictions/biomassmean_predictions.xlsx")


# Calculate Spearman's correlation for biomass as has strong outlier 
spearman_correlation_wstnz_control <- cor(wstnz_biomass_Control$true, wstnz_biomass_Control$predicted, method = "spearman")

spearman_correlation_wstnz_treatment <- cor(wstnz_biomass_Treatment$true, wstnz_biomass_Treatment$predicted, method = "spearman")


############

# Create the data frame
r_squared_data <- data.frame(
  Variable = rep(c("FT", "Leaves", "Biomass"), each = 2),
  Treatment = rep(c("Light", "Dark"), 3),
  R_squared = c(0.76, 0.83, 0.61, 0.72, 0.68, 0.70)
)

# Set the order for the 'Variable' factor
r_squared_data$Variable <- factor(r_squared_data$Variable, levels = c("FT", "Leaves", "Biomass"))

# Create the bar plot
ggplot(r_squared_data, aes(x = Variable, y = R_squared, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Light" = "skyblue", "Dark" = "steelblue")) +
  theme_minimal() +
  labs(
    title = "", 
    x = "Trait", 
    y = "R-squared", 
    fill = "Treatment"
  )





############################


#                       SCATTER  PLOT



# sctter plot of coefficients 

#Here change the file read in for different traits (meanft, meanleaves ect) and different conditions 
#Control
setwd("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/Elastic_Net/Transcriptome_elasticnet/")

########### try without the filtering 
EN_coef_ws_tnz_control <- read_excel("Results/DEG_Tnz_Ws_fixed/Control/analysis/coefficients_lists/meanft_coefs.xlsx")

EN_coef_ws_tnz_treatment <- read_excel("Results/DEG_Tnz_Ws_fixed/Treatment/analysis/coefficients_lists/meanft_coefs.xlsx")

#remove the bits at the top 
cols_to_keep_control <- which(EN_coef_ws_tnz_control[1, ] > 30)
cols_to_keep_treatment <- which(EN_coef_ws_tnz_treatment[1, ] > 30)


EN_coef_ws_tnz_filtered_control <- EN_coef_ws_tnz_control[, cols_to_keep_control]
EN_coef_ws_tnz_filtered_treatment <- EN_coef_ws_tnz_treatment[, cols_to_keep_treatment]
#remove the bits at the top 

EN_coef_ws_tnz_filtered_control <- EN_coef_ws_tnz_filtered_control[-c(1,2,3,4,5),]
EN_coef_ws_tnz_filtered_treatment <- EN_coef_ws_tnz_filtered_treatment[-c(1,2,3,4,5),]

names(EN_coef_ws_tnz_filtered_control)[names(EN_coef_ws_tnz_filtered_control) == '...1'] <- 'Genotype'
names(EN_coef_ws_tnz_filtered_treatment)[names(EN_coef_ws_tnz_filtered_treatment) == '...1'] <- 'Genotype'

#want  to take the average of the each column 
EN_coef_ws_tnz_filtered_control_mean <- apply(EN_coef_ws_tnz_filtered_control[,-1], 2, function(x) mean(x, na.rm = TRUE))
EN_coef_ws_tnz_filtered_treatment_mean <- apply(EN_coef_ws_tnz_filtered_treatment[,-1], 2, function(x) mean(x, na.rm = TRUE))

# with gene names as column names and one row of values.

# 1. Get the union of the gene names (column names)
all_genes <- union(names(EN_coef_ws_tnz_filtered_control_mean), names(EN_coef_ws_tnz_filtered_treatment_mean))

# 2. Add missing genes to control_df  and treatment df and fill with 0s
EN_coef_ws_tnz_filtered_control_mean[setdiff(all_genes, names(EN_coef_ws_tnz_filtered_control_mean))] <- 0

EN_coef_ws_tnz_filtered_treatment_mean[setdiff(all_genes, names(EN_coef_ws_tnz_filtered_treatment_mean))] <- 0

# 4. Reorder both dataframes to ensure columns match
control_df <- EN_coef_ws_tnz_filtered_control_mean[all_genes]
treatment_df <- EN_coef_ws_tnz_filtered_treatment_mean[all_genes]

# 5. Create a scatter plot: Control on x-axis, Treatment on y-axis

# Convert the dataframes to long format for ggplot
# Assuming control_df and treatment_df are named vectors

df_long <- data.frame(
  Gene = names(control_df),         # Gene names
  Control = as.numeric(control_df),  # Control values
  Treatment = as.numeric(treatment_df)  # Treatment values
)

# Create the scatter plot
ggplot(df_long, aes(x = Control, y = Treatment)) +
  geom_point(color = "black", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +  # Add y = x line
  labs(title = "",
       x = "Coefficients in Control", 
       y = "Coefficients in Treatment") +
  theme_classic() + 
  theme(axis.title.y = element_text(margin = margin(r = 10), size = 14),
        axis.title.x = element_text(margin = margin(t = 10), size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12, colour = "black"),
        axis.text.y = element_text(size = 14, color = "black"))
#









