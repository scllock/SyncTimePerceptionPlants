# Load physiological data =================

print('Loading phsyiological data')

#alldata_ADJUSTED1.csv is exp 1 done in July and alldata_ADJUSTED2.csv is exp 2 done in November
data_july<- read.csv("alldata_ADJUSTED1.csv")
data_nov<- read.csv("alldata_ADJUSTED2.csv")


#check str of data
str(data_nov)
str(data_july)


#pivot tables to be one column per genotype so can compare between experiments

#looking at Flowering time first.
FTonly1 <- data_july[,1:2]
FTonly2 <- data_nov[,1:2]

widerdata1<- FTonly1 %>%
  group_by(Genotype) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Days.to.flower) %>%
  dplyr::select(-row)

widerdata2<- FTonly2 %>%
  group_by(Genotype) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Days.to.flower) %>%
  dplyr::select(-row)

#then order new datafrmaes by column then can compare between 2 experiments easier
new_order1 <- sort(colnames(widerdata1))
orderexp1 <- widerdata1[, new_order1]
orderexp1_mean <- apply(orderexp1,2,mean, na.rm = TRUE)
orderexp1_sd <- apply(orderexp1,2,sd, na.rm = TRUE)
orderexp1_se <- apply(orderexp1,2,standard_error)


new_order2 <- sort(colnames(widerdata2))
orderexp2 <- widerdata2[, new_order2]
orderexp2_mean <- apply(orderexp2,2,mean, na.rm = TRUE)
orderexp2_sd <- apply(orderexp2,2,sd, na.rm = TRUE)
orderexp2_se <- apply(orderexp2,2,standard_error)


combine <- rbind(orderexp1,orderexp2)
all_mean <- apply(combine,2,mean, na.rm = TRUE)
all_mean <- as.data.frame(all_mean)
all_sd <- apply(combine,2,sd, na.rm = TRUE)
all_se <- apply(combine, 2, standard_error)
all_se <- as.data.frame(all_se)


all_mean_FT <- all_mean
all_se_FT <- all_se
# Create a data frame from your data
data <- data.frame(orderexp1_mean, orderexp2_mean)

# Fit the linear model
best <- lm(orderexp2_mean ~ orderexp1_mean, data = data)

# Plot using ggplot2

# Calculate the R^2 value
r_squared <- summary(best)$r.squared

# Format R^2 text with superscript
r_squared_text <- bquote(R^2 == .(round(r_squared, 2)))

# Plot using ggplot2
ggplot(data, aes(x = orderexp1_mean, y = orderexp2_mean)) +
  geom_point(size = 2) +  # Plot the points
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "black", linewidth = 1.2) +  # Add the linear model line
  labs(x = "Days to Flower Exp 1", y = "Days to Flower Exp 2") +
  theme_classic() +
  theme(
    # Adjust axis title margin
    axis.title.x = element_text(margin = margin(t = 10), size = 14),
    axis.title.y = element_text(margin = margin(r = 10), size = 14),
    # Increase axis text size
    axis.text = element_text(size = 14, color = "black")
  ) +
  # Add R^2 annotation with superscript
  annotate("text", x = 37, y = 25, label = as.expression(r_squared_text), 
           hjust = 1.1, vjust = 2, size = 5, color = "black")

#net do plot of standard error exp1 vs standard error exp2 

# Create a data frame from your data
data <- data.frame(orderexp1_se, orderexp2_se)

# Fit the linear model
best <- lm(orderexp2_se ~ orderexp1_se, data = data)

# Plot using ggplot2

# Calculate the R^2 value
r_squared <- summary(best)$r.squared

# Format R^2 text with superscript
r_squared_text <- bquote(R^2 == .(round(r_squared, 2)))

# Plot using ggplot2
ggplot(data, aes(x = orderexp1_se, y = orderexp2_se)) +
  geom_point(size = 2) +  # Plot the points
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "black", linewidth = 1.2) +  # Add the linear model line
  labs(x = "standard error Days to Flower Exp 1", y = " se Days to Flower Exp 2") +
  theme_classic() +
  theme(
    # Adjust axis title margin
    axis.title.x = element_text(margin = margin(t = 10), size = 14),
    axis.title.y = element_text(margin = margin(r = 10), size = 14),
    # Increase axis text size
    axis.text = element_text(size = 14, color = "black")
  ) +
  # Add R^2 annotation with superscript
  annotate("text", x = 0.9, y = 0.2, label = as.expression(r_squared_text), 
           hjust = 1.1, vjust = 2, size = 5, color = "black")

#next do plot of standard error exp1 vs standard error exp2 

#now look a the combined experment and do mean vs s.e 
data <- data.frame(all_mean, all_se)
# Fit the linear model
best <- lm(all_mean ~ all_se, data = data)

# Plot using ggplot2

# Calculate the R^2 value
r_squared <- summary(best)$r.squared

# Format R^2 text with superscript
r_squared_text <- bquote(R^2 == .(round(r_squared, 2)))


# Plot using ggplot2
ggplot(data, aes(x = all_mean, y = all_se)) +
  geom_point(size = 2) +  # Plot the points
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "black", linewidth = 1.2) +  # Add the linear model line
  labs(x = "Mean Days to Flower", y = " st.err Days to Flower") +
  theme_classic() +
  theme(
    # Adjust axis title margin
    axis.title.x = element_text(margin = margin(t = 10), size = 14),
    axis.title.y = element_text(margin = margin(r = 10), size = 14),
    # Increase axis text size
    axis.text = element_text(size = 14, color = "black")
  ) +
  # Add R^2 annotation with superscript
  annotate("text", x = 39, y = 0.5, label = as.expression(r_squared_text), 
           hjust = 1.1, vjust = 2, size = 5, color = "black")


#save combined experiments. 

#write.csv(all_mean, file = "FT_genotype_means_exp1_exp2.csv")
#write.csv(combine, file = "FT_genotype_exp1_exp2.csv")

#write.csv(all_sd, file = "FT_exp1_exp2_sd.csv")
#write.csv(all_se, file = "FT_exp1_exp2_se.csv")


#NUMBER OF LEAVES EXP1 VS EXP2 =============


Leavesonly1 <- data_july[,c(1,3)]
Leavesonly2 <- data_nov[,c(1,3)]

widerdata1<- Leavesonly1 %>%
  group_by(Genotype) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Leaf.no) %>%
  dplyr::select(-row)

widerdata2<- Leavesonly2 %>%
  group_by(Genotype) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Leaf.no) %>%
  dplyr::select(-row)

new_order1 <- sort(colnames(widerdata1))
orderexp1 <- widerdata1[, new_order1]
orderexp1_mean <- apply(orderexp1,2,mean, na.rm = TRUE)
orderexp1_sd <- apply(orderexp1,2,sd, na.rm = TRUE)
orderexp1_se <- apply(orderexp1,2,standard_error)

new_order2 <- sort(colnames(widerdata2))
orderexp2 <- widerdata2[, new_order2]
orderexp2_mean <- apply(orderexp2,2,mean, na.rm = TRUE)
orderexp2_sd <- apply(orderexp2,2,sd, na.rm = TRUE)
orderexp1_se <- apply(orderexp1,2,standard_error)


combine <- rbind(orderexp1,orderexp2)
all_mean <- apply(combine,2,mean, na.rm = TRUE)
all_mean <- as.data.frame(all_mean)
all_sd <- apply(combine,2,sd, na.rm = TRUE)
#ccombine the sd vales


# Apply the function to each column of the dataframe
all_se <- apply(combine, 2, standard_error)
all_se <- as.data.frame(all_se)

all_mean_Leaves <- all_mean
all_se_Leaves <- all_se

#plot of first exp vs second exp means

best <- lm(orderexp2_mean ~ orderexp1_mean)

plot(orderexp1_mean,orderexp2_mean, xlab = "Number of leaves Exp 1", ylab = "Number of leaves Exp 2", lwd = 2, pch = 19)
abline(best, lwd = 2)

#plot of first exp vs second exp standard deviaTION

best <- lm(orderexp2_se ~ orderexp1_se)

plot(orderexp1_se,orderexp2_se, xlab = "Days to Flower Exp 1", ylab = "Days to Flower Exp 2", lwd = 2, pch = 19)
abline(best, lwd = 2)

bestsd <- lm(orderexp2_sd ~ orderexp1_sd)
summary(bestsd)
plot(orderexp1_sd,orderexp2_sd, xlab = "SD of Days to Flower Exp 1", ylab = "SD of Days to Flower Exp 2", lwd = 2, pch = 19, ylim = c(0,5))
abline(bestsd, lwd = 2)
abline(v=20)

boxplot(orderexp1_sd,orderexp2_sd, xlab = "SD of Days to Flower Exp 1", ylab = "SD of Days to Flower Exp 2", lwd = 2, pch = 19)

#save files 

#write.csv(all_mean, file = "leaves_genotype_means_exp1_exp2.csv")
#write.csv(combine, file = "leaves_genotype_exp1_exp2.csv")

#write.csv(all_sd, file = "leaves_exp1_exp2_sd.csv")
#write.csv(all_se, file = "leaves_exp1_exp2_se.csv")

##LEAVES vs FT ================================

colnames(all_mean_FT) <- "all_mean_FT"

data <- data.frame(all_mean_FT, all_mean)
# Fit the linear model
best <- lm(all_mean_FT ~ all_mean, data = data)

# Plot using ggplot2

# Calculate the R^2 value
r_squared <- summary(best)$r.squared

# Format R^2 text with superscript
r_squared_text <- bquote(R^2 == .(round(r_squared, 2)))

# Plot using ggplot2
plot1 <- ggplot(data, aes(x = all_mean_FT, y = all_mean)) +
  geom_point(size = 2) +  # Plot the points
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "black", linewidth = 1.2) +  # Add the linear model line
  labs(x = "Mean Days to Flower", y = "Mean Number of Leaves") +
  theme_classic() +
  theme(
    # Adjust axis title margin
    axis.title.x = element_text(margin = margin(t = 10), size = 14),
    axis.title.y = element_text(margin = margin(r = 10), size = 14),
    # Increase axis text size
    axis.text = element_text(size = 14, color = "black")
  ) +
  # Add R^2 annotation with superscript
  annotate("text", x = 40, y = 12, label = as.expression(r_squared_text), 
           hjust = 1.1, vjust = 2, size = 5, color = "black")


#WEIGHTS EXP1 VS EXP2 ===================================


Biomassonly1 <- data_july[,c(1,4)]
Biomassonly2 <- data_nov[,c(1,4)]

widerdata1<- Biomassonly1 %>%
  group_by(Genotype) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = weight) %>%
  dplyr::select(-row)

widerdata2<- Biomassonly2 %>%
  group_by(Genotype) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = weight) %>%
  dplyr::select(-row)

#then order new datafrmaes by column then can compare between 2 experiments easier

new_order1 <- sort(colnames(widerdata1))
orderexp1 <- widerdata1[, new_order1]
orderexp1_mean <- apply(orderexp1,2,mean, na.rm = TRUE)
orderexp1_sd <- apply(orderexp1,2,sd, na.rm = TRUE)

new_order2 <- sort(colnames(widerdata2))
orderexp2 <- widerdata2[, new_order2]
orderexp2_mean <- apply(orderexp2,2,mean, na.rm = TRUE)
orderexp2_sd <- apply(orderexp2,2,sd, na.rm = TRUE)


combine <- rbind(orderexp1,orderexp2)
all_mean <- apply(combine,2,mean, na.rm = TRUE)
all_mean <- as.data.frame(all_mean)
all_sd <- apply(combine,2,sd, na.rm = TRUE)
#ccombine the sd vales

# Apply the function to each column of the dataframe
all_se <- apply(combine, 2, standard_error)

all_se <- as.data.frame(all_se)

all_mean_Biomass <- all_mean
all_se_Biomass <- all_se

#plot of first exp vs second exp means
best <- lm(orderexp2_mean ~ orderexp1_mean)

plot(orderexp1_mean,orderexp2_mean, xlab = "Dry biomass Exp 1", ylab = "Dry biomass Exp 2", lwd = 2, pch = 19)
abline(best, lwd = 2)
segments(0,0,1,1, lwd =2)
#plot of first exp vs second exp standard deviaTION


bestsd <- lm(orderexp2_sd ~ orderexp1_sd)
summary(bestsd)
plot(orderexp1_sd,orderexp2_sd, xlab = "SD of Days to Flower Exp 1", ylab = "SD of Days to Flower Exp 2", lwd = 2, pch = 19)
abline(bestsd, lwd = 2)
abline(v=20)

boxplot(orderexp1_sd,orderexp2_sd, xlab = "SD of Days to Flower Exp 1", ylab = "SD of Days to Flower Exp 2", lwd = 2, pch = 19)


#Weight vs FT =======================

data <- data.frame(all_mean_FT, all_mean)
# Fit the linear model
best <- lm(all_mean_FT ~ all_mean, data = data)

# Plot using ggplot2

# Calculate the R^2 value
r_squared <- summary(best)$r.squared

# Format R^2 text with superscript
r_squared_text <- bquote(R^2 == .(round(r_squared, 2)))

# Plot using ggplot2
plot2 <- ggplot(data, aes(x = all_mean_FT, y = all_mean)) +
  geom_point(size = 2) +  # Plot the points
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "black", linewidth = 1.2) +  # Add the linear model line
  labs(x = "Mean Days to Flower", y = "Mean Biomass") +
  theme_classic() +
  theme(
    # Adjust axis title margin
    axis.title.x = element_text(margin = margin(t = 10), size = 14),
    axis.title.y = element_text(margin = margin(r = 10), size = 14),
    # Increase axis text size
    axis.text = element_text(size = 14, color = "black")
  ) +
  # Add R^2 annotation with superscript
  annotate("text", x = 40, y = 0.01, label = as.expression(r_squared_text), 
           hjust = 1.1, vjust = 2, size = 5, color = "black")


(plot1|plot2)

#save ===
#write.csv(all_mean, file = "weight_genotype_means_exp1_exp2.csv")
#write.csv(combine, file = "weight_genotype_exp1_exp2.csv")

#write.csv(all_sd, file = "weight_exp1_exp2_sd.csv")
#write.csv(all_se, file = "weight_exp1_exp2_se.csv")


#combine FT, Leaves and Biomass and save for further analysis. 


all_data_combined <- cbind(all_mean_FT, all_se_FT,all_mean_Leaves,all_se_Leaves, all_mean_Biomass, all_se_Biomass)
colnames(all_data_combined) <- c("FT_mean",	"FT_se",	"leaves_mean",	"leaves_se",	"weights_mean",	"weights_se")

#save 
write.csv(all_data_combined, "FT_exp_combined.csv")

#
