

rm(list=ls())

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

exps <- data.frame(Exposures_covariates_recoded) 
gnm <- data.frame(Genomics_data_recoded)
exps$TL<- gnm$AdjTSRatio.0.0

# 1. Descriptive analyses of NA ---------------------------------------------------

## 1.1 NA numbers summary and barplots

### 1.1.1 the full exposure dataset
NA_Exposure <- as.data.frame(colSums(is.na(full_df)))
NA_Exposure <- as.data.frame(t(NA_Exposure))
rownames(NA_Exposure) <- c("NAnumber")
pdf("NA_Exposure.pdf", 24, 16)
barplot(as.matrix(NA_Exposure), names.arg = names(NA_Exposure), las = 2,
        cex.axis = 0.6, cex.names = 0.1,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "variables with NAs")
dev.off()

### 1.1.2 have a quick look at variables with n(NA) < 400,000 (about 80% of the total sample size)
tNA_Exposure <- as.data.frame(t(NA_Exposure))
df_lessNA_var <- as.data.frame(subset(tNA_Exposure, NAnumber<= 400000)) 
# only 208 of tthe 1090 variables left (in the dataset created on 09/02/2022)
library(dplyr)
df_lessNA <- full_df %>% select(colnames = rownames(df_lessNA_var))
colnames(df_lessNA) <- rownames(df_lessNA_var)
lessNA_Exposure <- as.data.frame(colSums(is.na(df_lessNA)))
lessNA_Exposure <- t(lessNA_Exposure)
rownames(lessNA_Exposure) <- c("NA_counts")
pdf("lessNA_Exposure.pdf", 24, 16)
barplot(as.matrix(lessNA_Exposure), names.arg = names(lessNA_Exposure), las = 2,
        cex.axis = 0.6, cex.names = 0.2,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "variables with NAs")
dev.off()
# just an explanation (further have to see the new dataset re-created via Step 4)
## NA pattern visualisation
library(naniar)
pdf("NA_Exposure_pattern.pdf", 24, 16)
vis_miss(full_df, cluster= T, sort_miss = T, warn_large_data = F)
dev.off()

# 2. Classic descriptive analysis --------------------------------------------
# TL and key covariates(age, sex, ethnicity, bmi, smoke, alcohol, stress, pollution, etc.)

# Outcome：TL
summary(full_df$TL)
ggplot(as.data.frame(full_df$TL), aes(x = full_df$TL)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4) +
  labs(title = "Telomere length", x = "Telomere length", y = "density") +
  xlim(c(0, 2))

# age
summary(full_df$AgeAssess.0.0)
ggplot(as.data.frame(full_df$BMI.0.0), aes(x = full_df$BMI.0.0)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4) +
  labs(title = "Age", x = "Age/kg·m^2", y = "density") +
  xlim(c(10, 55))

# bmi
summary(full_df$BMI.0.0)
ggplot(as.data.frame(full_df$BMI.0.0), aes(x = full_df$BMI.0.0)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4) +
  labs(title = "BMI", x = "BMI/kg·m^2", y = "density") +
  xlim(c(10, 50))

# alcohol
summary(full_df$ALC)
hist(full_df$ALC, main = "Alcohol", xlab = "Alcohol consumption/g", xlim = c(0,150))
lines(density(full_df$ALC), na.rm = T)
ggplot(as.data.frame(full_df$ALC), aes(x = full_df$ALC)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4) +
  labs(title = "Alcohol", x = "Alcohol consumption/g", y = "density") +
  xlim(c(0, 120)) +
  ylim(c(0, 0.125))

# smoking status
tab_smoke <- as.data.frame(table(full_df$TobSmok))
tab_smoke
colnames(tab_smoke) <- c("Smoking_status", "Frequency")
ggplot(tab_smoke, mapping = aes(Smoking_status, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  geom_text(aes(label = Frequency), vjust = 0)

# Ethnicity
tab_ethnicity <- as.data.frame(table(full_df$Ethnicity.0.0))
colnames(tab_ethnicity) <- c("Ethnicity_background", "Frequency")
ggplot(tab_ethnicity, mapping = aes(x = Ethnicity_background, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  geom_text(aes(label = Frequency), vjust = 0)

# qualification

# stress

# pollution 

# Heatmap showing correlation among exposures -----------------------------
suppressPackageStartupMessages(library(pheatmap))
# breaks = seq(, , length.out = ) # to be follow-uped
pheatmap(cor(exps), breaks = breaks,
         show_rownames = FALSE, show_colnames = FALSE)


# Univariate analyses -----------------------------------------------------


