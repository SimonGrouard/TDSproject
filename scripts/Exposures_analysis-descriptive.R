
rm(list=ls())

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

full_df <- data.frame(Exposures_covariates_recoded) 
gnm <- data.frame(Genomics_data_recoded)
full_df$TL<- gnm$AdjTSRatio.0.0

# 1. Descriptive analyses of NA ---------------------------------------------------

## 1.1 NA numbers summary and barplots

### 1.1.1 the full exposure dataset
NA_Exposure <- as.data.frame(colSums(is.na(full_df)))
NA_Exposure <- as.data.frame(t(NA_Exposure))
rownames(NA_Exposure) <- c("NAnumber")
pdf("NA_Exposure.pdf", 50, 30)
barplot(as.matrix(NA_Exposure), names.arg = names(NA_Exposure), las = 2,
        cex.axis = 0.6, cex.names = 0.1,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "variables with NAs")
dev.off()

### 1.1.2 variables with n(NA) < 400,000 (about 80% of the total sample size)
tNA_Exposure <- as.data.frame(t(NA_Exposure))
df_lessNA_var <- as.data.frame(subset(tNA_Exposure, NAnumber<= 400000)) 
# only 216 of the 817 variables remained in our dataset
df_lessNA <- full_df %>% select(colnames = rownames(df_lessNA_var))
colnames(df_lessNA) <- rownames(df_lessNA_var)
lessNA_Exposure <- as.data.frame(colSums(is.na(df_lessNA)))
lessNA_Exposure <- t(lessNA_Exposure)
rownames(lessNA_Exposure) <- c("NA_counts")
pdf("lessNA_Exposure.pdf", 30, 16)
barplot(as.matrix(lessNA_Exposure), names.arg = names(lessNA_Exposure), las = 2,
        cex.axis = 0.6, cex.names = 0.2,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "variables with NAs")
dev.off()
# just an explanation (further have to see the new dataset re-created via Step 4)

### 1.1.3 variables with n(NA) < 300,000
tNA_Exposure <- as.data.frame(t(NA_Exposure))
df_lessNA_var <- as.data.frame(subset(tNA_Exposure, NAnumber<= 300000)) 
# 146 of the 817 variables remained
library(dplyr)
df_lessNA <- full_df %>% select(colnames = rownames(df_lessNA_var))
colnames(df_lessNA) <- rownames(df_lessNA_var)
lessNA_Exposure <- as.data.frame(colSums(is.na(df_lessNA)))
lessNA_Exposure <- t(lessNA_Exposure)
rownames(lessNA_Exposure) <- c("NA_counts")
pdf("lessNA_Exposure.pdf", 30, 16)
barplot(as.matrix(lessNA_Exposure), names.arg = names(lessNA_Exposure), las = 2,
        cex.axis = 0.6, cex.names = 0.2,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "variables with NAs")
dev.off()

### 1.1.4 variables with n(NA) < 200,000 
tNA_Exposure <- as.data.frame(t(NA_Exposure))
df_lessNA_var <- as.data.frame(subset(tNA_Exposure, NAnumber<= 200000)) 
# 128 of the 817 variables remained
library(dplyr)
df_lessNA <- full_df %>% select(colnames = rownames(df_lessNA_var))
colnames(df_lessNA) <- rownames(df_lessNA_var)
lessNA_Exposure <- as.data.frame(colSums(is.na(df_lessNA)))
lessNA_Exposure <- t(lessNA_Exposure)
rownames(lessNA_Exposure) <- c("NA_counts")
pdf("lessNA_Exposure.pdf", 30, 16)
barplot(as.matrix(lessNA_Exposure), names.arg = names(lessNA_Exposure), las = 2,
        cex.axis = 1, cex.names = 1,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "variables with NAs")
dev.off()

### 1.1.5 variables with n(NA) < 100,000 
tNA_Exposure <- as.data.frame(t(NA_Exposure))
df_lessNA_var <- as.data.frame(subset(tNA_Exposure, NAnumber<= 100000)) 
# 122 variables remained
library(dplyr)
df_lessNA <- full_df %>% select(colnames = rownames(df_lessNA_var))
colnames(df_lessNA) <- rownames(df_lessNA_var)
lessNA_Exposure <- as.data.frame(colSums(is.na(df_lessNA)))
lessNA_Exposure <- t(lessNA_Exposure)
rownames(lessNA_Exposure) <- c("NA_counts")
pdf("lessNA_Exposure_100,000.pdf", 30,16)
barplot(as.matrix(lessNA_Exposure), names.arg = names(lessNA_Exposure), las = 2,
        cex.axis = 1, cex.names = 0.6,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "variables with NAs")
dev.off()

## NA pattern visualisation
library(naniar)
pdf("NA_Exposure_pattern.pdf", 30, 21)
vis_miss(df_lessNA, cluster= T, sort_miss = T, warn_large_data = F)
dev.off()

# Final dataset 
mydf <- na.omit(df_lessNA) 
dim(mydf) # new sample size = 315,024
summary(mydf)
mydf_type<-as.data.frame(apply(mydf,2,typeof))

# mydf[,1:75]:categorical
mydf[,1:75] <- apply(mydf[,1:75], 2, as.factor)

# mydf[,76:122]:numeric:
## replace 'Prefer not to answer' and 'Do not know' with -2
## replace 'Less than one' with -1
mydf[,76:122][mydf[,76:122]== "Less than one"] <- "-1"
mydf[,76:122][mydf[,76:122] == "Prefer not to answer"] <- "-2"
mydf[,76:122][mydf[,76:122] == "Do not know"] <- "-2"
mydf[,76:122] <- apply(mydf[,76:122], 2, as.numeric)


# 2. Classic descriptive analysis --------------------------------------------
# key covariates(sex, age, ethnicity, bmi, smoke, alcohol, stress, etc.)
# sex
summary_sex <- summary(mydf$Sex.0.0)
tab_sex <- as.data.frame(table(mydf$Sex.0.0))
colnames(tab_sex) <- c("Sex", "Frequency")

perc_sex <- paste0(round(tab_sex$Frequency/sum(tab_sex$Frequency)*100, 2),"%")
lab_sex <- paste0(tab_sex$Frequency, " (", perc_sex, ")")

pdf("Dsc_Sex.pdf", 10, 6)
ggplot(tab_sex, mapping = aes(Sex, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5)) +
  geom_text(aes(label = lab_sex, vjust = 0))
dev.off()

# Ethnicity
summary_ethnicity <- summary(full_df$full_df$Ethnicity.0.0)
tab_ethnicity <- as.data.frame(table(full_df$Ethnicity.0.0))
colnames(tab_ethnicity) <- c("Ethnicity_background", "Frequency")

perc_ethnicity <- paste0(round(tab_ethnicity$Frequency/sum(tab_ethnicity$Frequency)*100, 2),"%")
lab_ethnicity <- paste0(tab_ethnicity$Frequency, " (", perc_ethnicity, ")")

pdf("Dsc_Ethnicity.pdf", 10,6)
ggplot(tab_ethnicity, mapping = aes(x = Ethnicity_background, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  geom_text(aes(label = lab_ethnicity), size = 1.6, vjust = 0)
dev.off()

# smokePast
summary_smokePast <- summary(mydf$SmokePast.0.0)
tab_smokePast <- as.data.frame(table(mydf$SmokePast.0.0))
colnames(tab_smokePast) <- c("Past_Smoking_Status", "Frequency")

perc_smokePast <- paste0(round(tab_smokePast$Frequency/sum(tab_smokePast$Frequency)*100, 2),"%")
lab_smokePast <- paste0(tab_smokePast$Frequency, " (", perc_smokePast, ")")

pdf("Dsc_SmokePast.pdf", 10, 6)
ggplot(tab_smokePast, mapping = aes(Past_Smoking_Status, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5)) +
  geom_text(aes(label = lab_smokePast), vjust = 0)
dev.off()

# smokeCurr
summary_smokeCurr <- summary(mydf$SmokeCurr.0.0)
tab_smokeCurr <- as.data.frame(table(mydf$SmokeCurr.0.0))
colnames(tab_smokeCurr) <- c("Current_Smoking_Status", "Frequency")

perc_smokeCurr <- paste0(round(tab_smokeCurr$Frequency/sum(tab_smokeCurr$Frequency)*100, 2),"%")
lab_smokeCurr <- paste0(tab_smokeCurr$Frequency, " (", perc_smokeCurr, ")")

pdf("Dsc_SmokeCurr.pdf", 10, 6)
ggplot(tab_smokeCurr, mapping = aes(Current_Smoking_Status, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5)) +
  geom_text(aes(label = lab_smokeCurr), vjust = 0)
dev.off()

# age
summary_age <- summary(mydf$AgeAssess.0.0)
pdf("Dsc_Age.pdf", 10, 6)
ggplot(as.data.frame(mydf$AgeAssess.0.0), aes(x = mydf$AgeAssess.0.0)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4) +
  labs(title = "Age", x = "Age/years", y = "density") +
  xlim(c(35, 75))
dev.off()

# bmi
summary_bmi <- summary(mydf$BMI.0.0)
summary_age <- summary(mydf$AgeAssess.0.0)
pdf("Dsc_BMI.pdf", 10, 6)
ggplot(as.data.frame(mydf$BMI.0.0), aes(x = mydf$BMI.0.0)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4) +
  labs(title = "BMI", x = "BMI/kg·m^2", y = "density") +
  xlim(c(10, 50))
dev.off()

# alcohol intake frequency
summary_AlcFreq <- summary(mydf$AlcFreq.0.0)
tab_AlcFreq <- as.data.frame(table(mydf$AlcFreq.0.0))
colnames(tab_AlcFreq) <- c("Alcohol_Intake_Frequency", "Frequency")
perc_AlcFreq <- paste0(round(tab_AlcFreq$Frequency/sum(tab_AlcFreq$Frequency)*100, 2),"%")
lab_AlcFreq <- paste0(tab_AlcFreq$Frequency, " (", perc_AlcFreq, ")")

pdf("Dsc_AlcFreq.pdf", 10,6)
ggplot(tab_AlcFreq, mapping = aes(x = Alcohol_Intake_Frequency, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  geom_text(aes(label = lab_AlcFreq), size =3, vjust = 0)
dev.off()

# Qualifications
summary_qlf<- summary(mydf$Qualifications.0.0)
tab_qlf <- as.data.frame(table(mydf$Qualifications.0.0))
colnames(tab_qlf) <- c("Qualifications", "Frequency")
perc_qlf <- paste0(round(tab_qlf$Frequency/sum(tab_qlf$Frequency)*100, 2),"%")
lab_qlf <- paste0(tab_qlf$Frequency, " (", perc_qlf, ")")

pdf("Dsc_Qualifications.pdf", 10,6)
ggplot(tab_qlf, mapping = aes(x = Qualifications, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5)) + 
  geom_text(aes(label = lab_qlf), size =3, vjust = 0)
dev.off()

# Stress
summary_MS2Yrs <- summary(mydf$IllMournStress2Yrs.0.0)
tab_MS2Yrs <- as.data.frame(table(mydf$IllMournStress2Yrs.0.0))
colnames(tab_MS2Yrs) <- c("Illness_injury_bereavement_stress_in_last_2_years", "Frequency")
perc_MS2Yrs <- paste0(round(tab_MS2Yrs$Frequency/sum(tab_MS2Yrs$Frequency)*100, 2),"%")
lab_MS2Yrs <- paste0(tab_MS2Yrs$Frequency, " (", perc_MS2Yrs, ")")

pdf("Dsc_MS2Yrs.pdf", 10, 6)
ggplot(tab_MS2Yrs, mapping = aes(x = Illness_injury_bereavement_stress_in_last_2_years, y = Frequency)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  geom_text(aes(label = lab_MS2Yrs), size =3, vjust = 0)
dev.off()

# 3. Correlation heatmap for continuous variables--------------------------------------------------
suppressPackageStartupMessages(library(pheatmap))

scaledMydf <- as.data.frame(scale(mydf[76:122]))

breaks = seq(-1, 1, length.out = 100)
pdf("corHmpNumeric.pdf", 15, 15)
pheatmap(cor(scaledMydf), breaks =breaks,
         show_rownames = T, show_colnames = T)
dev.off()

# 4. MCA --------------------------------------------------------------------
install.packages(c("FactoMineR", "factoextra"))
suppressPackageStartupMessages(library("FactoMineR", "factoextra"))
res.mca <- MCA(mydf[1:100,1:75], ncp = 5, graph = TRUE)
print(res.mca)

## Eigenvalues/Variances
library("factoextra")
eig.val <- get_eigenvalue(res.mca)
pdf("screeplot.pdf", 10, 6)
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))
dev.off()

## Biplot
pdf("Biplot.pdf", 10, 6)
fviz_mca_biplot(res.mca, 
                repel = TRUE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal())
dev.off()

## Result
var <- get_mca_var(res.mca)
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

## Correlation between variables and principal dimensions
pdf("Corrplot.pdf", 10, 6)
fviz_mca_var(res.mca, choice = "mca.cor", 
             repel = TRUE, # Avoid text overlapping (slow)
             ggtheme = theme_minimal())
dev.off()

## Coordinates of variable categories
head(round(var$coord, 2), 4)
pdf("VarCatPlot.pdf", 10, 6)
fviz_mca_var(res.mca, 
             repel = TRUE, # Avoid text overlapping (slow)
             ggtheme = theme_minimal()) # visualise only variable categories
dev.off()

## Quality of representation of variable categories
head(var$cos2, 4)
# Color by cos2 values: quality on the factor map
pdf("VarCatCos2Value.pdf", 10, 6)
fviz_mca_var(res.mca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal())
dev.off()
# Cos2 Corrplot
library("corrplot")
pdf("Cos2CorrPlot.pdf", 24, 16)
corrplot(var$cos2, is.corr=FALSE)
dev.off()

# Cos2 of variable categories on Dim 1 and 2
pdf("VarCatCos2ValueBarplot.pdf", 45, 27)
fviz_cos2(res.mca, choice = "var", axes = 1:2)
dev.off()

### to be follow-uped







