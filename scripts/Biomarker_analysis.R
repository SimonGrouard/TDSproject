
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
install.packages('glmnet')
suppressPackageStartupMessages(library(glmnet))


rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#--------------------------------------------------------
#reading in data and combining the telomere length
data <- readRDS('../extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds')
data <- data.frame(data)

Genomics_data_recoded <- readRDS('../extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded

data$tlen <- gen$AdjTSRatio.0.0

#--------------------------------------------------------------
#univariate analysis and plotting
pval <- rep(0,ncol(data)) 
for (j in 1:ncol(data)){
  model <- lm(tlen~.,data=data[,c(j,ncol(data))])
  coefs <- summary(model)$coefficients
  pval[j] <- coefs[2,'Pr(>|t|)']
}

plot(-log10(pval[1:ncol(data)]),pch=16,ylab="-log(p-values)",xlab="TS Ratio") # Manhattan plot
abline(h=-log10(0.05/ncol(data)), col='red',lty=2) # p-val = 0.05/p 
colnames(df[pval<=0.05/ncol(data)])

#using a subset of the data
df <- data[sample(1:nrow(data), 100000),]

pval <- rep(0,ncol(df)) 
for (j in 1:ncol(df)){
model <- lm(tlen~.,data=df[,c(j,ncol(df))])
  coefs <- summary(model)$coefficients
  pval[j] <- coefs[2,'Pr(>|t|)']
}

plot(-log10(pval[1:ncol(df)]),pch=16,ylab="-log(p-values)",xlab="TS Ratio") # Manhattan plot
abline(h=-log10(0.05/ncol(df)), col='red',lty=2) # p-val = 0.05/p 
colnames(df[pval<=0.05/ncol(df)])

#there are 9 more significant biomarkers when using the full dataset

#------------------------------------
#descriptive analysis and filtering of NA values

#assessment of NA values across rows
table(rowSums(is.na(data))>10)
#when setting the threshold at ~20% missingness results in loss of 14.7% of participants

data_fil <- data %>% filter(rowSums(is.na(data))<=10.2)

install.packages('pheatmap')
suppressPackageStartupMessages(library(pheatmap))

#heatmap of NAs across the variables
pheatmap(cor(is.na(data_fil)),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)
#can see clusters of correlation

summary((data))
# looking at the summary there are some variables with very large numbers of NAs so further analysis should be performed with 
#a fraction of these filtered out

#setting a 20% threshold for NA
table(colSums(is.na(data_fil))>100492.2)
#5 columns have more than 20% of the participants with no record so is best to remove them
#LPA has 126859 NA
#E2 has 425806 NA 
#DirBil has 103886 NA
#RF has 461151 NA
#MiAlbUr 349574 NA

##filtering out columns with more than 20% NA
df_new <- data_fil %>% select(-contains("LPA"))
df_new <- df_new %>% select(-contains("E2"))
df_new <- df_new %>% select(-contains("DirBil"))
df_new <- df_new %>% select(-contains("RF"))
df_new <- df_new %>% select(-contains("MiAlbUr"))

table(is.na(df_new))

#------------------------------------------------------
#descriptive analysis

pheatmap(cor((df_new), use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)
#can see strong positive correlations between APOB, cholesterol, LDLDir
#strong pos correlation between APOA and HDL Cholesterol
#see strong neg correlation between lymphocytes and neutrophils
#pos correlation between AST and ALT (aspartate and alanine transaminase)
#obvious correlations between variables that have both a count and percentage metric- perhaps would be best to filter the counts and keep percentages?

df_new <- df_new %>% select(-contains("Count"))



##running PCA for descriptive analysis
pcaX = prcomp(df_new)
ev = with(pcaX, sdev**2/sum(sdev**2))
plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
     ylab = "Proportion of EV", ylim = c(0, 1.2), cex = 0.3)
points(cumsum(ev), pch = 19, col = "skyblue", cex = 0.3)
legend("top", pch = 19, col = c("navy", "skyblue"),
       legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)

sum(!cumsum(ev)>0.9)




