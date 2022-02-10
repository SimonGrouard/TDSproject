suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
install.packages('glmnet')
suppressPackageStartupMessages(library(glmnet))

df <- data.frame(biomarker_nNightingale_recoded)
genomics <- data.frame(Genomics_data_recoded)
df$tlen <- genomics$AdjTSRatio.0.0
pval <- rep(0,ncol(df)) 
for (j in 1:ncol(df)){
  model <- lm(tlen~.,data=df[,c(j,ncol(df))])
  coefs <- summary(model)$coefficients
  pval[j] <- coefs[2,'Pr(>|t|)']
}

plot(-log10(pval[1:ncol(df)]),pch=16,ylab="-log(p-values)",xlab="TS Ratio") # Manhattan plot
abline(h=-log10(0.05/ncol(df)), col='red',lty=2) # p-val = 0.05/p 
colnames(df[pval<=0.05/ncol(df)])
