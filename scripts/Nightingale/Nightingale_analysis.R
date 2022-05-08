suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
install.packages('glmnet')
suppressPackageStartupMessages(library(glmnet))
library(tidyverse)

#-----------------
#reading in data
rm(list=ls())

nightingale <- readRDS('../extraction_and_recording/outputs/final/Nightingale_metabolomics_final.rds')
df <- nightingale
Genomics_data_recoded <- readRDS('../extraction_and_recording/outputs/final/Genomics_data_final.rds')

genomics <- data.frame(Genomics_data_recoded)
df$TSR <- genomics$AdjTSRatio.0.0

##filtering out all those who did not take part in the metabolomics study (checked and for the first 10 or so columns it was the exact same number of those with/without NA so omitted those rows specifically from the first column for ease)
df_fil <- df %>%
  filter_at(vars(TotChol.0.0), all_vars(!is.na(.)))  

##filtering out any columns relating to technical covariates
df_new <- df_fil %>% select(-contains("QCF"))
df_new <- df_new %>% select(-contains("Spectr"))
df_new <- df_new %>% select(-contains("MQF"))
df_new <- df_new %>% select(-contains("Ship"))
df_new <- df_new %>% select(-contains("HLactate"))
df_new <- df_new %>% select(-contains("HPyruvate"))
df_new <- df_new %>% select(-contains("LGlucose"))
df_new <- df_new %>% select(-contains("LProtein"))

#-----------------------------------
#reading in univariate results

univ <- readRDS('../Results_nightingale_biomarker_univ/univ_pval_nightingale_1.rds')

univ <- data.frame(univ)
univ_2 <- univ[c(-166),]
univ_2$pval <- as.numeric(levels(univ_2$pval))[univ_2$pval]


plot(-log10(univ_2$pval),pch=16,ylab="-log(p-values)")
abline(h=-log10(0.05/ncol(univ)), col='red')
abline(h=-log10(0.05), col='blue')

table(univ_2$pval<=(0.05))
#145 significant at 0.05

univ_2$pvalBF <- p.adjust(univ_2$pval, method="bonf")
table(univ_2$pvalBF<=(0.05))
#123 significant for BF 0.05

univ_2$pvalBH <- p.adjust(univ_2$pval, method="BH")
table(univ_2$pval<=(0.05))
#145 significant for BH

#which are most significant
univ_night_order <- univ_2[order(univ_2$pval), ]


univ_fil <- univ %>% filter(pvalBH < 0.05)
saveRDS(univ_fil, 'Univ_nightingale_df.rds')

names <- univ_fil[,1] #a vector of the names of biomarkers found to be statistically significant at BH 0.05 level

#filtering the dataframe based on the non significant biomarkers

df_new_2 <- df_new[ ,colnames(df_new) %in% names]

#-------------------------------------------------------------
#NA analysis

install.packages('pheatmap')
suppressPackageStartupMessages(library(pheatmap))

#heatmap of NA values
pheatmap(cor(is.na(df_new_2)),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)

table(is.na(df_new_2$TSR))

#filtering rows that have more than 20% missingness
df_na_fil <- df_new_2 %>% filter(rowSums(is.na(df_new_2))<= 0.2*(ncol(df_new_2)))

#filtering out those with no telomere length
df_tl_fil <- df_na_fil %>% filter(!is.na(TSR))

#filtering columns that have more than 20% missingness
#there are no columns with more than 20% missingness or even 10% missingness so levels are low


#-----------------------------------------------------------
##plotting a heatmap

cor <- cor(df_tl_fil, use = 'pairwise.complete.obs')
saveRDS(cor, 'Nightingale_correlation.RDS')
pheatmap(cor(df_tl_fil, use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)
##can see some clusters of positive correlation mostly

##running PCA for descriptive analysis

#getting rid of all other na values because it cant deal with them
df_om <- na.omit(df_tl_fil)

pcaX = prcomp(df_om)
ev = with(pcaX, sdev**2/sum(sdev**2))
plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
     ylab = "Proportion of EV", ylim = c(0, 1.2), cex = 0.3)
points(cumsum(ev), pch = 19, col = "skyblue", cex = 0.3)
legend("top", pch = 19, col = c("navy", "skyblue"),
       legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)

sum(!cumsum(ev)>0.9)
#only three principal components are needed to explain 90% of the variance

#looking at negative correlations
n_cor <- data.frame(Nightingale_correlation)

n_cor_neg <-  plot(n_cor[,] < 0)
