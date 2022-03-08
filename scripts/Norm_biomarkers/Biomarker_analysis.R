### THE OVERALL FILE I HAVE USED FOR CREATING THE APPROPRIATE BITS OF CODE TO BE RUN ON THE SERVER IN ANOTHER FILE
#AND GRAPH CREATION
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
install.packages('glmnet')
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(tidyverse))

rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
 
#--------------------------------------------------------
#reading in data and combining the telomere length (twice for the different separations of data later in the script)
data_1 <- readRDS('../extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds')
data_1 <- data.frame(data_1)

Genomics_data_recoded <- readRDS('../../extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded

data_1$tlen <- gen$AdjTSRatio.0.0

data <- readRDS('../../extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds')
data <- data.frame(data)

data$tlen <- gen$AdjTSRatio.0.0

#--------------------------------------------------------------
#univariate analysis and plotting
univ_b <- readRDS('../Results_nightingale_biomarker_univ/univ_pval_biomarker_1.rds')

univ_b <- data.frame(univ_b)
univ_b$pval <- as.numeric(levels(univ_b$pval))[univ_b$pval]

univ_b_2 <- univ_b[-1,]

plot(-log10(univ_b_2$pval),pch=16,ylab="-log(p-values)")
abline(h=-log10(0.05/ncol(univ_b_2)), col='red')

table(univ_b_2$pval<=(0.05))
#46 significant at 0.05

univ_b_2$pvalBF <- p.adjust(univ_b_2$pval, method="bonf")
table(univ_b_2$pvalBF<=(0.05))
#45 significant for BF 0.05

univ_b_2$pvalBH <- p.adjust(univ_b_2$pval, method="BH")
table(univ_b_2$pval<=(0.05))
#46 significant for BH

#which are most significant

univ_sig_order <- univ_b_2[order(univ_b_2$pval), ]

plot(-log10(univ_2$pval),pch=16,ylab="-log(p-values)")
abline(h=-log10(0.05/ncol(univ)), col='red')
abline(h=-log10(0.05), col='blue')

#------------------------------------
#descriptive analysis and filtering of NA values

#assessment of NA values across rows
table(rowSums(is.na(data))>10)
#when setting the threshold at ~20% missingness results in loss of 14.7% of participants

install.packages('pheatmap')
suppressPackageStartupMessages(library(pheatmap))

#heatmap of NAs across the variables
pheatmap(cor(is.na(data),use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)

#heatmap of all correlations
pheatmap(cor(data, use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)


#removing highly correlated variables where includes both percentage and count of blood cells
data <- data %>% select(-contains("Perc"))

#plotting heatmap after removing highly correlated variables
pheatmap(cor(is.na(data_new)),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)

summary((data))
# looking at the summary there are some variables with very large numbers of NAs so further analysis should be performed with 
#a fraction of these filtered out

#filtering rows that have more than 20% missingness
data_fil <- data %>% filter(rowSums(is.na(data))<= 0.2*(ncol(data)))

#filtering out those with no telomere length
data_fil <- data_fil %>% filter(!is.na(tlen))

#filtering columns that have more than 20% missingness
#data_fil_3 <- data_fil_2 %>% filter(colSums(is.na(data_fil_2))<= 0.2*(nrow(data_fil_2))) #for some reason this line of code isnt working not sure why so will filter manually below
#running colSums(is.na(data_fil_2))<= 0.2*(nrow(data_fil_2)) in the terminal reveals that there are four variables that have more than 20% missingness

##filtering out columns with more than 20% NA

df_fin <- data_fil %>% select(-contains("LPA"))
df_fin <- df_fin %>% select(-contains("E2"))
df_fin <- df_fin %>% select(-contains("RF"))
df_fin <- df_fin %>% select(-contains("MiAlbUr"))

table(is.na(df_fin))

#------------------------------------------------------
#descriptive analysis

pheatmap(cor((df_fin), use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)
#can see strong positive correlations between APOB, cholesterol, LDLDir
#strong pos correlation between APOA and HDL Cholesterol
#see strong neg correlation between lymphocytes and neutrophils
#pos correlation between AST and ALT (aspartate and alanine transaminase)

####checking if correlations change with percentage instead of count
data_1 <- data_1 %>% select(-contains("Count"))
data_fil_1 <- data_1 %>% filter(rowSums(is.na(data_1))<= 0.2*(ncol(data_1)))


pheatmap(cor((df_fin_1), use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)


##running PCA for descriptive analysis
df_om <- na.omit(df_fin)
pcaX = prcomp(df_om)
ev = with(pcaX, sdev**2/sum(sdev**2))
plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
     ylab = "Proportion of EV", ylim = c(0, 1.2), cex = 0.3)
points(cumsum(ev), pch = 19, col = "skyblue", cex = 0.3)
legend("top", pch = 19, col = c("navy", "skyblue"),
       legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)

sum(!cumsum(ev)>0.9)

#assessment of neutrophil to lymphocyte ratio---------------------------------------------

summary(df_fin$LymphCount.0.0)
summary(df_fin$NeutCount.0.0)

set.seed(1)
fin_sample <-df_fin[sample(1:nrow(df_fin), 100000), ]


vals <- fin_sample$NeutCount.0.0/fin_sample$LymphCount.0.0

plot(density(vals, na.rm=TRUE)) 

summary(vals)
#from this small sample can see the max value is 29.518- massive outlier and maybe a technical error

boxplot(vals)

summary(fin_sample$NeutCount.0.0)

#are these associations caused by cancer 

#recoding the cancer variable to yes and no

Exposures_covariates_recoded_combined_final <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds")

Exposures_covariates_recoded_combined_final$CancerDiag <- ifelse(Exposures_covariates_recoded_combined_final$CancerDiag == 'Yes - you will be asked about this later by an interviewer', 1, 0)

comb <- cbind(data, Exposures_covariates_recoded_combined_final$CancerDiag)

comb <- comb %>% select(-contains("Perc"))

#filtering rows that have more than 20% missingness
comb_fil <- comb %>% filter(rowSums(is.na(comb))<= 0.2*(ncol(comb)))

#filtering out those with no telomere length
comb_fil <- comb_fil %>% filter(!is.na(tlen))

##filtering out columns with more than 20% NA

comb_fin <- comb_fil %>% select(-contains("LPA"))
comb_fin <- comb_fin %>% select(-contains("E2"))
comb_fin <- comb_fin %>% select(-contains("RF"))
comb_fin <- comb_fin %>% select(-contains("MiAlbUr"))

#analysis for counts
comb <- comb %>% select(-contains("Count"))

set.seed(1)
comb_sample <-comb_fin[sample(1:nrow(comb_fin), 100000), ]

boxplot(comb_fin$LymphPerc.0.0~comb_fin$`Exposures_covariates_recoded_combined_final$CancerDiag`, xlab = 'Cancer status', ylab = 'Lymphocyte Percentage')
boxplot(comb_fin$NeuPerc.0.0~comb_fin$`Exposures_covariates_recoded_combined_final$CancerDiag`, xlab = 'Cancer status', ylab = 'Neutrophil Percentage')

plot(comb_fin$LymphCount.0.0~comb_fin$NeutCount.0.0, xlab= 'Neutrophil count', ylab='Lymphocyte count')

comb_sample %>%
  ggplot(aes(x=LymphCount.0.0,y=NeutCount.0.0, color=`Exposures_covariates_recoded_combined_final$CancerDiag`)) +
  geom_point(alpha=0.5, size=2) +
  labs(y='Neutrophil count', x='Lymphocyte count')

#analysis for perc
boxplot(comb_fin$LymphCount.0.0~comb_fin$`Exposures_covariates_recoded_combined_final$CancerDiag`, xlab = 'Cancer status', ylab = 'Lymphocyte count')
boxplot(comb_fin$NeutCount.0.0~comb_fin$`Exposures_covariates_recoded_combined_final$CancerDiag`, xlab = 'Cancer status', ylab = 'Neutrophil count')

plot(comb_fin$LymphPerc.0.0~comb_fin$NeuPerc.0.0, xlab= 'Neutrophil percentage', ylab='Lymphocyte percentage', col = comb_fin$`Exposures_covariates_recoded_combined_final$CancerDiag`)

library(ggplot2)

set.seed(2)
comb_sample <-comb_fin[sample(1:nrow(comb_fin), 100000), ]

comb_sample %>%
  ggplot(aes(LymphPerc.0.0,NeuPerc.0.0, color=`Exposures_covariates_recoded_combined_final$CancerDiag`)) +
  geom_point(alpha=0.5, size=2) +
  labs(y='Neutrophil percentage', x='Lymphocyte percentage')

##Adjusted univariate analysis----------------------

adj <- readRDS('../Results_nightingale_biomarker_univ/Results_biomarker_univ_3/univ_pval_adj_biomarker.rds')
adj <- data.frame(adj)
adj <- adj[c(-53),]
adj <- adj[c(-52),]
adj <- adj[c(-51),]

adj$pval <- as.numeric(levels(adj$pval))[adj$pval]

plot(-log10(adj$pval),pch=16,ylab="-log(p-values)")
abline(h=-log10(0.05/ncol(adj)), col='red')
abline(h=-log10(0.05), col='blue')


table(adj$pval<=(0.05))
#39 significant at 0.05

adj$pvalBF <- p.adjust(adj$pval, method="bonf")
table(adj$pvalBF<=(0.05))
#33 significant for BF 0.05

adj$pvalBH <- p.adjust(adj$pval, method="BH")
table(adj$pval<=(0.05))
#39 significant for BH

#which are most significant
adj_order <- adj[order(adj$pval), ]

saveRDS(adj_order, 'Biomarker_adj_df.rds')


#Adjusted count data analysis----------------------------

adj_count <- readRDS('../../Results_nightingale_biomarker_univ/Results_biomarker_univ_3/univ_pval_adj_biomarker_count.rds')

adj_count <- data.frame(adj_count)

adj_count$pval <- as.numeric(levels(adj_count$pval))[adj_count$pval]

adj_count <- adj_count[c(-46),]
adj_count <- adj_count[c(-45),]
adj_count <- adj_count[c(-44),]

plot(-log10(adj_count$pval),pch=16,ylab="-log(p-values)")
abline(h=-log10(0.05/ncol(adj_count)), col='red')
abline(h=-log10(0.05), col='blue')


table(adj_count$pval<=(0.05))
#32 significant at 0.05

adj_count$pvalBF <- p.adjust(adj_count$pval, method="bonf")
table(adj_count$pvalBF<=(0.05))
#26 significant for BF 0.05

adj_count$pvalBH <- p.adjust(adj_count$pval, method="BH")
table(adj_count$pval<=(0.05))
#32 significant for BH

adj_count_order <- adj_count[order(adj_count$pval),] 

#Adjusted percentage data analysis-------------------

adj_perc <- readRDS('../../Results_nightingale_biomarker_univ/Results_biomarker_univ_3/univ_pval_adj_biomarker_perc.rds')

adj_perc <- data.frame(adj_perc)
