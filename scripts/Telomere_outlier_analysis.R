rm(list=ls())
library(tidyverse)
library(ggplot2)
#reading in genomics data
Genomics_data_recoded <- readRDS('../extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded

#any association between gender and telomere length?
#filter out those with no genetic sex
gen2 <- gen2 %>% filter(!is.na(gen2$GeneticSex.0.0))
gen2$GeneticSex.0.0 <- ifelse(gen2$GeneticSex.0.0 == 'Female', 0, 1)

ggplot(data=gen2, aes(x=AdjTSRatio.0.0)) + geom_density() + facet_grid('GeneticSex.0.0')
plot(gen2$GeneticSex.0.0)

fit <- glm(GeneticSex.0.0~AdjTSRatio.0.0, family = binomial(link=logit), data = gen2)
summary(fit)
#no significant difference

#reading in normal biomarker data
data <- readRDS('../extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds')

data_full <- cbind(data, Genomics_data_recoded$AdjTSRatio.0.0)

data_full$TSR <- data_full$`Genomics_data_recoded$AdjTSRatio.0.0`

data_full <- data_full %>% select(-`Genomics_data_recoded$AdjTSRatio.0.0`)

#filtering rows that have more than 20% missingness
data_fil <- data_full %>% filter(rowSums(is.na(data_full))<= 0.2*(ncol(data_full)))

#filtering out those with no telomere length
data_fil <- data_fil %>% filter(!is.na(TSR))

##filtering out columns with more than 20% NA

data_fin <- data_fil %>% select(-contains("LPA"))
data_fin <- data_fin %>% select(-contains("E2"))
data_fin <- data_fin %>% select(-contains("RF"))
data_fin <- data_fin %>% select(-contains("MiAlbUr")

#filtering the data to be all above the 1.5* quantile to see if there are any associations with certain biomarkers
###data_fin <- data_fin %>% filter(!is.na(data_fin$TSR))
Q <- quantile(data_fin$TSR, probs=c(.25, .75))
                                
iqr <- IQR(data_fin$TSR)
data_fin2 <- data_fin %>% filter(TSR > (Q[2] + 1.5*iqr))  
                                
plot(density(data_fin$TSR))

#filtering based on outlier status

#1.5* outlier
data_full_2 <- data_fin %>% filter(TSR > (Q[2] + 1.5*iqr))

pheatmap(cor((data_full_2), use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)

plot(density(data_full_2$TSR))

#not outlier
data_full_3 <- data_fin %>% filter(TSR < (Q[2] + 1.5*iqr))

pheatmap(cor((data_full_3), use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)

plot(density(data_full_3$TSR))

#3* outlier
data_full_4 <- data_fin %>% filter(TSR > (Q[2] + 3*iqr))

pheatmap(cor((data_full_4), use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)

plot(density(data_full_4$TSR))

#6* outlier
data_full_5 <- data_fin %>% filter(TSR > (Q[2] + 6*iqr))

pheatmap(cor((data_full_5), use = 'pairwise.complete.obs'),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = TRUE, show_colnames = TRUE)

plot(density(data_full_5$TSR))

#reading


