# Libraries ---------------------------------------------------------------

#library(lme4)
library(dplyr)
library(glmnet)

# Data --------------------------------------------------------------------

telomere_length <- readRDS(("../../extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- readRDS(("../../extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds"))
bio <- readRDS(("../../extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds"))
test_telomere <- head(telomere_length, 100)
test_biomarker <- head(bio, 100)

# Data Cleaning -----------------------------------------------------------
## Adding telomere length variable to biomarker dataframe
bio_telomere <- bio %>%
  mutate(Age = exposures$AgeAssess,
         AdjTSRatio = telomere_length$AdjTSRatio.0.0)
#test_bio_telomere <- head(bio_telomere, 1000)

## Removal of NAs

#filtering rows that have more than 20% missingness
#bio_fil <- test_bio_telomere %>% filter(rowSums(is.na(test_bio_telomere))<= 0.2*(ncol(test_bio_telomere)))
bio_fil <- bio_telomere %>% filter(rowSums(is.na(bio_telomere))<= 0.2*(ncol(bio_telomere)))


#filtering out those with no telomere length
bio_fil <- bio_fil %>% filter(!is.na(AdjTSRatio))

##filtering out columns with more than 20% NA

bio_fin <- bio_fil %>% select(-contains("LPA"))
bio_fin <- bio_fin %>% select(-contains("E2"))
bio_fin <- bio_fin %>% select(-contains("RF"))
bio_fin <- bio_fin %>% select(-contains("MiAlbUr"))
bio_fin = na.omit(bio_fin)
# Analysis ----------------------------------------------------------------
## LASSO Models
x = as.matrix(bio_fin[,1:47])
y= bio_fin$AdjTSRatio
set.seed(1)
model.lasso <- cv.glmnet(x = x, y = y, alpha = 1)

betas=coef(model.lasso, s='lambda.1se')[-1]
names(betas)=rownames(coef(model.lasso))[-1]

#pdf('lasso_plot.pdf')
#plot(model.lasso)
#dev.off()

best_lam <- model.lasso$lambda.1se
best_model <- glmnet(x = x, 
                     y = y, lambda = best_lam)

table(coef(model.lasso, s='lambda.1se')[-1]!=0)
betas_2=coef(best_model, s='lambda.1se')[-1]
names(betas_2)=rownames(coef(best_model, s='lambda.1se'))[-1]

#pdf('lasso_plot_best.pdf')
#plot(best_model)
#dev.off

#saveRDS(betas, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_bio_pt1_3.rds')
#saveRDS(betas_2, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_bio_pt1_2.rds')

#plotting beta values-----------------------

#'best' lambda beta values
#pdf('Beta_plot.pdf')
#plot(lasso_bio_pt1_2[lasso_bio_pt1_2!=0], type = 'h', col='navy', lwd=2, xaxt='n', xlab='', ylab='Beta')
#axis(side = 1, at = 1:sum(lasso_bio_pt1_2!=0), labels = names(lasso_bio_pt1_2)[lasso_bio_pt1_2!=0], las=2, cex.axis = 0.5)
#abline(h=0, lty=2)
#dev.off()