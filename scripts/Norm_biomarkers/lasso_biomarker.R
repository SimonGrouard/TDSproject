# Libraries ---------------------------------------------------------------

library(lme4)
library(dplyr)
library(glmnet)

# Data --------------------------------------------------------------------

telomere_length <- readRDS(here::here("extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- readRDS(here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds"))
bio <- readRDS(here::here("extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds"))
test_telomere <- head(telomere_length, 1000)
test_exposures <- head(exposures, 1000)

# Data Cleaning -----------------------------------------------------------
## Adding telomere length variable to exposure dataframe
bio_telomere <- bio %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0,
         ZAdjTSRatio = telomere_length$ZAdjTSRatio.0.0,
         Sex = exposures$Sex,
         Age = exposures$AgeAssess)
test_bio_telomere <- head(bio_telomere, 1000)

## Removal of NAs




# Analysis ----------------------------------------------------------------
## LASSO Models
set.seed(1)
lasso_model <- cv.glmnet(x = test_exposures, 
                         y = test_telomere$AdjTSRatio.0.0)
plot(lasso_model)

# Work in Progress
lassomodels.1se = lassomodels = NULL
for (k in 1:ncol(proteins)) {
  print(k)
  model.lasso = cv.glmnet(x = tr, y = proteins[,
                                               k])
  lassomodels.1se = cbind(lassomodels.1se, coef(model.lasso,
                                                s = "lambda.1se")[-1])
  lassomodels = cbind(lassomodels, coef(model.lasso,
                                        s = "lambda.min")[-1])
}
t1 = Sys.time()
print(t1 - t0)
colnames(lassomodels) = colnames(lassomodels.1se) = colnames(proteins)
rownames(lassomodels) = rownames(lassomodels.1se) = colnames(tr)


