
# Libraries ---------------------------------------------------------------

suppressMessages(library(lme4))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))


# parameters --------------------------------------------------------------

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])

# Data --------------------------------------------------------------------

telomere_length <- readRDS(here::here("extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- readRDS(here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds"))
test_telomere <- head(telomere_length, 1000)
test_exposures <- head(exposures, 1000)

# Data Cleaning -----------------------------------------------------------
## Adding telomere length variable to exposure dataframe
exposures_telomere <- exposures %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0,
         ZAdjTSRatio = telomere_length$ZAdjTSRatio.0.0)
test_exposures_telomere <- head(exposures_telomere, 1000)

## Removal of NAs: handled by imputation


# Analysis ----------------------------------------------------------------

## LASSO Models
set.seed(500)
train <- sample(1:nrow(test_exposures_telomere), 0.7 * nrow(test_exposures_telomere))
test <- seq(1, nrow(test_exposures_telomere))[-train]

set.seed(20)
lasso_model <- cv.glmnet(x = test_exposures, 
                         y = test_telomere$AdjTSRatio.0.0, alpha = 1)
plot(lasso_model)

## Regularization Path
lasso_path = glmnet(x = test_exposures, y = test_telomere$AdjTSRatio.0.0, alpha = 1, model = "gaussian")
plot(lasso_path)

# Work in Progress

## get lambda min and lambda.1se
lasso_model$lambda.min
lasso_model$ lambda.1se

bestlam_lasso <- model_lasso$lambda.1se

lassomodels.1se = cbind(lassomodels.1se, coef(model.lasso,
                                              s = "lambda.1se")[-1])


# summary of exposures used to predict telomere length------------------------------------------------

# numbers of non-zero beta coefs
table(coef(lasso_model, s = "lambda.1se")[-1] != 0)

# visualise non-zero beta coefficients
betas = coef(lasso_model, s = "lambda.1se")[-1]
names(betas) = rownames(coef(lasso_model, s = "lambda.1se"))[-1]
plot(betas[betas != 0], type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(beta))
axis(side = 1, at = 1:sum(betas != 0), labels = names(betas)[betas != 0], las = 2)
abline(h = 0, lty = 2)


## extract the beta coefs for calibrated model
lassomodels.1se = lassomodels = NULL
for (k in 1:ncol(proteins)) {
  print(k)
  model.lasso = cv.glmnet(x = test_exposures, y = test_telomere$AdjTSRatio.0.0[,
                                               k], alpha = 1, family = "gaussian")
  lassomodels.1se = cbind(lassomodels.1se, coef(model.lasso,
                                                s = "lambda.1se")[-1])
  lassomodels = cbind(lassomodels, coef(model.lasso,
                                        s = "lambda.min")[-1])
}
t1 = Sys.time()
print(t1 - t0)
colnames(lassomodels) = colnames(lassomodels.1se) = colnames(test_exposures)
rownames(lassomodels) = rownames(lassomodels.1se) = colnames(test_telomere$AdjTSRatio.0.0)

# save the lasso models ---------------------------------------------------

ifelse(dir.exists("/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure"),"",dir.create("/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure"))
saveRDS(lassomodels.1se, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure/lasso_1se.rds")
saveRDS(lassomodels, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure/lasso_min.rds")




