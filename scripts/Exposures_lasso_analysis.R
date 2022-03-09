
# Libraries ---------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(lme4))
suppressMessages(library(glmnet))
suppressMessages(library(focus))
suppressMessages(library(pheatmap))

# parameters --------------------------------------------------------------

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])

# Data --------------------------------------------------------------------
telomere_length <- readRDS(here::here("extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- readRDS(here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds")) # dataset name to be confirmed
test_telomere <- head(telomere_length, 1000)
test_exposures <- head(exposures, 1000)

## Adding telomere length variable to exposure dataframe
exposures_telomere <- exposures %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0,
         ZAdjTSRatio = telomere_length$ZAdjTSRatio.0.0)
test_exposures_telomere <- head(exposures_telomere, 1000)

## Removal of NAs: will be handled by imputation


# Analysis ----------------------------------------------------------------

## LASSO Models
set.seed(10)
t0 = Sys.time()
lasso_model <- cv.glmnet(x = test_exposures, 
                         y = test_telomere$AdjTSRatio.0.0, alpha = 1)
t1 = Sys.time()
print(t1-t0)

plot(lasso_model)

## Regularization Path
lasso_path = glmnet(x = test_exposures, y = test_telomere$AdjTSRatio.0.0, alpha = 1)
plot(lasso_path)

## get lambda min and lambda.1se
lasso_model$lambda.min
lasso_model$ lambda.1se

## final model
best_lam <- model_lasso$lambda.1se # in case of including too many variables
best_model <- glmnet(x = test_exposures, 
                     y = test_telomere$AdjTSRatio.0.0, lambda = best_lam)

## selected exposure variables
table(coef(best_model,"lambda.1se"))[2:ncol(test_exposures)+1,]
beta_lasso = coef(best_model, s = "lambda.1se")[2:ncol(test_exposures)+1,]
selected_lasso = names(beta_lasso)[which(beta_lasso != 0)]
print(paste0(length(selected_lasso)), " exposures are selected")
print(selected_lasso)

## visualise selected variables with their non-zero beta coefficients
plot(beta_lasso[best_lasso != 0],type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(best_lasso))
axis(side = 1, at = 1:sum(best_lasso != 0), labels = selected_lasso, las = 2)
abline(h = 0, lty = 2)

## colnames and rownames
colnames(lassomodels) = colnames(lassomodels.1se) = colnames(test_exposures)
rownames(lassomodels) = rownames(lassomodels.1se) = colnames(test_telomere$AdjTSRatio.0.0)

# save the lasso models ---------------------------------------------------

ifelse(dir.exists("/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure"),"",dir.create("/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure"))
saveRDS(lasso_models, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure/lasso_min.rds")
saveRDS(best_model, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_Lasso_Exposure/lasso_1se.rds")


# Stability Selection -----------------------------------------------------

## running analysis
t0 = Sys.time()
out = VariableSelection(xdata = test_exposures, ydata = test_telomere$AdjTSRatio.0.0, 
                        verbose = F, penalty.factor = c(rep(1, ncol(test_exposures))),
                        family = "gaussian")
t1 = Sys.time()
print(t1-t0)
CalibrationPlot(out)

## calibrated selection proportions
selprop = SelectionProportions(out)
print(selprop)
hat_params = Argmax(out) ## calibrated parameters
print(hat_params)

## Visualistion of selection proportions
par(mar = c(10, 5, 5, 1))
plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"),
     cex.lab = 1.5)
abline(h = hat_params[2], lty = 2, col = "darkred")
for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprops)[i], las = 2, 
       col = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"),
       col.axis = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"))
}
    

# Comparison of exposure selection performed by different approaches ----------------------------

## Stability Selection vs. Lasso
par(mar = c(10, 5, 5, 1))
plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(selprop) %in% selected_lasso, yes = "blue", no = "grey"),
     cex.lab = 1.5)
abline(h = hat_params[2], lty = 2, col = "darkred")
for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprops)[i], las = 2, 
       col = ifelse(names(selprop) %in% selected_lasso, yes = "blue", no = "grey"),
       col.axis = ifelse(names(selprop) %in% selected_lasso, yes = "blue", no = "grey"))
}

## Univariate vs. Multivariate
