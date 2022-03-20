LoadPackages=function(packages){
  for (i in 1:length(packages)){
    suppressPackageStartupMessages(library(packages[i], character.only=TRUE))
  }
}

LoadPackages(c("focus","igraph","glmnet", "pheatmap"))

df <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/Exposures_imputed_full.rds')

## test with first 1000 rows
# test_df <- as.matrix(head(df, 1000)) 
# X = test_df[,-52]
# Y = test_df[,52]

X = df[,-52]
Y = df[,52]
print(dim(X))
# Stability Selection -----------------------------------------------------

## running analysis
t0 = Sys.time()
out = VariableSelection(xdata = X, ydata = Y, 
                        verbose = F, penalty.factor =  c(rep(1, 49), 0, rep(1, 13), 0, rep(1, 106)), # age and sex are always selected
                        family = "gaussian")
t1 = Sys.time()
print(t1-t0)
saveRDS(out, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_imputed/out.RDS')

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_imputed/CalibrationPlot.pdf')
CalibrationPlot(out, cex.lab = 0.5, cex.axis = 0.2)
dev.off()

## calibrated selection proportions
selprop = SelectionProportions(out)
print(selprop)
SelVar = names(selprop)[which(selprop != 0)]
print(paste0(length(SelVar), " variables were selected by stability selection"))
saveRDS(selprop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_imputed/SelProp.RDS')
hat_params = Argmax(out) ## calibrated parameters
print(hat_params)
saveRDS(hat_params, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_imputed/hat_params.RDS')

## Visualisation of selection proportions
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_imputed/SelProp.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(selprop >= hat_params[2], yes = "blue", no = "grey"),
     cex.lab = 1.5)
abline(h = hat_params[2], lty = 2, col = "darkred")
for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprop)[i], las = 2, 
       cex.axis = 0.25,
       col = ifelse(selprop >= hat_params[2], yes = "blue", no = "grey"), 
       col.axis = ifelse(selprop[i] >=hat_params[2], yes = "blue", no = "grey"))
}
dev.off()

# Comparison of exposure selection performed by different approaches ----------------------------

## Stability Selection vs. Lasso
## get the lasso results
selected_lasso_lambda.min <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Lasso/lasso_imputed/selected_beta_lasso.min.rds")
selected_lasso_lambda.1se <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Lasso/lasso_imputed/selected_beta_lasso.1se.rds")

## plot

### selected_lasso_lambda.min
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_imputed/SSvsLambda.min.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(selprop) %in% selected_lasso_lambda.min, yes = "blue", no = "grey"),
     cex.lab = 1.5,
     main = "Comparison between Stability Selection and Lambda.min")
abline(h = hat_params[2], lty = 2, col = "darkred")
for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprop)[i], las = 2, cex.axis = 0.25,
       col = ifelse(names(selprop) %in% selected_lasso_lambda.min, yes = "blue", no = "grey"),
       col.axis = ifelse(names(selprop)[i] %in% selected_lasso_lambda.min, yes = "blue", no = "grey"))
}
dev.off()

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_imputed/SSvsLambda.1se.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(selprop) %in% selected_lasso_lambda.1se, yes = "blue", no = "grey"),
     cex.lab = 1.5,
     main = "Comparison between Stability Selection and Lambda.1se")
abline(h = hat_params[2], lty = 2, col = "darkred")
for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprop)[i], las = 2, cex.axis = 0.25,
       col = ifelse(names(selprop) %in% selected_lasso_lambda.1se, yes = "blue", no = "grey"),
       col.axis = ifelse(names(selprop)[i] %in% selected_lasso_lambda.1se, yes = "blue", no = "grey"))
}
dev.off()

