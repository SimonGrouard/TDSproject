library(devtools)
untar("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/focus_1.0.1.tar.gz")
install("focus", upgrade = "always")


LoadPackages=function(packages){
  for (i in 1:length(packages)){
    suppressPackageStartupMessages(library(packages[i], character.only=TRUE))
  }
}

LoadPackages(c("focus","igraph","glmnet", "pheatmap"))

df <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/dataForlasso.rds')

## test with first 1000 rows
# test_df <- as.matrix(head(df, 1000)) 
# X = test_df[,1:169]
# Y = test_df[,170]

X = df[,1:169]
Y = df[,170]

# Stability Selection -----------------------------------------------------

## running analysis
t0 = Sys.time()
out = VariableSelection(xdata = X, ydata = Y, 
                        verbose = F, penalty.factor =  c(rep(1, 39), 0, rep(1, 22), 0, rep(1, 106)), # age and sex are always selected
                        family = "gaussian")
out_all = VariableSelection(xdata = X, ydata = Y, 
                            verbose = F, 
                            family = "gaussian")
t1 = Sys.time()
print(t1-t0)
saveRDS(out, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/out.RDS')
saveRDS(out_all, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/out_all.RDS')

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/CalibrationPlot.pdf')
CalibrationPlot(out, cex.lab = 0.5, cex.axis = 0.2)
dev.off()

## calibrated selection proportions
selprop = SelectionProportions(out)
print(selprop)
saveRDS(selprop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SelProp.RDS')
hat_params = Argmax(out) ## calibrated parameters
print(hat_params)
saveRDS(hat_params, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/hat_params.RDS')

selprop_all = SelectionProportions(out_all)
print(selprop)
saveRDS(selprop_all, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SelProp_all.RDS')
hat_params_all = Argmax(out_all) ## calibrated parameters
print(hat_params_all)
saveRDS(hat_params_all, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/hat_params_all.RDS')


## Visualisation of selection proportions
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SelProp.pdf')
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

## Visualisation of selection proportions of all variables
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SelProp_all_new.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop_all, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(selprop_all >= hat_params_all[2], yes = "blue", no = "grey"),
     cex.lab = 1.5)
abline(h = hat_params_all[2], lty = 2, col = "darkred")
for (i in 1:length(selprop_all)){
  axis(side = 1, at = i, labels = names(selprop_all)[i], las = 2, 
       cex.axis = 0.25,
       col = ifelse(selprop_all >= hat_params_all[2], yes = "blue", no = "grey"), 
       col.axis = ifelse(selprop_all[i] >=hat_params_all[2], yes = "blue", no = "grey"))
}
dev.off()
# Comparison of exposure selection performed by different approaches ----------------------------

## Stability Selection vs. Lasso
## get the lasso results
selected_lasso_lambda.min <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Lasso/Lasso2/selected_beta_lasso.min.rds")
selected_lasso_lambda.1se <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Lasso/Lasso2/selected_beta_lasso.1se.rds")

## plot

### selected_lasso_lambda.min
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SSvsLambda.min.pdf')
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

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SSvsLambda.1se.pdf')
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



# plot for all variables --------------------------------------------------

### selected_lasso_lambda.min
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SSvsLambda.min.all.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop_all, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(selprop_all) %in% selected_lasso_lambda.min, yes = "blue", no = "grey"),
     cex.lab = 1.5,
     main = "Comparison between Stability Selection and Lambda.min")
abline(h = hat_params_all[2], lty = 2, col = "darkred")
for (i in 1:length(selprop_all)){
  axis(side = 1, at = i, labels = names(selprop_all)[i], las = 2, cex.axis = 0.25,
       col = ifelse(names(selprop_all) %in% selected_lasso_lambda.min, yes = "blue", no = "grey"),
       col.axis = ifelse(names(selprop_all)[i] %in% selected_lasso_lambda.min, yes = "blue", no = "grey"))
}
dev.off()

### selected_lasso_lambda.min
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SSvsLambda.1se.all.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop_all, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(selprop_all) %in% selected_lasso_lambda.1se, yes = "blue", no = "grey"),
     cex.lab = 1.5,
     main = "Comparison between Stability Selection and Lambda.1se")
abline(h = hat_params_all[2], lty = 2, col = "darkred")
for (i in 1:length(selprop_all)){
  axis(side = 1, at = i, labels = names(selprop_all)[i], las = 2, cex.axis = 0.25,
       col = ifelse(names(selprop_all) %in% selected_lasso_lambda.1se, yes = "blue", no = "grey"),
       col.axis = ifelse(names(selprop_all)[i] %in% selected_lasso_lambda.1se, yes = "blue", no = "grey"))
}
dev.off()

### SS vs univariate analyses
selected_uni <- subset(uni_res, uni_res$pval < 0.05/340) # Bonferroni multiple test correction
selected_uni_varnames <- rownames(selected_uni)

SPforPlot <- subset(SelProp_all, SelProp_all>= hat_params_all[2])
uni_SelProp <- subset(SelProp_all, names(SelProp_all)%in%rownames(selected_uni))


pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SSvsUni_new.all.pdf')
par(mar = c(10, 5, 5, 1))
plot(uni_SelProp, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(uni_SelProp)%in%names(SPforPlot), yes = "blue", no = "grey"),
     cex.lab = 1.5,
     main = "Comparison between Stability Selection and Univariate Analyses")
abline(h = hat_params_all[2], lty = 2, col = "darkred")
for (i in 1:length(uni_SelProp)){
  axis(side = 1, at = i, labels = names(uni_SelProp)[i], las = 2, cex.axis = 0.3,
       col = ifelse(uni_SelProp>= hat_params_all[2], yes = "blue", no = "grey"),
       col.axis = ifelse(uni_SelProp[i] >= hat_params_all[2], yes = "blue", no = "grey"))
}
dev.off()


# extract betas for stability selection
mean_beta <- apply(out_all$Beta, 2, mean)
saveRDS(mean_beta, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/mean_betas.rds')
## extract out those with >90% selection proportion:
sel_prop_90 <- SelProp_all[SelProp_all >= hat_params_all[2]]
mean_betas_sel_prop_90 <- mean_betas[names(mean_betas) %in% names(sel_prop_90)]

## Visualisation of selection proportions
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/StabilitySelection/SS_beta.pdf')
par(mar = c(10, 5, 5, 1))
plot(mean_betas_sel_prop_90, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Beta Coefficients of Stability Selection", xaxt = "n", 
     col = "blue",
     cex.lab = 1.5)
for (i in 1:length(mean_betas_sel_prop_90)){
  axis(side = 1, at = i, labels = names(mean_betas_sel_prop_90)[i], las = 2, 
       cex.axis = 0.6,
       col = "blue",
       col.axis = "blue")
}
dev.off()
