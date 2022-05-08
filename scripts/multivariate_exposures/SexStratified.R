# Libraries ---------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(lme4))
suppressMessages(library(glmnet))

# dataset -----------------------------------------------------------------

df <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/dataForlasso.rds')
female <- subset(df, df$Sex_Female == "1")
male <- subset(df, df$Sex_Female == "0")
Xf = as.matrix(female[,1:169])
Yf = as.matrix(female[,170])
Xm = as.matrix(male[,1:169])
Ym = as.matrix(male[,170])

match("AgeAssess",names(df)) # 40
match("Sex_Female",names(df)) # 63

# female Lasso------------------------------------------------------------------
## LASSO Model: always keeping Sex_Female and AgeAssess
set.seed(10)
t0 = Sys.time()
lasso_model <- cv.glmnet(x = Xf, y = Yf,  
                         penalty.factor = c(rep(1, 39), 0, rep(1, 22), 0, rep(1, 106)))
t1 = Sys.time()
print(t1-t0)

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/lasso_model.pdf')
plot(lasso_model)
dev.off()

## Regularization Path
lasso_path = glmnet(x = Xf, y = Yf, penalty.factor = c(rep(1, 39), 0, rep(1, 22), 0, rep(1, 106)))
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/lasso_path.pdf')
plot(lasso_path)
dev.off()

## get lambda min and lambda.1se
lasso_model$lambda.min
lasso_model$ lambda.1se

### lambda.1se
beta_lasso_lambda.1se = coef(lasso_model, s = "lambda.1se")[2:(ncol(Xf) + 1), ]
saveRDS(beta_lasso_lambda.1se, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/beta_lasso.1se.rds")
selected_lasso_lambda.1se = names(beta_lasso_lambda.1se)[which(beta_lasso_lambda.1se != 0)]
saveRDS(selected_lasso_lambda.1se, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/selected_beta_lasso.1se.rds")
print(paste0(length(selected_lasso_lambda.1se), " exposures are selected when using lambda.1se"))

## visualise selected variables with their non-zero beta coefficients

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/lambda.1se_lasso_selected_exposures.pdf')
plot(beta_lasso_lambda.1se[beta_lasso_lambda.1se != 0], type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(beta_lasso_lambda.1se))
axis(side = 1, at = 1:sum(beta_lasso_lambda.1se != 0), labels = names(beta_lasso_lambda.1se)[beta_lasso_lambda.1se != 0], las = 2, cex.axis = 0.5)
abline(h = 0, lty = 2)
dev.off()


### lambda.min
beta_lasso_lambda.min = coef(lasso_model, s = "lambda.min")[2:(ncol(Xf) + 1), ]
saveRDS(beta_lasso_lambda.min, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/beta_lasso.min.rds")
selected_lasso_lambda.min = names(beta_lasso_lambda.min)[which(beta_lasso_lambda.min != 0)]
saveRDS(selected_lasso_lambda.min, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/selected_beta_lasso.min.rds")
print(paste0(length(selected_lasso_lambda.min), " exposures are selected when using lambda.min"))

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/lambda.min_lasso_selected_exposures.pdf')
plot(beta_lasso_lambda.min[beta_lasso_lambda.min != 0], type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(beta_lasso_lambda.min))
axis(side = 1, at = 1:sum(beta_lasso_lambda.min != 0), labels = names(beta_lasso_lambda.min)[beta_lasso_lambda.min != 0], las = 2, cex.axis = 0.35)
abline(h = 0, lty = 2)
dev.off()

# save the lasso models 
saveRDS(lasso_model, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/lasso_model.rds")


# female Stability --------------------------------------------------------
LoadPackages=function(packages){
  for (i in 1:length(packages)){
    suppressPackageStartupMessages(library(packages[i], character.only=TRUE))
  }
}

LoadPackages(c("focus","igraph","glmnet", "pheatmap"))

## running analysis
t0 = Sys.time()
out = VariableSelection(xdata = Xf, ydata = Yf, 
                        verbose = F,
                        family = "gaussian")
t1 = Sys.time()
print(t1-t0)
saveRDS(out, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Stability/out.RDS')

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Stability/CalibrationPlot.pdf')
CalibrationPlot(out, cex.lab = 0.5, cex.axis = 0.3)
dev.off()

## calibrated selection proportions
selprop = SelectionProportions(out)
print(selprop)
saveRDS(selprop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Stability/SelProp.RDS')
hat_params = Argmax(out) ## calibrated parameters
print(hat_params)
saveRDS(hat_params, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Stability/hat_params.RDS')

## Visualisation of selection proportions
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Stability/SelProp.pdf')
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
selected_lasso_lambda.min <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/selected_beta_lasso.min.rds")
selected_lasso_lambda.1se <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Lasso/selected_beta_lasso.1se.rds")

## plot

### selected_lasso_lambda.min
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Stability/SSvsLambda.min.pdf')
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

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Female/Stability/SSvsLambda.1se.pdf')
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


# male Lasso --------------------------------------------------------------------

## LASSO Model: always keeping Sex_Female and AgeAssess
set.seed(10)
t0 = Sys.time()
lasso_model <- cv.glmnet(x = Xm, y = Ym,  
                         penalty.factor = c(rep(1, 39), 0, rep(1, 22), 0, rep(1, 106)))
t1 = Sys.time()
print(t1-t0)

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/lasso_model.pdf')
plot(lasso_model)
dev.off()

## Regularization Path
lasso_path = glmnet(x = Xm, y = Ym, penalty.factor = c(rep(1, 39), 0, rep(1, 22), 0, rep(1, 106)))
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/lasso_path.pdf')
plot(lasso_path)
dev.off()

## get lambda min and lambda.1se
lasso_model$lambda.min
lasso_model$ lambda.1se

### lambda.1se
beta_lasso_lambda.1se = coef(lasso_model, s = "lambda.1se")[2:(ncol(Xm) + 1), ]
saveRDS(beta_lasso_lambda.1se, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/beta_lasso.1se.rds")
selected_lasso_lambda.1se = names(beta_lasso_lambda.1se)[which(beta_lasso_lambda.1se != 0)]
saveRDS(selected_lasso_lambda.1se, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/selected_beta_lasso.1se.rds")
print(paste0(length(selected_lasso_lambda.1se), " exposures are selected when using lambda.1se"))

## visualise selected variables with their non-zero beta coefficients

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/lambda.1se_lasso_selected_exposures.pdf')
plot(beta_lasso_lambda.1se[beta_lasso_lambda.1se != 0], type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(beta_lasso_lambda.1se))
axis(side = 1, at = 1:sum(beta_lasso_lambda.1se != 0), labels = names(beta_lasso_lambda.1se)[beta_lasso_lambda.1se != 0], las = 2, cex.axis = 0.5)
abline(h = 0, lty = 2)
dev.off()


### lambda.min
beta_lasso_lambda.min = coef(lasso_model, s = "lambda.min")[2:(ncol(Xm) + 1), ]
saveRDS(beta_lasso_lambda.min, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/beta_lasso.min.rds")
selected_lasso_lambda.min = names(beta_lasso_lambda.min)[which(beta_lasso_lambda.min != 0)]
saveRDS(selected_lasso_lambda.min, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/selected_beta_lasso.min.rds")
print(paste0(length(selected_lasso_lambda.min), " exposures are selected when using lambda.min"))

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/lambda.min_lasso_selected_exposures_new.pdf')
plot(beta_lasso.min[beta_lasso.min != 0], type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(beta_lasso.min))
axis(side = 1, at = 1:sum(beta_lasso.min != 0), labels = names(beta_lasso.min)[beta_lasso.min != 0], las = 2, cex.axis = 0.35)
abline(h = 0, lty = 2)
dev.off()

# save the lasso models 
saveRDS(lasso_model, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/lasso_model.rds")


# male Stability ----------------------------------------------------------

## running analysis
t0 = Sys.time()
out = VariableSelection(xdata = Xm, ydata = Ym, 
                        verbose = F, family = "gaussian")
t1 = Sys.time()
print(t1-t0)
saveRDS(out, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Stability/out.RDS')

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Stability/CalibrationPlot.pdf')
CalibrationPlot(out, cex.lab = 0.5, cex.axis = 0.2)
dev.off()

## calibrated selection proportions
selprop = SelectionProportions(out)
print(selprop)
saveRDS(selprop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Stability/SelProp.RDS')
hat_params = Argmax(out) ## calibrated parameters
print(hat_params)
saveRDS(hat_params, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Stability/hat_params.RDS')

## Visualisation of selection proportions
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Stability/SelProp.pdf')
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
selected_lasso_lambda.min <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/selected_beta_lasso.min.rds")
selected_lasso_lambda.1se <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Lasso/selected_beta_lasso.1se.rds")

## plot

### selected_lasso_lambda.min
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Stability/SSvsLambda.min.pdf')
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

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Sensitivity/Male/Stability/SSvsLambda.1se.pdf')
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



