
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

df <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/dataForlasso.rds')

## test with first 1000 rows
test_df <- as.matrix(head(df, 1000))
X = test_df[,1:169]
Y = test_df[,170]


# Analysis ----------------------------------------------------------------

## LASSO Models
set.seed(10)
t0 = Sys.time()
lasso_model <- cv.glmnet(x = X, 
                         y = Y, alpha = 1)
t1 = Sys.time()
print(t1-t0)

plot(lasso_model)

## Regularization Path
lasso_path = glmnet(x = X, y = Y, alpha = 1)
plot(lasso_path)

## get lambda min and lambda.1se
lasso_model$lambda.min # 0.0086
lasso_model$ lambda.1se # 0.0240

## final model
best_lam <- lasso_model$lambda.min # if use 1se, all beta_coefs will be 0
best_model <- glmnet(x = X, y = Y, lambda = best_lam)
lambda1se_model <- glmnet(x = X, y = Y, lambda = lasso_model$ lambda.1se)

## selected exposure variables
beta_lasso = matrix(coef(best_model, s = "lambda.1se")[-1], ncol = 169)
colnames(beta_lasso) <- colnames(X)
selected_lasso = colnames(beta_lasso)[which(beta_lasso != 0)]
print(selected_lasso)

## visualise selected variables with their non-zero beta coefficients
plot(beta_lasso[beta_lasso != 0],type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(best_lasso))
axis(side = 1, at = 1:sum(best_lasso != 0), labels = selected_lasso, las = 2)
abline(h = 0, lty = 2)


# save the lasso models ---------------------------------------------------
saveRDS(lasso_model, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_model.rds")
saveRDS(lambda1se_model, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_1se.rds")
saveRDS(best_model, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_min.rds")


# Stability Selection -----------------------------------------------------

## running analysis
t0 = Sys.time()
out = VariableSelection(xdata = X, ydata = Y, 
                        verbose = F, penalty.factor = c(rep(1, ncol(X))),
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

## traditional stability selection
#Stability analysis-----------------------
LassoSub = function(k=1, Xdata, Ydata, family="gaussian", penalty.factor=NULL) {
  if (is.null(penalty.factor)){
    penalty.factor=rep(1,ncol(Xdata))
  }
  set.seed(k)
  s = sample(nrow(Xdata), size = 0.8 * nrow(Xdata))
  Xsub = Xdata[s, ]
  Ysub = Ydata[s]
  model.sub = cv.glmnet(x = Xsub, y = Ysub, alpha = 1, family=family, penalty.factor=penalty.factor)
  coef.sub = coef(model.sub, s = "lambda.1se")[-1]
  return(coef.sub)
}

niter=100
lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=X, Ydata=Y)

apply(lasso.stab, 2, FUN=function(x){sum(x!=0)})

lasso.prop=apply(lasso.stab, 1, FUN=function(x){sum(x!=0)/length(x)})
names(lasso.prop)=colnames(X)

lasso.prop=sort(lasso.prop, decreasing = TRUE)
plot(lasso.prop[lasso.prop > 0], type = "h", col = "navy",
     lwd = 3, xaxt = "n", xlab = "", ylab = expression(beta),
     ylim = c(0, 1.2), las = 1)
text(lasso.prop[lasso.prop > 0.2] + 0.07, labels = names(lasso.prop[lasso.prop >
                                                                      0.2]), pos = 3, srt = 90, cex = 0.7)
