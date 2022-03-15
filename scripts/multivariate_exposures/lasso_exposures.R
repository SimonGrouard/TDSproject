# Libraries ---------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(lme4))
suppressMessages(library(glmnet))

df <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/dataForlasso.rds')
X = as.matrix(df[,1:169])
Y = as.matrix(df[,170])



# Analysis ----------------------------------------------------------------
## before setting the penalty factor
match("AgeAssess",names(df)) # 40
match("Sex_Female",names(df)) # 63
## LASSO Model: always keeping Sex_Female and AgeAssess
set.seed(10)
t0 = Sys.time()
lasso_model <- cv.glmnet(x = X, y = Y,  
                         penalty.factor = c(rep(1, 39), 0, rep(1, 19), 0, rep(1, 109)))
t1 = Sys.time()
print(t1-t0)

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_model.pdf')
plot(lasso_model)
dev.off()

## Regularization Path
lasso_path = glmnet(x = X, y = Y, penalty.factor = c(rep(1, 39), 0, rep(1, 19), 0, rep(1, 109)))
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_path.pdf')
plot(lasso_path)
dev.off()

## get lambda min and lambda.1se
lasso_model$lambda.min
lasso_model$ lambda.1se

### lambda.1se
beta_lasso_lambda.1se = coef(lasso_model, s = "lambda.1se")[2:(ncol(X) + 1), ]
saveRDS(beta_lasso_lambda.1se, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/beta_lasso.1se.rds")
selected_lasso_lambda.1se = names(beta_lasso_lambda.1se)[which(beta_lasso_lambda.1se != 0)]
saveRDS(selected_lasso_lambda.1se, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/selected_beta_lasso.1se.rds")
print(paste0(length(selected_lasso_lambda.1se), " exposures are selected when using lambda.1se"))

## visualise selected variables with their non-zero beta coefficients

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lambda.1se_lasso_selected_exposures.pdf')
plot(beta_lasso_lambda.1se[beta_lasso_lambda.1se != 0], type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(beta_lasso_lambda.1se))
axis(side = 1, at = 1:sum(beta_lasso_lambda.1se != 0), labels = names(beta_lasso_lambda.1se)[beta_lasso_lambda.1se != 0], las = 2, cex.axis = 0.5)
abline(h = 0, lty = 2)
dev.off()

### lambda.min
beta_lasso_lambda.min = coef(lasso_model, s = "lambda.min")[2:(ncol(X) + 1), ]
saveRDS(beta_lasso_lambda.min, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/beta_lasso.min.rds")
selected_lasso_lambda.min = names(beta_lasso_lambda.min)[which(beta_lasso_lambda.min != 0)]
saveRDS(selected_lasso_lambda.min, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/selected_beta_lasso.min.rds")
print(paste0(length(selected_lasso_lambda.min), " exposures are selected when using lambda.min"))

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lambda.min_lasso_selected_exposures.pdf')
plot(beta_lasso_lambda.min[beta_lasso_lambda.min != 0], type = "h", col = "navy", lwd = 3,
     xaxt = "n", xlab = "", ylab = expression(beta_lasso_lambda.min))
axis(side = 1, at = 1:sum(beta_lasso_lambda.min != 0), labels = names(beta_lasso_lambda.min)[beta_lasso_lambda.min != 0], las = 2, cex.axis = 0.15)
abline(h = 0, lty = 2)
dev.off()

# save the lasso models ---------------------------------------------------
saveRDS(lasso_model, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_model.rds")

# Stability analysis -----------------------
## Lambda.1se
LassoSub = function(k=1, Xdata, Ydata, family="gaussian", 
                    penalty.factor=c(rep(1, 39), 0, rep(1, 19), 0, rep(1, 109))) {
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
saveRDS(lasso.prop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso.1se_Prop.rds')

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Results_multivariate_exposures/Lasso/lambda.1se_PropPlot.pdf')
plot(lasso.prop[lasso.prop>0.2], type = 'h', col='navy', lwd=3, xaxt='n', 
     xlab='', ylab=expression(beta), ylim=c(0,1.2), las=1)
text(lasso.prop[lasso.prop>0.2]+0.07, labels = names(lasso.prop[lasso.prop>0.2]), 
     pos=3, srt=90, cex=0.7)
dev.off()




## Lambda.min
LassoSub = function(k=1, Xdata, Ydata, family="gaussian", 
                    penalty.factor=c(rep(1, 39), 0, rep(1, 19), 0, rep(1, 109))) {
  if (is.null(penalty.factor)){
    penalty.factor=rep(1,ncol(Xdata))
  }
  set.seed(k)
  s = sample(nrow(Xdata), size = 0.8 * nrow(Xdata))
  Xsub = Xdata[s, ]
  Ysub = Ydata[s]
  model.sub = cv.glmnet(x = Xsub, y = Ysub, alpha = 1, family=family, penalty.factor=penalty.factor)
  coef.sub = coef(model.sub, s = "lambda.min")[-1]
  return(coef.sub)
}

niter=100
lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=X, Ydata=Y)
apply(lasso.stab, 2, FUN=function(x){sum(x!=0)})
lasso.prop=apply(lasso.stab, 1, FUN=function(x){sum(x!=0)/length(x)})
names(lasso.prop)=colnames(X)

lasso.prop=sort(lasso.prop, decreasing = TRUE)
saveRDS(lasso.prop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso.min_Prop.rds')

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Results_multivariate_exposures/Lasso/lambda.min_PropPlot.pdf')
plot(lasso.prop[lasso.prop>0.2], type = 'h', col='navy', lwd=3, xaxt='n', 
     xlab='', ylab=expression(beta), ylim=c(0,1.2), las=1)
text(lasso.prop[lasso.prop>0.2]+0.07, labels = names(lasso.prop[lasso.prop>0.2]), 
     pos=3, srt=90, cex=0.7)
dev.off()




