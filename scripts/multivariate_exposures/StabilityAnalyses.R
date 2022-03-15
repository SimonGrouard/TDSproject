suppressMessages(library(dplyr))
suppressMessages(library(lme4))
suppressMessages(library(glmnet))


df <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/dataForlasso.rds')
x = as.matrix(df[,1:169])
y = as.matrix(df[,170])

# lasso models ---------------------------------------------------
lasso_model <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_model.rds")

#Stability analysis for lambda.min-----------------------
LassoSub = function(k=1, Xdata, Ydata, family="gaussian", penalty.factor=NULL) {
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
lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=x, Ydata=y)

apply(lasso.stab, 2, FUN=function(x){sum(x!=0)})

lasso.prop=apply(lasso.stab, 1, FUN=function(x){sum(x!=0)/length(x)})
names(lasso.prop)=colnames(x)

lasso.prop=sort(lasso.prop, decreasing = TRUE)

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Results_multivariate_exposures/Lasso/lambda.min_prop.pdf')
plot(lasso.prop[lasso.prop>0.2], type = 'h', col='navy', lwd=3, xaxt='n', 
     xlab='', ylab=expression(beta), ylim=c(0,1.2), las=1)
text(lasso.prop[lasso.prop>0.2]+0.07, labels = names(lasso.prop[lasso.prop>0.2]), 
     pos=3, srt=90, cex=0.7)
dev.off()

saveRDS(lasso.prop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_prop_lambda.min.rds')

#Stability analysis for lambda.1se-----------------------
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
lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=x, Ydata=y)

apply(lasso.stab, 2, FUN=function(x){sum(x!=0)})

lasso.prop=apply(lasso.stab, 1, FUN=function(x){sum(x!=0)/length(x)})
names(lasso.prop)=colnames(x)

lasso.prop=sort(lasso.prop, decreasing = TRUE)

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Results_multivariate_exposures/Lasso/lambda.1se_prop.pdf')
plot(lasso.prop[lasso.prop>0.2], type = 'h', col='navy', lwd=3, xaxt='n', 
     xlab='', ylab=expression(beta), ylim=c(0,1.2), las=1)
text(lasso.prop[lasso.prop>0.2]+0.07, labels = names(lasso.prop[lasso.prop>0.2]), 
     pos=3, srt=90, cex=0.7)
dev.off()

saveRDS(lasso.prop, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/Lasso/lasso_prop_lambda.1se.rds')

