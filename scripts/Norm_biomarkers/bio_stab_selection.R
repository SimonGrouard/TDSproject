
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
lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=x, Ydata=y)

apply(lasso.stab, 2, FUN=function(x){sum(x!=0)})

lasso.prop=apply(lasso.stab, 1, FUN=function(x){sum(x!=0)/length(x)})
names(lasso.prop)=colnames(x)

lasso.prop=sort(lasso.prop, decreasing = TRUE)

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/plot_prop_imp.pdf')
plot(lasso.prop[lasso.prop>0.2], type = 'h', col='navy', lwd=3, xaxt='n', 
     xlab='', ylab=expression(beta), ylim=c(0,1.2), las=1)
text(lasso.prop[lasso.prop>0.2]+0.07, labels = names(lasso.prop[lasso.prop>0.2]), 
     pos=3, srt=90, cex=0.7)
dev.off()

saveRDS(lasso.prop, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_prop_imp.rds')




#stability selection 2-------------------------

out=VariableSelection(xdata=x, ydata=y, verbose=FALSE, 
                      penalty.factor=c(rep(1,ncol(x)),0,0), 
                      family="gaussian")

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/stab_selection.pdf')
CalibrationPlot(out)
dev.off()

# Calibrated selection proportions 
selprop=SelectionProportions(out)
selprop <- data.frame(selprop)

saveRDS(selprop, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/selection_proportions.rds')

# Calibrated parameters
hat_params=Argmax(out)
print(hat_params)

# Visualisation of selection proportions
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/vis_sel_prop.pdf')
par(mar=c(10,5,1,1))
plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(selprop>=hat_params[2], yes="red", no="grey"), cex.lab=1.5)
abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2, 
       col=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"),
       col.axis=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"))
}
dev.off()
