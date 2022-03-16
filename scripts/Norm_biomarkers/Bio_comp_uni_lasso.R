uni <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_nightingale_biomarker_univ/univ_pval_beta_adj_biomarker.rds')
uni <- data.frame(uni)
uni <- arrange(uni, desc(pval))

uni_select <- subset(uni, uni$pval <= 0.05)
uni_select <- arrange(uni_select, desc(beta_coef))
selected_uni = rownames(uni_select)

## plot

### SS vs selected_uni
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Biomarker_lasso/Stab_vs_uni.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(selprop) %in% selected_uni, yes = "blue", no = "grey"),
     cex.lab = 1.5,
     main = "Comparison between Stability Selection and Univariate Analysis")
abline(h = hat_params[2], lty = 2, col = "darkred")
for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprop)[i], las = 2, cex.axis = 0.25,
       col = ifelse(names(selprop) %in% selected_uni, yes = "blue", no = "grey"),
       col.axis = ifelse(names(selprop)[i] %in% selected_uni, yes = "blue", no = "grey"))
}
dev.off()

# Lambda.1se vs Univariate ---------------------------------------------
## visualise selected variables with their non-zero beta coefficients
selected_beta_lasso.1se = names(beta_lasso.1se)[which(beta_lasso.1se != 0)]
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ComparisonWithUni/beta_Lambda.1se_vs_Uni.pdf')
par(mar = c(10, 5, 5, 1))
plot(uni_select$beta_coef, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "univariate beta_coefs", xaxt = "n", 
     col = ifelse(rownames(uni_select) %in% selected_beta_lasso.1se, yes = "blue", no = "grey"),
     cex.lab = 1.5, 
     main = "Comparison between exposures selected by Lasso_Lambda.1se and Univariate Analysis", cex.main = 0.8,
     sub = "(Both univariate and multivariate analyses were adjusted for age and sex)", cex.sub = 0.5)
for (i in 1:length(uni_select$beta_coef)) {
  axis(side = 1, at = i, labels = rownames(uni_select)[i], las = 2, cex.axis = 0.2, col = ifelse(rownames(uni_select)[i] %in% selected_beta_lasso.1se, yes = "blue", no = "grey"),
       col.axis = ifelse(rownames(uni_select)[i] %in% selected_beta_lasso.1se,
                         yes = "blue", no = "grey"))
}
dev.off()

uni_select$inLasso1se <- ifelse(rownames(uni_select)%in% selected_beta_lasso.1se, yes = "1", no = "0")
print(paste0(dim(subset(uni_select,uni_select$inLasso1se == "1"))[1], " exposures identified by univariate analysis were also selected by multivariate analysis (lambda.1se)"))
