uni <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_nightingale_biomarker_univ/univ_pval_beta_adj_biomarker_2.rds')
#uni <- data.frame(uni)
#uni$V1 <- as.character(uni$V1)
#uni$V1 <- substr(uni$V1,1,nchar(uni$V1[1:50])-4)
#saveRDS(uni, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_nightingale_biomarker_univ/univ_pval_beta_adj_biomarker_2.rds')

#uni <- data.frame(uni)
uni$pval <- as.numeric(as.character(uni$pval))
uni <- arrange(uni, desc(pval))
uni$beta <- as.numeric(as.character(uni$beta))
uni <- uni[-c(53),]
uni <- uni[-c(52),]
uni <- uni[-c(51),]
uni <- uni[-c(8),]
uni <- uni[-c(8),]
uni <- uni[-c(2),]
uni <- uni[-c(19),]

uni_select <- uni %>% filter(pval <= 0.05)
uni_select <- arrange(uni_select, desc(beta))
selected_uni = (uni_select$V1)


selected_beta_lasso.1se <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/Lasso/Lasso2/selected_beta_lasso.1se.rds')


## plot

### SS vs selected_uni
selprop <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Biomarker_lasso/selection_proportions_bio_imp_new_2.rds')
hat_params <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Biomarker_lasso/hat_params_bio_imp_new_2.rds')


pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Biomarker_lasso/Sel_vs_univ_2.pdf')
par(mar = c(10, 5, 5, 1))
plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(names(selprop) %in% selected_uni, yes = "red", no = "grey"),
     cex.lab = 1.5,
     main = "Comparison between Stability Selection and Univariate Analysis")
abline(h = hat_params[2], lty = 2, col = "darkred")
for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprop)[i], las = 2, cex.axis = 0.7,
       col = ifelse(names(selprop) %in% selected_uni, yes = "red", no = "grey"),
       col.axis = ifelse(names(selprop)[i] %in% selected_uni, yes = "red", no = "grey"))
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

# Lambda.min vs Univariate ---------------------------------------------
## visualise selected variables with their non-zero beta coefficients
selected_beta_lasso.min = names(beta_lasso.min)[which(beta_lasso.min != 0)]
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ComparisonWithUni/beta_Lambda.min_vs_Uni.pdf')
par(mar = c(10, 5, 5, 1))
plot(uni_select$beta_coef, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "univariate beta_coefs", xaxt = "n", 
     col = ifelse(rownames(uni_select) %in% selected_beta_lasso.min, yes = "blue", no = "grey"),
     cex.lab = 1.5, 
     main = "Comparison between exposures selected by Lasso_Lambda.min and Univariate Analysis", cex.main = 0.8,
     sub = "(Both univariate and multivariate analyses were adjusted for age and sex)", cex.sub = 0.5)
for (i in 1:length(uni_select$beta_coef)) {
  axis(side = 1, at = i, labels = rownames(uni_select)[i], las = 2, cex.axis = 0.2, col = ifelse(rownames(uni_select)[i] %in% selected_beta_lasso.1se, yes = "blue", no = "grey"),
       col.axis = ifelse(rownames(uni_select)[i] %in% selected_beta_lasso.min,
                         yes = "blue", no = "grey"))
}
dev.off()


uni_select$inLassomin <- ifelse(rownames(uni_select)%in% selected_beta_lasso.min, yes = "1", no = "0")
uni_select$inLasso1se <- ifelse(rownames(uni_select)%in% selected_beta_lasso.1se, yes = "1", no = "0")
print(paste0(dim(subset(uni_select,uni_select$inLassomin == "1"))[1], " exposures identified by univariate analysis were also selected by multivariate analysis (lambda.min)"))
print(paste0(dim(subset(uni_select,uni_select$inLasso1se == "1"))[1], " exposures identified by univariate analysis were also selected by multivariate analysis (lambda.1se)"))


# y(p-values) ~ x(beta_coef)

## lambda.min
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ComparisonWithUni/Lambda.min_vs_Uni.pdf')
uni_select$Lambda.min <- ifelse(uni_select$name %in% selected_beta_lasso.min, yes = "selected by lasso.min", no = "not selected by lasso.min")
uni_select %>%
  ggplot(aes(x= beta_coef, y=log_pval, col=Lambda.min)) +
  geom_point() + 
  geom_vline(xintercept=c(-0.005, 0.005), col="red", lty = 2) +
  geom_hline(yintercept=-log10(0.05/dim(uni_select)[1]), col="red", lty = 2) +
  ggtitle("Volcano plot of comparing variable selected by univariate and Lasso models (Boferroni corrected)")
dev.off()
print(paste0(table(uni_select$Lambda.min)[2], " exposures selected by univariate models are also selected by Lasso model (lambda.min)"))

## lambda.1se
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ComparisonWithUni/Lambda.1se_vs_Uni.pdf')
uni_select$Lambda.1se <- ifelse(uni_select$name %in% selected_beta_lasso.1se, yes = "selected by lasso.1se", no = "not selected by lasso.1se")
uni_select %>%
  ggplot(aes(x= beta_coef, y=log_pval, col=Lambda.1se)) +
  geom_point() + 
  geom_vline(xintercept=c(-0.005, 0.005), col="red", lty = 2) +
  geom_hline(yintercept=-log10(0.05/dim(uni_select)[1]), col="red", lty = 2) +
  ggtitle("Volcano plot of comparing variable selected by univariate and Lasso models (Boferroni corrected)")
print(paste0(0, " exposures selected by univariate models are also selected by Lasso model (lambda.1se)"))
dev.off()


# y(selprop) ~ x(beta_coef)
sp <- data.frame(selprop)
sp$V1 <- rownames(sp)
colnames(sp) <- c("Selection_Proportion","V1")

join <- right_join(sp, uni_select, by= "V1")
join[is.na(join)] <- 0

s1 <-subset(join, join$Selection_Proportion >= 0.9)
join$significant_bio_in_univariate_analyses <- ifelse(join$V1 %in% s1$V1, yes = "Selected", no = "Not selected")
join <- join[-c(40),]
join$beta <- as.numeric(as.character(join$beta))

join$V1 <- substr(join$V1,1,nchar(uni$V1[1:40])-4)
rownames(join)<-join$V1
#library(ggplot2)
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Biomarker_lasso/Bio_comp_uni_multi_2.pdf')

  ggplot(data=join,aes(x= beta, y=Selection_Proportion, col=significant_bio_in_univariate_analyses)) +
  scale_color_discrete(name='Stability Selection Status', type=getOption('Dark2'))+
  geom_point() + 
  geom_hline(yintercept=0.9, col="red", lty = 2) +
  ggtitle("Comparison of biomarkers selected by univariate analyses and stability selection") +
  scale_x_continuous(name = 'Beta', breaks = c(-0.05, -0.030, -0.010, 0, 0.01), limits= c(-0.06,0.02) ) +
  scale_y_continuous(name = 'Selection Proportion') +
  geom_text(label = rownames(join), check_overlap = TRUE, size = 3.5, nudge_y = 0.02) +
  theme_bw()
dev.off()
print(paste0(table(join$significant_exposures_in_univariate_analyses)[2], " exposures selected by univariate models are also selected by stability selection"))


#betas vs pvalues

rownames(uni) <- uni$V1
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Biomarker_lasso/Bio_comp_uni_volc.pdf')

  ggplot(uni, aes(x= beta, y=-log(pval))) +
  geom_point() + 
  geom_hline(yintercept=c(2.995732,4.60517), col="red", lty = 2) +
  ggtitle("P-values vs Beta Values") +
  scale_x_continuous(name = 'Beta', breaks = c(-0.05, -0.030, -0.010, 0, 0.01)) +
  scale_y_continuous(name = '-Log P-value', limits=c(0,450)) +
  geom_text(label = rownames(uni), check_overlap = TRUE, size = 2.5, nudge_y = 12) +
  theme_bw()
dev.off()

?geom_text
