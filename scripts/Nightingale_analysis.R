suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
install.packages('glmnet')
suppressPackageStartupMessages(library(glmnet))

df <- data.frame(Nightingale_metabolomics_recoded)
genomics <- data.frame(Genomics_data_recoded)
df$TSR <- genomics$AdjTSRatio.0.0
pval <- rep(0,ncol(df_new)) 

##filtering out all those who did not take part in the metabolomics study (checked and for the first 10 or so columns it was the exact same number of those with/without NA so omitted those rows specifically from the first column for ease)
df_fil <- data_f %>%
  filter_at(vars(TotChol.0.0), all_vars(!is.na(.)))  

##filtering out any columns relating to technical covariates
df_new <- df_fil %>% select(-contains("QCF"))
df_new <- df_new %>% select(-contains("Spectr"))
df_new <- df_new %>% select(-contains("MQF"))
df_new <- df_new %>% select(-contains("Ship"))
df_new <- df_new %>% select(-contains("HLactate"))
df_new <- df_new %>% select(-contains("HPyruvate"))
df_new <- df_new %>% select(-contains("LGlucose"))
df_new <- df_new %>% select(-contains("LProtein"))

##running univariate analysis as seen in the comp epi wk 3 practical
for (j in 1:ncol(df_new)){
  model <- lm(TSR~.,data=df_new[,c(j,ncol(df_new))])
  coefs <- summary(model)$coefficients
  pval[j] <- coefs[2,'Pr(>|t|)']
}

plot(-log10(pval[1:ncol(df_new)]),pch=16,ylab="-log(p-values)",xlab="TS Ratio") # Manhattan plot
abline(h=-log10(0.05/ncol(df_new)), col='red',lty=2) # p-val = 0.05/p 
colnames(df_new[pval<=0.05/ncol(df_new)])
##125 significant markers

##plotting a heatmap
install.packages('pheatmap')
suppressPackageStartupMessages(library(pheatmap))

table(is.na(df_new$Glucose.0.0))

#getting rid of all other na values because it cant deal with them
df_om <- na.omit(df_new)

pheatmap(cor(df_om),  breaks = seq(-1, 1, length.out = 100),
         show_rownames = FALSE, show_colnames = FALSE)
##can see some clusters of positive correlation mostly

##running PCA for descriptive analysis

pcaX = prcomp(df_om)
ev = with(pcaX, sdev**2/sum(sdev**2))
plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
     ylab = "Proportion of EV", ylim = c(0, 1.2), cex = 0.3)
points(cumsum(ev), pch = 19, col = "skyblue", cex = 0.3)
legend("top", pch = 19, col = c("navy", "skyblue"),
       legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)

sum(!cumsum(ev)>0.9)
#only three principal components are needed to explain 90% of the variance

#distribution of Telomere length variable
summary(df_new$TSR)
##appears there are some high outliers as third quartile is 0.9077 but max is 5.36
(table(df_new$TSR>1))
#1 person with TSR bigger than 5, 1 person between 4-5, 3 people between 3-4, 7 people between 2-3, over 1000 between 1-2

ggplot(df_new, aes(x = TSR )) +
  geom_density()

fil <- df_new %>% filter(TSR >0.9077)
plot(density(fil$TSR))

fil2 <- df_new %>% filter(TSR > 2)
plot(density(fil2$TSR))

