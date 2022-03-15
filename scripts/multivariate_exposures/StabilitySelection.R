install.packages("devtools")
library(devtools)
untar("/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/focus_1.0.1.tar.gz")
install("focus", upgrade = "always")

suppressMessages(library(dplyr))
suppressMessages(library(lme4))
suppressMessages(library(glmnet))
suppressMessages(library(focus))
suppressMessages(library(pheatmap))
suppressMessages(library(igraph))

df <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/dataForlasso.rds')

## test with first 1000 rows
test_df <- as.matrix(head(df, 1000))
X = test_df[,1:169]
Y = test_df[,170]


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
