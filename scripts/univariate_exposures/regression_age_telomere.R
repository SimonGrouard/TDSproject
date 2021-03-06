rm(list=ls())
path<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path) # the data is not is this folder
getwd()

library(tidyverse)
library(patchwork)
library(kableExtra) # for pretty tables
library(gridExtra) # for pretty tables
library(pixiedust) # for pretty tables
library(jtools) # to use summ for glm

# Extract relevant data
extracted <- readRDS("../extraction_and_recording/outputs/ukb_extracted_AgeTSratio.rds")

dim(extracted)
extracted %>%
  select(AgeAssess.0.0, AdjTSRatio.0.0, ZadjTSRatio.0.0) %>% 
  lapply(is.na) %>%
  sapply(sum,simplify=TRUE)

AgeAdjTS <- extracted %>% 
  select(AgeAssess.0.0, AdjTSRatio.0.0) %>% 
  filter_at(vars(AgeAssess.0.0,AdjTSRatio.0.0),all_vars(!is.na(.)))
dim(AgeAdjTS)

# Run these lines to exclude outliers (> 1.5 times IQR)
Q <- quantile(AgeAdjTS$AdjTSRatio.0.0, probs=c(.25, .75), na.rm = T)
iqr <- IQR(AgeAdjTS$AdjTSRatio.0.0, na.rm = T)

AgeAdjTS <- AgeAdjTS %>% filter(AdjTSRatio.0.0 > (Q[1] - 1.5*iqr) & 
                                  AdjTSRatio.0.0 < (Q[2] + 1.5*iqr))  
dim(AgeAdjTS)
# Explanatory analysis before filtering observations
a.1<-extracted %>% #histogram for telomere length
  filter(AdjTSRatio.0.0 < 1.5,
         AdjTSRatio.0.0 > 0.4) %>% 
  ggplot(.) +
  geom_histogram(aes(AdjTSRatio.0.0), fill="Navy", bins = 30) +
  geom_vline(xintercept=quantile(extracted$AdjTSRatio.0.0, c(0.25, 0.75), na.rm=TRUE), colour="red")
summary(extracted$AdjTSRatio.0.0, digits=2)

b.1<-ggplot(data=extracted) + #histogram for age 
  geom_histogram(aes(AgeAssess.0.0), fill="orange",bins = 18) +
  geom_vline(xintercept=quantile(extracted$AgeAssess.0.0, c(0.25, 0.75), na.rm=TRUE), colour="red")
length(unique(extracted$AgeAssess.0.0))
summary(extracted$AgeAssess.0.0, digits=2)

c.1<-ggplot(data=extracted) + #boxplot for age
  geom_boxplot(aes(y=AgeAssess.0.0), fill="orange")

# Explanatory analysis after removing NAs
a.2<-AgeAdjTS %>% #histogram for telomere length
  filter(AdjTSRatio.0.0 < 1.5,
         AdjTSRatio.0.0 > 0.4) %>% 
  ggplot(.) +
  geom_histogram(aes(AdjTSRatio.0.0), fill="Navy", bins = 30) +
  geom_vline(xintercept=quantile(AgeAdjTS$AdjTSRatio.0.0, c(0.25, 0.75)), colour="red")
summary(AgeAdjTS$AdjTSRatio.0.0, digits=2)

b.2<-ggplot(data=AgeAdjTS) + #histogram for age 
  geom_histogram(aes(AgeAssess.0.0), fill="orange",bins = 18) +
  geom_vline(xintercept=quantile(AgeAdjTS$AgeAssess.0.0, c(0.25, 0.75)), colour="red")
length(unique(AgeAdjTS$AgeAssess.0.0))
summary(AgeAdjTS$AgeAssess.0.0, digits=2)

c.2<-ggplot(data=AgeAdjTS) + #boxplot for age
  geom_boxplot(aes(y=AgeAssess.0.0), fill="orange")

d<-ggplot(data=AgeAdjTS) + #boxplot for telomere length
  geom_boxplot(aes(y=AdjTSRatio.0.0), fill="navy")

# Compare before and after removing NAs
explanatory <- (b.1 | c.1) /
  (b.2 | c.2) # it doesn't affect the distribution of age to remove the NAs

ggsave("age_telomere_explanatory.png", plot = explanatory)

# Linear regression between telomere length and age
fit <- lm(AdjTSRatio.0.0 ~ AgeAssess.0.0, data = AgeAdjTS)

plot(fit, 1:2) # check for assumptions

reg <- ggplot(AgeAdjTS, aes(x = AgeAssess.0.0, y = AdjTSRatio.0.0)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") + 
  ylim(0, 2)

ggsave("age_telomere_regression.png", plot = reg)

# Table for estimates of the regression
summary_table_of_regression <- function(fit_reg, all_terms, row_var_interest, equation, title) { 
  a1<-data.frame(summ(fit_reg, confint = TRUE, ci.width = .95)$coeftable)
  CI<-paste(round(a1$X2.5.,digits=5),round(a1$X97.5.,digits=5),sep=" ; ")
  terms<-all_terms 
  cofint<-cbind(terms,data.frame(a1,CI))%>% relocate(CI,.after=terms)
  rownames(cofint)<-NULL
  cofint$terms<-as.character(cofint$terms)
  cofint$CI<-as.character(cofint$CI)
  # now create and plot the beautiful table
  kable(
    dust(cofint %>% select(col=c(1,2,3,7))) %>%
      sprinkle(col=3, round=5) %>%
      sprinkle(col=4, fn=quote(pvalString(value)))%>%
      sprinkle_colnames(col1="Terms",
                        col2="CI",
                        col3="Estimate",
                        col4="P-value") %>% 
      sprinkle_print_method("html"),
    align='c',booktabs = TRUE) %>%
    kable_styling(font_size = 15) %>%
    row_spec(0, color = 'black', background = 'lightgray') %>%
    row_spec(row_var_interest, bold = TRUE) %>% 
    add_header_above(equation, background = "lightblue") %>% 
    add_header_above(title, underline = TRUE, background = "lightblue")
}

all_terms <- c("Intercept", "Age")
equation <- c("Adjusted TS Ratio ~ Age"=4)
title <- c("Telomere length modeled by age of participant"=4)
summary_table_of_regression(fit, all_terms, c(2), equation, title)

# Linear regression between telemore length (z adjusted) and age
AgeZadjTS <- extracted %>% 
  select(AgeAssess.0.0, ZadjTSRatio.0.0) %>% 
  filter_at(vars(AgeAssess.0.0,ZadjTSRatio.0.0),all_vars(!is.na(.)))
dim(AgeZadjTS)

e.1<-AgeZadjTS %>% #histogram for telomere length
  filter(ZadjTSRatio.0.0 < 4,
         ZadjTSRatio.0.0 > -4) %>% 
  ggplot(.) +
  geom_histogram(aes(ZadjTSRatio.0.0), fill="Navy", bins = 30) +
  geom_vline(xintercept=quantile(extracted$ZadjTSRatio.0.0, c(0.25, 0.75), na.rm=TRUE), colour="red")

e.2<-ggplot(data=AgeZadjTS) + #boxplot for telomere length
  geom_boxplot(aes(y=ZadjTSRatio.0.0), fill="navy")

explanatory2 <- (a.1 | d) /
  (e.1 | e.2)
ggsave("age_telomere_explanatory2.png", plot = explanatory2)

fit2 <- lm(ZadjTSRatio.0.0 ~ AgeAssess.0.0, data = AgeZadjTS)
plot(fit2, 1:2) # check for assumptions
summary(fit2) # pretty much the same results
summary_table_of_regression(fit2, all_terms, c(2), equation, title)

reg2 <- ggplot(AgeZadjTS, aes(x = AgeAssess.0.0, y = ZadjTSRatio.0.0)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red")

ggsave("age_Ztelomere_regression.png", plot = reg2)




