rm(list=ls())

# Libraries
library(data.table)
library(ggplot2)

# Data
mydata <- fread("/rds/general/project/hda_21-22/live/TDS/General/Data/ukb47946.csv", data.table=FALSE, nrows=500, select=1:1000)
# data_extract <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/ukb_extracted.rds")
data_exploratory <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/ukb_extracted.rds", 
                            data.table=FALSE, nrows=500)

# Reading arguments: For sh script to load data in future
args <- commandArgs(trailingOnly = TRUE)
ukb_path <- as.character(args[1])


# Analysis ----------------------------------------------------------------
# Distribution of Sex
tab_sex <- table(data_exploratory$Sex.0.0)
sex_barplot <- barplot(tab_sex, beside=TRUE, ylim=c(0, max(tab_sex) + 100000), main = "Distribution of Sex")
text(sex_barplot, tab_sex + 25000, tab_sex)

# Distribution of DOB
summary(data_exploratory$YOB.0.0)
hist(data_exploratory$YOB.0.0, main = "Distribution of Years of Birth")

# Distribution of Telomere Length
