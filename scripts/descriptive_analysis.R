rm(list=ls())

# Libraries
library(data.table)

# Data
mydata <- fread("/rds/general/project/hda_21-22/live/TDS/General/Data/ukb47946.csv", data.table=FALSE, nrows=500, select=1:1000)
data_extract <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/ukb_extracted.rds")

# Reading arguments: For sh script to load data in future
args <- commandArgs(trailingOnly = TRUE)
ukb_path <- as.character(args[1])