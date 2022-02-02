mydata <- fread("/rds/general/project/hda_21-22/live/TDS/General/Data/ukb47946.csv", data.table=FALSE, nrows=500, select=1:1000)

# Reading arguments: For sh script to load data in future
args <- commandArgs(trailingOnly = TRUE)
ukb_path <- as.character(args[1])