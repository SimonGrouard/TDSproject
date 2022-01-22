rm(list=ls())
library(data.table)

# Loading an extract of the data (500 rows and 1,000 columns)
mydata <- fread("/rds/general/project/hda_21-22/live/TDS/General/Data/ukb47946.csv", data.table=FALSE, nrows=500, select=1:1000)
summary(mydata)
head(mydata)
