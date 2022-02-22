suppressMessages(library(tidyverse))

expo_path <- "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined.rds"
expo<-readRDS(expo_path)

all_arrays <- function(all_expo, aim_expo){
  # inputs: all_expo = all the colnames of expo / aim_expo = the CodingName of the outcome variable
  # output: all the arrays related to the aim_expo exposure
  
  return(all_expo[startsWith(all_expo, aim_expo)])
}

### Highest level

highest_level <- function(expo_truncated){
  level <- NA
  for (i in 1:length(expo_truncated)){
    # could do it quicker with max(level, expo_truncated[[i]]) but would have to recode low, intermediate, high is numbers (and NAs would fuck it up)
    if(!is.na(expo_truncated[[i]])){
      if (is.na(level)){ # then the ith value is better than level
        level <- expo_truncated[[i]]
      }else if (level == "Low" && !is.na(expo_truncated[[i]])){ # then the ith value is better or equal to Low
        level <- expo_truncated[[i]]
      }else if (level == "Intermediate" && expo_truncated[[i]] == "High"){ # then the ith value is better than Intermediate
        level <- expo_truncated[[i]]
      }
    }
  }  
  return(level)
}

qualif_arrays <- all_arrays(colnames(expo), "Qualifications")
all_highest <- apply(expo[,qualif_arrays], 1, highest_level)
expo <- data.frame(expo, Qualifications = all_highest) # add new column in expo dataframe
expo <- expo[,!(names(expo) %in% qualif_arrays)] # remove previous qualification columns

### One hot encoding 


### Pure binary 




### Remove the columns we didn't highlight

remove_col <- c("HomeHeatType", "CurrEmployment", "JobCommuteType", "VascHeartProbDiag", "ManicSymp", "PainLastMonth", "LeisureActivity", "PhysType4Wk", "WH_LC", "WPnoisy", "WPcold", "WPhot", "WPdusty", "ShiftWork", "DaySW", "WH_exact", "CancerFirst", "NonCancerFirst", "OpFirst")

# save the result
saveRDS(pvalues, expo_path)