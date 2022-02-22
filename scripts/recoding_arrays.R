suppressMessages(library(tidyverse))

expo_path <- "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined.rds"
expo<-readRDS(expo_path)

print(dim(expo))

all_arrays <- function(aim_expo, all_expo){
  # inputs: all_expo = all the colnames of expo / aim_expo = the CodingName of the variable
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

qualif_arrays <- all_arrays(colnames(expo), "Qualifications") # 6 variables
all_highest <- apply(expo[,qualif_arrays], 1, highest_level) 
expo <- data.frame(expo, Qualifications = all_highest) # add new column in expo dataframe
expo <- expo[,!(names(expo) %in% qualif_arrays)] # remove previous qualification columns

### One hot encoding 



### Pure binary 

pure_binary <- function(expo_truncated){
  level <- NA
  for (i in 1:length(expo_truncated)){
    if (expo_truncated[[i]] == "None of the above"){
      level <- 0
    } else if(expo_truncated[[i]] %in% c("Prefer not to answer", "Do not know")){
      level <- NA
    } else if(!is.na(expo_truncated[[i]])){
      level <- 1
    }
  }
  return(level)
}

to_binary <- c("GasSolidFuelCookHeat", "VitMinSupp", "DietarySupp", "FatherIllness", "MotherIllness", "SibIllness", "WPchemical", "WPsmok")
to_binary <- sapply(to_binary, all_arrays, all_expo = colnames(expo)) # 129 variables
all_binary <- 

### Remove the columns we didn't highlight

remove_col <- c("HomeHeatType", "CurrEmployment", "JobCommuteType", "VascHeartProbDiag", "ManicSymp", "PainLastMonth", "LeisureActivity", "PhysType4Wk", "WH_LC", "WPnoisy", "WPcold", "WPhot", "WPdusty", "ShiftWork", "DaySW", "WH_exact", "CancerFirst", "NonCancerFirst", "OpFirst")
remove_col <- unlist(sapply(remove_col, all_arrays, all_expo = colnames(expo))) # 419 variables!
expo <- expo[,!(names(expo) %in% remove_col)]

# save the result
print(dim(expo))
saveRDS(pvalues, expo_path)