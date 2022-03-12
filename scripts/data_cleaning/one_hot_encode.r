suppressMessages(library(tidyverse))
suppressMessages(library(mltools))
suppressMessages(library(data.table))

expo_path <- "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/"
final_path <- paste(expo_path,"Exposures_covariates_recoded_combined_banded.rds",sep="")
expo<-readRDS(final_path)
print(dim(expo))

# change some factor variables to numeric/integer because they're not supposed to be factors

contains_letters <- function(x) {
  new <- sub("\\.","",x) # replace the first occurence of "." in a string by ""
  # so all doubles will be like integers
  return(grepl("\\D", new)) # true if contains any character that is not a number
  # including "." so that is why we removed them 
  # (otherwise doubles would have been TRUE)
}

convert_type <- function(column){ 
  # convert the colname variable to the right class and type, takes a dataframe as input
  # it filters rows containing letters if the class is supposed to be integer or numeric
  # it converts to factors if it is not integer or numeric (may be redundant)
  
  if (is.character(column)){ # then we convert to factor before moving on
    column <- as.factor(column)
  }
  if (is.factor(column)){ 
    column <- droplevels(column) # then we drop unused levels
  }
  
  col <- names(table(column)) # all the different categories of the variable X
  cat_letters <- col[sapply(col, FUN = contains_letters)] # all the categories containing letters
    
  if (sum(length(cat_letters)) < length(col)){ # ie it is not supposed to be a categorical variable
  
    without_letters <- column[!(column %in% cat_letters)][1] # the first element of the column (after removing all rows containing letters)
    if (!grepl("\\.", without_letters)) { # if the true class is integer
      # convert categories with letters to NA and the rest to integers
      if (is.factor(column)){ # it is currently factor instead of integer
        column <- suppressWarnings(as.integer(levels(column))[column]) 
      }else if(is.numeric(column)){ # it is currently numeric instead of integer
        column <- suppressWarnings(as.integer(column))
      }
    }
    else if(grepl("\\.", without_letters)){ # if the true class is numeric
      # then convert categories with letters to NA and the rest to numeric
      if (is.factor(column)){ # it is currently factor instead of numeric
        column <- droplevels(column)
        column <- suppressWarnings(as.numeric(levels(column))[column]) 
      }else if(is.integer(column)){ # it is currently integer instead of integer
        column <- suppressWarnings(as.numeric(column))
      }
    }
  }
  return(column)
}

expo[] <- lapply(expo, convert_type) # expo[] is used to keep the output as dataframe (because lapply returns a list normally, and sapply a vector)
                                     # and the apply function converts to matrix, which converts factors to character which is not wanted either

# take the subset of columns which are really factors
is_factor <- sapply(expo, is.factor)
expo_factors <- expo[,which(is_factor)]
print(dim(expo_factors))

# divide in two parts: equal to two levels, more than two levels
multilevel_factors <- sapply(expo_factors, function(x) {length(levels(x))>2})
expo_factors_superior2 <- expo_factors[,which(multilevel_factors)] # more than two levels
expo_factors_2 <- expo_factors[,which(!multilevel_factors)] # equal to two levels
print(dim(expo_factors_superior2))
print(dim(expo_factors_2))

#for (i in 1:dim(expo_factors_superior2)[2]){ # what are the factor variables with more than 10 levels
#  if (length(levels(expo_factors_superior2[,i])) > 10){
#    print(colnames(expo_factors_superior2)[i])
#    print(length(levels(expo_factors_superior2[,i])))
#    j = j+1
#  }
#}

# deal with factors which have two levels
change_colname <- function(x){
  if (levels(expo_factors_2[[x]])[1] != '0'){ # if it's not already one hot encoded
    x <- paste0(x,'_',levels(expo_factors_2[[x]])[1])
  }
  return(x)
}
colnames(expo_factors_2) <- sapply(colnames(expo_factors_2), change_colname)

change_column <- function(x){
  x <- as.factor(x)
  if (levels(x)[1] != '0'){ # if it's not already one hot encoded
    levels(x) <- c(1,0)
  }
  return(as.integer(levels(x))[x])
}
expo_factors_2 <- data.frame(apply(expo_factors_2, 2, change_column))

# deal with factors which have more than two levels
expo_factors_superior2 <- one_hot(as.data.table(expo_factors_superior2))

# put them back in expo
expo_final <- data.frame(expo, expo_factors_2, expo_factors_superior2)
expo_final <- expo_final %>% select(!colnames(expo_factors))
print(dim(expo_final))

saveRDS(expo_final, paste0(expo_path, "Exposures_covariates_recoded_combined_final_onehot.rds"))
