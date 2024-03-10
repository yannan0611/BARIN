###############################################################################
############################# Import Libraries ################################
###############################################################################
library(dplyr)


###############################################################################
########################### Define Input Parameters ###########################
###############################################################################

### Set up working directory
dir = "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.4/Final_Output/Output_Table"
setwd(dir)

### Predicted file name
name_suffix = "_BEAST_NewLarge.tif"

### Define the threshold for impoundment
impoundment_threshold = 1
impoundment_threshold =impoundment_threshold * 900

### Output accuracy table name
output_name = "AccuracyTable_BEAST4_03_4.csv"

###############################################################################
############################ Data Pre-Processing ##############################
###############################################################################

### load the dataset
match = read.csv("file_match.csv", header = TRUE, stringsAsFactors = FALSE)
reference = read.csv("planet_validation.csv", header = TRUE, stringsAsFactors = FALSE)
predicted = read.csv("LargeDam_Count_BEAST_4_03_4_Pos.csv", skip = 1, header = TRUE)
### Convert "#VALUE!" and flooding events less than the threshold into NA
reference[reference == "#VALUE!"] <- NA

### Convert tables into numeric values
reference = as.data.frame(lapply(reference, as.numeric))
match = as.data.frame(lapply(match, as.numeric))

### Subset the match dataset
list = reference$OGF_ID
subset_match = subset(match, OGF_ID %in% list)

###############################################################################
############################# For Loop Starts #################################
###############################################################################

output_header = list("OGF_ID", "Match_ID", 2009, 2010, 2011, 2012, 2013, 2014, 
                     2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)
accuracy_table = as.data.frame(do.call(cbind, output_header))

for (i in (1 : length(subset_match$OGF_ID))) {
  file_count = subset_match$File_Number[i]
  OGF_ID = subset_match$OGF_ID[i]
  
  ### Retrieve observed reference values
  reference_value = reference[reference$OGF_ID == OGF_ID, ]
  
  ### Retrieve predicted values from BEAST/BFAST
  file_name = paste(file_count,name_suffix, sep = "")
  predicted_value = predicted[predicted$File_Name == file_name, ]
  
  reference_value = reference_value[, (ncol(reference_value) - 13):ncol(reference_value)]
  predicted_value = predicted_value[, (ncol(predicted_value) - 13):ncol(predicted_value)]
  
  ### Iterate through the row
  new_row = c(OGF_ID, file_count)
  
  for (x in (1 : length(reference_value))){
    
    ### Assign initial accuracy label
    if (is.na(reference_value[,x]) | is.na(predicted_value[,x])){
      accuracy = NA
    } else if (reference_value[,x] >= impoundment_threshold & predicted_value[,x] > 0){
      accuracy = "TP"
    } else if (reference_value[,x] < impoundment_threshold & predicted_value[,x] == 0){
      accuracy = "TN"
    } else if (reference_value[,x] < impoundment_threshold & predicted_value[,x] > 0){
      accuracy = "FP"
    } else if (reference_value[,x] >= impoundment_threshold & predicted_value[,x] == 0){
      accuracy = "FN"
    } 
    
    ### 1-year search window
    if (x <= 13){
      if (!is.na(reference_value[,x+1])){
        if (reference_value[,x+1] >= impoundment_threshold & predicted_value[,x] > 0){
          accuracy = "TP"
          print("+1")
        }
      }
    }
    
    if (x >= 2){
      if (!is.na(reference_value[,x-1])){
        if (reference_value[,x-1] >= impoundment_threshold & predicted_value[,x] > 0){
          accuracy = "TP"
          print(-1)
        }
      }
    }
    
    
    new_row = c(new_row, accuracy)
  }
  
  accuracy_table = rbind(accuracy_table, new_row)
  
  ### Write Progress
  progress = i / length(subset_match$OGF_ID) * 100
  cat("Current Progress: ", progress, "%", "\n")
}
  
write.csv(accuracy_table, output_name)

print("All Done!")
