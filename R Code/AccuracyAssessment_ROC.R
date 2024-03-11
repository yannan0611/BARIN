###############################################################################
############################# Import Libraries ################################
###############################################################################

library(magrittr)
library(rgdal)
library(raster)
library(terra)
library(dplyr)
library(bfastSpatial)


###############################################################################
############################# Define Parameter ################################
###############################################################################

setwd("/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Reference_Table")

##################################
### Cluster Counting Parameters:
##################################

image_dir = "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BEAST_1"
image_list = list.files(image_dir, "*.tif", full.names = TRUE)

dam_dir = "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/BeaverDamPoint/AlgonquinDamBuffer_60m.shp"
# dam_location = vect(dam_dir)

years <- c(1984:2022)

# X-Axis value for ROC
probability_thresholds = seq(0, 1, by = 0.1)

size_threshold = 2
pixel_resolution = 30 * 30 # Landsat Resolution


##################################
### Accuracy Assessment Parameters:
##################################

### Predicted file name
name_suffix = "_BEAST_Large.tif"

### Define the threshold for impoundment
impoundment_threshold = 1
impoundment_threshold =impoundment_threshold * 900

### ROC Output Name
ROC_outfile = "BEAST1_ROC.csv"

###############################################################################
############################## Define Function ################################
###############################################################################

countChangePixels <- function(X, year, p_thd = probability_threshold, size_thd = size_threshold, ID = ID){
  probs = X[[c(1:5)]]
  dates = floor(X[[c(6:10)]])
  dates[probs < p_thd] <- NA
  dates[dates != year] <- NA
  dates[!is.na(dates)] <- 1
  nchanges <- calc(dates, sum, na.rm = TRUE)
  nchanges[nchanges > 1] <- 1
  clump_id = clump(nchanges, directions = 8, gaps = FALSE)
  clump_size <- clumpSize(nchanges)
  clump_id[clump_size <= size_thd] <- NA
  # plot(clump_size)
  clump_id = rasterToPolygons(clump_id, f = function(x){x != 0}, dissolve = TRUE, na.rm = TRUE)
  # plot(clump_size)
  
  ## Beaver Dam Point Feature
  target_crs <- "+proj=utm +zone=17 +datum=WGS84 + units=m + no_defs"
  
  ### Find the corresponding polygon
  identifier = match$OGF_ID[match$File_Number == ID]
  dams = vect(dam_dir)
  dam_poly = dams[dams$OGF_ID == identifier]
  dam_poly = project(dam_poly, target_crs)
  # plot(dam_poly, add = TRUE)
  
  ### Iterate through all features
  total_impoundment = 0
  if (length(clump_id)!= 0){
    for (x in c(1:length(clump_id))){
      feature_x = clump_id[x,]
      vector_x = vect(feature_x)
      intersection = intersect(vector_x, dam_poly)
      impoundment_area = 0
      if (length(intersection$OGF_ID) > 0){
        impoundment_area = area(feature_x)
      }
      total_impoundment = total_impoundment + impoundment_area
      impoundment_area = 0
    }
  }
  # print(impoundment_area)
  return(total_impoundment)
}

###############################################################################
############################# For Loop Starts #################################
###############################################################################

### Output Rows
TP_list = c()
TN_list = c()
FP_list = c()
FN_list = c()



for (a in c(1:length(probability_thresholds))){
  
  #######################
  ####### Rename ########
  #######################
  # 1-10 represents 0 to 1 probability threshold
  name = as.character(a)
  ClusterCounting_Output = paste("BEAST1_ClusterCount_",name,".csv", sep = "")
  AccuracyTable_Output = paste("BEAST1_AccuracyTable_",name,".csv", sep = "")
  
  #######################
  
  probability_threshold = probability_thresholds[a]
  
  ### load the dataset
  match = read.csv("file_match.csv", header = TRUE, stringsAsFactors = FALSE)
  reference = read.csv("planet_validation.csv", header = TRUE, stringsAsFactors = FALSE)
  ### Convert "#VALUE!" and flooding events less than the threshold into NA
  reference[reference == "#VALUE!"] <- NA
  
  ### Convert tables into numeric values
  reference = as.data.frame(lapply(reference, as.numeric))
  match = as.data.frame(lapply(match, as.numeric))
  
  ### Subset the match dataset
  list = reference$OGF_ID
  subset_match = subset(match, OGF_ID %in% list)
  
  #################################################
  ######### Step 1: Cluster Counting ##############
  #################################################
  
  list = list("File_Name", 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992,
              1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003,
              2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014,
              2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)
  
  df = as.data.frame(do.call(cbind, list))
  
  for (i in c(1:length(image_list))){
    input_raster = image_list[i]
    image = brick(input_raster)
    
    list = str_split(input_raster, "/")
    ID = list[[1]][length(list[[1]])]
    # cat("Currently working on:", ID, "\n")
    new_row = c(ID)
    
    file_number = str_split(ID, "_")
    file_number = as.numeric(file_number[[1]][1])
    
    for (year in years){
      cluster = countChangePixels(image,
                                  year,
                                  p_thd = probability_threshold,
                                  size_thd = size_threshold, 
                                  ID = file_number)
      new_row = c(new_row, cluster)
    }
    new_row = as.data.frame(t(new_row))
    df = rbind(df, new_row)
    
    ### Write Progress
    progress = i / length(image_list) * 100
    cat("Threshold:", probability_threshold, ",", "Current Progress: ", progress, "%", "\n")
  }
  
  write.csv(df, ClusterCounting_Output)
  
  #################################################
  ############ Step 2: Accuracy table #############
  #################################################
  
  predicted = read.csv(ClusterCounting_Output, skip = 1, header = TRUE)
  
  
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
        accuracy = "NA"
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
            print("-1")
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
  

  write.csv(accuracy_table, AccuracyTable_Output)
  
  #################################################
  ### Step 3: Receiver Operating Characteristic ###
  #################################################
  
  accuracy_table = read.csv(AccuracyTable_Output)

  TP_freq = sum(accuracy_table == "TP", na.rm = TRUE)
  TN_freq = sum(accuracy_table == "TN", na.rm = TRUE)
  FP_freq = sum(accuracy_table == "FP", na.rm = TRUE)
  FN_freq = sum(accuracy_table == "FN", na.rm = TRUE)
  
  TP_list = c(TP_list, TP_freq)
  TN_list = c(TN_list, TN_freq)
  FP_list = c(FP_list, FP_freq)
  FN_list = c(FN_list, FN_freq)
  
}


###############################################################################
################################ Output ROC ###################################
###############################################################################

ROC_df = data.frame(
  probability = probability_thresholds,
  TP = TP_list,
  TN = TN_list,
  FP = FP_list,
  FN = FN_list
)

write.csv(ROC_df, ROC_outfile)

print("All done!")