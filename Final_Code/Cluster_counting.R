###############################################################################
# for each year:
# gather all pixels with changepoints in that year and probability > thd from all 5 rasters
# compute the clump size of the pixel clusters that have been gathered
# remove any pixels belonging to clusters with size < size_thd (# of pixels)
# count the remaining pixels and assign to the pixel count table
###############################################################################



###############################################################################
############################# Import Libraries ################################
###############################################################################

library(magrittr)
library(rgdal)
library(raster)
library(bfastSpatial)



###############################################################################
############################# Define Parameter ################################
###############################################################################

setwd("/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.4/Final_Output/Output_Table")
image_dir = "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.4/Final_Output/BEAST_4"
image_list = list.files(image_dir, "*.tif", full.names = TRUE)

years <- c(1984:2022)
probability_threshold = 0.3
size_threshold = 4
pixel_resolution = 30 * 30 # Landsat Resolution

output_file_name = "LargeDam_Count_BEAST_4_03_4_Pos.csv"


###############################################################################
############################## Define Function ################################
###############################################################################

countChangePixels <- function(X, year, p_thd = probability_threshold, size_thd = size_threshold) {
  
  # X: change prob/magn + date stack (10 layers)
  probs <- X[[c(1:5)]]
  dates <- floor(X[[c(6:10)]])
  dates[probs < p_thd] <- NA
  dates[dates != year] <- NA
  dates[!is.na(dates)] <- 1
  nchanges <- calc(dates, sum, na.rm = TRUE)
  nchanges[nchanges > 1] <- 1
  
  clump_size <- clumpSize(nchanges)
  nchanges[clump_size <= size_thd] <- 0
  
  count = cellStats(nchanges, 'sum')
  change_area = count * pixel_resolution
  return(change_area)
}

###############################################################################
############################## For-Loop Starts ################################
###############################################################################

### Create data frame for appending rows:
list = list("File_Name", 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992,
            1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003,
            2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014,
            2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)

df = as.data.frame(do.call(cbind, list))

for (i in c(1:length(image_list))){
  input_raster = image_list[i]
  image = stack(input_raster)
  
  list = str_split(input_raster, "/")
  ID = list[[1]][length(list[[1]])]
  # cat("Currently working on:", ID, "\n")
  new_row = c(ID)
  for (year in years){
    cluster = countChangePixels(image,
                                year,
                                p_thd = probability_threshold,
                                size_thd = size_threshold)
    new_row = c(new_row, cluster)
  }
  new_row = as.data.frame(t(new_row))
  df = rbind(df, new_row)
  
  ### Write Progress
  progress = i / length(image_list) * 100
  cat("Current Progress: ", progress, "%", "\n")
  
}

write.csv(df, output_file_name)

print("All Done!")