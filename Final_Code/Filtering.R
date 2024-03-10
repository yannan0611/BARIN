###############################################################################
########################### Import Libraries ##################################
###############################################################################

library(rgdal)
library(ggplot2)
library(Rbeast)
library(bfast)
library(zoo)
library(magrittr)
library(bfastPlot)
library(raster)
library(doParallel)
library(foreach)
library(dplyr)
library(RColorBrewer)
library(devtools)
library(bfastSpatial)


###############################################################################
########################### Inputs Directory ##################################
###############################################################################

### Model Output Folder
model_output = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.4/LargeDam_Data/Output/BFAST_4"

### Lake Mask Folder
lake_mask_folder = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.4/TWI_Filter/LakeMask_warp"

### TWI Mask Folder
twi_mask_folder = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.4/TWI_Filter/TWImask_warp"

### Final Output Folder Path
output_folder_path = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.4/Final_Output/BFAST_4"
output_prefix = "BFAST_Model_4.tif"

### TWI Parameter
twi_parameter = 10

### Image List
model_list = list.files(model_output, "*.tif", full.names = TRUE)
lake_list = list.files(lake_mask_folder, "*.tif", full.names = TRUE)
twi_list = list.files(twi_mask_folder, "*.tif", full.names = TRUE)

model_name_list = list.files(model_output, "*.tif")


###############################################################################
############################### For Loop ######################################
###############################################################################

for (i in c(1 : length(model_list))){
  
  model_image = stack(model_list[i])
  lake_filter = raster(lake_list[i])
  twi_filter = raster(twi_list[i])
  
  ### Convert TWI Raster Into Binary Mask
  twi_filter[twi_filter < twi_parameter] <- 0
  twi_filter[twi_filter >= twi_parameter] <- 1
  
  final_output = model_image * lake_filter * twi_filter
  
  ### Write Outputs
  output_name = model_name_list[i]
  output_dir = paste(output_folder_path, output_name, sep = "/")
  writeRaster(final_output, output_dir, overwrite = TRUE)
  
  ### Write Progress
  progress = i / length(model_list) * 100
  cat("Current Progress: ", progress, "%", "\n")
  
}

print("All Done!")
