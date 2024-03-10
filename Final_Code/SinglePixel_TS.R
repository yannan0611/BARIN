
###############################################################################
############################# Import Libraries ################################
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
############################# Define Parameters ###############################
###############################################################################

nbreaks = 5
time = "1 year"

# Input dataset directory
text_dir <- "C:/Users/wang25\OneDrive - University of Guelph/Desktop/M1.5/ROC_Updated/MNDWI_Input/Bands"
image_dir <- "C:/Users/wang25\OneDrive - University of Guelph/Desktop/M1.5/ROC_Updated/MNDWI_Input/Images"

# List all input files
text_list <- list.files(text_dir, "*.txt", full.names = TRUE)
image_list <- list.files(image_dir, "*.tif", full.names = TRUE)


for (i in c(1:length(text_list))){
  
  input_text <- text_list[i]
  input_raster <- image_list[i]
  
  # Read band information
  f = file(input_text)
  band_names = readLines(f)
  close(f)
  
  split_names <- strsplit(band_names, "_")
  dates <- sapply(split_names, FUN=function(x) {
    return(rev(x)[2])
  })
  
  dates <- as.Date(dates, "%Y%m%d")
  
  wetness_TS_stack = brick(input_raster)
  plot(wetness_TS_stack)
  
}
### check a pixel given a lon/lat coordinate pair

## 51, 143
# lon <- -133.6515744
# lat <- 68.7143048
46.095310,-78.864225

lon <- -78.864225
lat <- 46.095310

targ_crs <- projection(wetness_TS_stack)
pt <- SpatialPoints(t(c(lon, lat)), proj4string = CRS("+init=epsg:4269"))%>% spTransform(CRS = targ_crs)

# cellFromXY returns cell number
cell <- cellFromXY(wetness_TS_stack, pt)
# cell = 2212
X <- wetness_TS_stack[cell][1,]
X[X == 0] <- NA
plot(dates, X)
out = beast.irreg(X, 
                  time = dates,
                  deltat = time,
                  period = 1,
                  mcmc.seed = 123,
                  print.options = FALSE,
                  print.progress = FALSE,
                  quiet = TRUE,
                  gui = FALSE,
                  tcp.minmax = c(0, nbreaks),
                  season = "none",
                  # deseasonalize = TRUE,
                  tseg.min = 1
)
plot(out)


##########################
# Load Output Raster #

example_dir = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.3/Arctic_Sample/Positive_1000m/Sample_Outputs/BEAST/Deltat_1/7_BEAST.tif"
example_image = brick(example_dir)

probs = example_image[[1:5]]
dates = example_image[[6:10]]

dates[probs < 0.4] <- NA
probs[probs < 0.4] <- NA

plot(probs[[1]])
plot(dates[[1]])

probs_output = "C:/Users/wang25/Downloads/Example/Arctic/probs.tif"
dates_output = "C:/Users/wang25/Downloads/Example/Arctic/dates.tif"

writeRaster(probs, probs_output, overwrite = TRUE)
writeRaster(dates, dates_output, overwrite = TRUE)












##################################



library(ggplot2)
library(reshape2)
library(gridExtra)

# given a BEAST results called "out":


pos_cp <- out$trend$pos_cp
pos_cp_prob <- out$trend$pos_cpPr
pos_cp[pos_cp_prob < 0.3] <- NA
pos_cp <- pos_cp[!is.na(pos_cp)]

plotdf <- data.frame(time = out$time, mndwi = out$data, predicted = out$trend$Y, probability = out$trend$pos_cpOccPr)
plotdf <- na.omit(plotdf)


p1 <- ggplot(data = plotdf, aes(x = time, y = mndwi)) +
  geom_point() +
  geom_line(aes(y = predicted), alpha = 0.1, lwd = 1.2) +
  geom_vline(xintercept = pos_cp, col = 'blue', lty = 2, lwd = 1) +
  lims(y = c(-2000, 1000)) +
  labs(x = "", y = "MNDWI") +
  theme_bw()

p2 <- ggplot(data = plotdf, aes(x = time, y = probability)) +
  #geom_line(col = 'blue') +
  geom_bar(stat = "identity", position = "dodge", fill = 'blue', width = 0.5) +
  #geom_vline(xintercept = pos_cp, col = 'blue', lty = 2, lwd = 1) +
  #geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5, col = 'red', lty = 2) +
  lims(y = c(0, 1)) +
  labs(y = "Positive Change Probability") +
  theme_bw()

  
grid.arrange(p1, p2)
