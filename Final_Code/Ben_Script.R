r <- raster(nrows = 10, ncols = 10)
values(r) <- 1:ncell(r)
plot(r)
f <- matrix(c(1,1,1,1,0,1,1,1,1), 3,3)
Moran(r, f)

input_raster = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.3/BEAST_test/Algonquin_VisibleDam/Point2/TS_Characterization.tif"
stack = brick(input_raster)
layer1 = stack[[1]]
Moran(layer1,f)
layer2 = stack[[2]]
plot(layer2)
Moran(layer2, f)
plot(stack[[6]])
Moran(stack[[6]],f)





####
library(bfastSpatial)

highest_prob = stack[[1]]
highest_prob_date = stack[[6]]
thd <- 0.5
highest_prob_date[highest_prob < thd] <- NA
clump_size <- clumpSize(highest_prob_date)
plot(highest_prob_date)

highest_prob_date[clump_size <= 2] <- NA
plot(highest_prob_date)



###

# for each year:
  # gather all pixels with changepoints in that year and probability > thd from all 5 rasters
  # compute the clump size of the pixel clusters that have been gathered
  # remove any pixels belonging to clusters with size < size_thd (# of pixels)
  # count the remaining pixels and assign to the pixel count table


library(magrittr)

years <- c(1985:2022)

countChangePixels <- function(X, year, p_thd = 0.3, size_thd = 2) {
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
  plot(nchanges)
  
  return(cellStats(nchanges, 'sum'))
}

for(year in years){cat(year, ":", countChangePixels(stack, year, p_thd = 0.3, size_thd = 5), "\n")}

