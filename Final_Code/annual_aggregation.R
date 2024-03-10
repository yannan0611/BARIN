library(raster)

BEAST_raster = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.3/Arctic_Sample/Positive_1000m/Sample_Outputs/BEAST/Deltat_1/7_BEAST.tif"
BFAST_raster = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/M1.3/BEAST_test/Algonquin_VisibleDam/Point2/BFAST_Output.tif"

BEAST = brick(BEAST_raster)
BFAST = brick(BFAST_raster)

plot(BEAST[[6]])
plot(BFAST[[6]])

year = 2007
p_thd = 0.3

annual_aggregation <- function(X, year){
  
  ### Filter out break change year
  probs <- X[[c(1:5)]]
  dates <- floor(X[[c(6:10)]])
  
  # dates[dates != year] <- NA
  dates[dates < year - 1 | dates > year + 1] <- NA
  dates[probs < p_thd] <- NA
  dates[!is.na(dates)] <- 1
  nchanges <- calc(dates, sum, na.rm = TRUE)
  nchanges[nchanges > 1] <- 1
  
  ### Filter out break change probability
  
  probs <- X[[c(1:5)]]
  dates <- floor(X[[c(6:10)]])
  
  probs[dates < year - 1 | dates > year + 1] <- NA
  # probs[dates != year] <- NA
  probs[probs < p_thd] <- NA
  max_prob = calc(probs, max, na.rm = TRUE)
  
}

output_date = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/Ben/BEAST_date_2018-2020.tif"
output_prob = "C:/Users/wang25/OneDrive - University of Guelph/Desktop/Ben/BEAST_prob_2006-2008.tif"


writeRaster(nchanges, output_date, overwrite = TRUE)
writeRaster(max_prob, output_prob, overwrite = TRUE)
