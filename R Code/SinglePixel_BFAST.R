### MNDWI BFAST segmentation
library(bfast)
library(zoo)
library(magrittr)
library(bfastPlot)


# setwd("D:/Yannans_Stuff/Semester_3/TimeSeries_Excel/TVC_TS")


# load data
df <- read.csv("D:/Yannans_Stuff/Semester_1/Data/Visible Dam/Landsat5_7_8_MNDWI_TS_2.csv")
df$date <- as.Date(df$date, "%Y-%m-%dT%H:%M:%S")


# remove positive outliers
df$MNDWI[df$MNDWI > 2000] <- NA
df <- na.omit(df)

# rescale MNDWI
df$MNDWI <- df$MNDWI / 10000


# create a 'regularized' time series from the MNDWI and dates columns
bts <- bfastts(df$MNDWI, df$date, type = 'irregular')


# initialize a data frame for harmonic and trend model fitting purposes
# we will use a harmonic model of the 1st order --> one season per cycle (year)
bpp <- bfastpp(bts, order = 1)

# detect breakpoints?
# `h` is the moving window size used to compute model residuals, and ultimately select the "best" segmentation model
# a smaller window size seems to allow for smaller segments
# it is a fraction between 0 and 1
# You can also change `breaks` to be an integer, representing the maximum allowable number of breakpoints
breaks <- breakpoints(response ~ harmon+trend, data = bpp, breaks = "BIC", h = 0.05)

# add a segment ID to the data frame according to the breakpoints
#bpp$segment <- breakfactor(breaks) ### not working as I would expect
bpp$segment <- sprintf("segment%s", length(breaks$breakpoints) + 1)

i <- length(breaks$breakpoints)
for(b in rev(breaks$breakpoints)) {
  bpp$segment[c(1:b)] <- sprintf("segment%s", i)
  i <- i - 1
}



# get formatted break dates
if(length(unique(bpp$segment)) == 1) {
  breakdates <- NA
} else {
  breakdates <- bpp$time[breaks$breakpoints]
  
}


# fit a piecewise harmonic+trend model to the data frame
# unless there are no breaks, then just fit the whole model at once
if(any(is.na(breakdates))) {
  model <- lm(response ~ harmon+trend, data = bpp)
} else {
  model <- lm(response ~ segment/(harmon+trend), data = bpp)
}


# extract the model coefficients
coef <- coefficients(model)

# trend slope for each segment
trends <- coef[which(grepl('trend', names(coef)))]
names(trends) <- unique(bpp$segment)

# model predictions for each time step
bpp$prediction <- predict(model, newdata = bpp)

# residuals
bpp$residual <- bpp$response - bpp$prediction

# smooth prediction
tt <- bfastts(bpp$response, time2date(bpp$time), type = "irregular") %>% time()
new_bpp <- bfastpp(tt, order = 1)
new_bpp$segment <- NA
for(i in c(length(breakdates):1)) {
  idx <- which(new_bpp$time < breakdates[i])
  new_bpp[idx,'segment'] <- sprintf('segment%s', i)
}

idx <- which(new_bpp$time >= breakdates[length(breakdates)])
new_bpp[idx,'segment'] <- sprintf('segment%s', length(breakdates) + 1)
new_bpp$prediction <- predict(model, new_bpp)

magnitude <- c()

for (breakdate in breakdates){
  for (i in 1:length(bpp$time)){
    if (bpp$time[i] == breakdate){
      break_magnitude = bpp$response[i + 1] - bpp$response[i]
      magnitude = append(magnitude, break_magnitude)
    }
  }
}


### Plots ###
lo <- matrix(c(1:2), nr=2, nc=1)
layout(lo)
op <- par(mar = c(0, 5, 0, 5), oma = c(3, 0, 2, 0))


# 1. Data overlaid with smooth model predictions
plot(prediction ~ time, data = new_bpp, type = 'l', col = '#0099FF', xaxt = 'n', lwd = 1.2, ylab = "MNDWI")
points(response ~ time, data = bpp, cex = 0.5, col = "#000000")
abline(v = breakdates, lty = 2)

# 2. Trend component only (modelled harmonic component subtracted from prediction)
plot(prediction ~ time, data = new_bpp, cex = 0, xaxt = 'n', ylab = "MNDWI Trend")
for (i in 1:length(unique(new_bpp$segment))) {
  seg <- unique(new_bpp$segment)[i]
  trend <- subset(new_bpp, segment == seg)
  trend$harmon[] <- 0
  trend$prediction <- predict(model, newdata = trend)
  lines(trend$time[trend$segment == seg], trend$prediction[trend$segment == seg], lwd = 2.5)
}

abline(v = breakdates, lty = 2)
axis(1, at = seq(1985, 2022, by = 5))
