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
library(gridExtra)


time = "1 months"
nbreaks <- 5

df <- read.csv("D:/Yannans_Stuff/Semester_1/Data/Visible Dam/Landsat5_7_8_MNDWI_TS_2.csv")
df$date <- as.Date(df$date, "%Y-%m-%dT%H:%M:%S")

# Because Rbeast only compute based on Time-Series value, we here need to create
# a time-serie object using bfastts() function

df$MNDWI[df$MNDWI > 2000] <- NA
df <- na.omit(df)

out = beast.irreg(df$MNDWI, time = df$date,
                  deltat = time,
                  tcp.minmax = c(0,nbreaks),
                  mcmc.seed = 123,
                  deseasonalize = TRUE,
                  tseg.min = 12)

plot(out)

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
  lims(y = c(-3000, 1000)) +
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

