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
################################ User Inputs ##################################
###############################################################################

# Input dataset directory
image_dir <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.4/LargeDam_Data/TS_Stack"
text_dir <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.4/LargeDam_Data/text_file"

# List all input files
text_list <- list.files(text_dir, "*.txt", full.names = TRUE)
image_list <- list.files(image_dir, "*.tif", full.names = TRUE)

# Output folder path
bfast_output_folder_path_1 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BFAST_1"
bfast_output_folder_path_2 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BFAST_2"
bfast_output_folder_path_3 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BFAST_3"
bfast_output_folder_path_4 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BFAST_4"

beast_output_folder_path_1 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BEAST_1"
beast_output_folder_path_2 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BEAST_2"
beast_output_folder_path_3 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BEAST_3"
beast_output_folder_path_4 <- "/Users/yannanwang/Library/CloudStorage/OneDrive-UniversityofGuelph/Desktop/M1.5/ROC_Updated/Output_Folder/BEAST_4"

# Output file name prefix (CHANGE THIS FOR EACH CLASS!!!!!)
bfast_prefix <- "BFAST_Large.tif"
beast_prefix <- "BEAST_Large.tif"

# NA Threshold for BEAST and BFAST
pr_threshold <- 0 # Probability threshold used for BEAST
nathresh <- 900 # 1200 is used for Algonquin Site for BFAST
nbreaks <- 5 # Number of positive breaks extracted

###############################################################################
############################ Multi-Core Detection #############################
###############################################################################
UseCores <- detectCores() - 1
cl <- makeCluster(UseCores)
registerDoParallel(cl)

##############################################################################
########################## BEAST & BFAST Definition ##########################
##############################################################################


################################## BFAST_1 ###################################
### "BIC" + "Harmonic & Trend" + "h = 0.05"
MNDWI_Characterization_BFAST_1 <- function (x) {
  
  x[x == 0] <- NA
  
  sumNA <- sum(is.na(x))
  if(sumNA > nathresh) {
    return(rep(NA, nbreaks*2))
  }
  
  x <- x / 10000
  
  bts <- bfastts(x, dates, type = "irregular")
  bpp <- bfastpp(bts, order = 1)
  
  breaks <- breakpoints(response ~ harmon+trend, data = bpp, breaks = "BIC", h = 0.05)
  
  bpp$segment <- sprintf("segment%s", length(breaks$breakpoints) + 1)
  
  if(is.na(breaks$breakpoints[1])){
    return(rep(NA, nbreaks*2))
  } else {
    i <- length(breaks$breakpoints)
    for(b in rev(breaks$breakpoints)) {
      bpp$segment[c(1:b)] <- sprintf("segment%s", i)
      i <- i - 1
    }
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
  
  
  ##########################################
  getMagnitude <- function(breakdate, trendOnly = FALSE) {
    data_bpp <- new_bpp
    
    if(trendOnly) {
      data_bpp$harmon[] <- 0
      data_bpp$prediction <- predict(model, data_bpp)
    }
    
    beforeBreak_bpp <- subset(data_bpp, time < breakdate)
    afterBreak_bpp <- subset(data_bpp, time > breakdate)
    
    pre_predicted <- beforeBreak_bpp$prediction[nrow(beforeBreak_bpp)]
    post_predicted <- afterBreak_bpp$prediction[1]
    
    return(post_predicted - pre_predicted)
  }
  
  magnitude <- sapply(breakdates, getMagnitude, trendOnly = TRUE)
  
  
  df <- data.frame(break_magntidue = c(magnitude), breakdate = c(breakdates))
  filtered_df <- df[df$break_magntidue >= 0, ]
  sorted_df <- filtered_df[order(filtered_df$break_magntidue, decreasing = TRUE),]
  
  final_magnitude_list <- c(sorted_df$break_magntidue, rep(NA, nbreaks))
  final_date_list <- c(sorted_df$breakdate, rep(NA, nbreaks))
  
  final_magnitude <- final_magnitude_list[1:nbreaks]
  final_date <- final_date_list[1:nbreaks]
  
  res <- c(final_magnitude, final_date)
  return(res)
}

################################## BFAST_2 ###################################
### "5" + "Harmonic & Trend" + "h = 0.05"

MNDWI_Characterization_BFAST_2 <- function (x) {
  
  x[x == 0] <- NA
  
  sumNA <- sum(is.na(x))
  if(sumNA > nathresh) {
    return(rep(NA, nbreaks*2))
  }
  
  x <- x / 10000
  
  bts <- bfastts(x, dates, type = "irregular")
  bpp <- bfastpp(bts, order = 1)
  
  breaks <- breakpoints(response ~ harmon+trend, data = bpp, breaks = 5, h = 0.05)
  
  bpp$segment <- sprintf("segment%s", length(breaks$breakpoints) + 1)
  
  if(is.na(breaks$breakpoints[1])){
    return(rep(NA, nbreaks*2))
  } else {
    i <- length(breaks$breakpoints)
    for(b in rev(breaks$breakpoints)) {
      bpp$segment[c(1:b)] <- sprintf("segment%s", i)
      i <- i - 1
    }
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
  
  
  ##########################################
  getMagnitude <- function(breakdate, trendOnly = FALSE) {
    data_bpp <- new_bpp
    
    if(trendOnly) {
      data_bpp$harmon[] <- 0
      data_bpp$prediction <- predict(model, data_bpp)
    }
    
    beforeBreak_bpp <- subset(data_bpp, time < breakdate)
    afterBreak_bpp <- subset(data_bpp, time > breakdate)
    
    pre_predicted <- beforeBreak_bpp$prediction[nrow(beforeBreak_bpp)]
    post_predicted <- afterBreak_bpp$prediction[1]
    
    return(post_predicted - pre_predicted)
  }
  
  magnitude <- sapply(breakdates, getMagnitude, trendOnly = TRUE)
  
  
  df <- data.frame(break_magntidue = c(magnitude), breakdate = c(breakdates))
  filtered_df <- df[df$break_magntidue >= 0, ]
  sorted_df <- filtered_df[order(filtered_df$break_magntidue, decreasing = TRUE),]
  
  final_magnitude_list <- c(sorted_df$break_magntidue, rep(NA, nbreaks))
  final_date_list <- c(sorted_df$breakdate, rep(NA, nbreaks))
  
  final_magnitude <- final_magnitude_list[1:nbreaks]
  final_date <- final_date_list[1:nbreaks]
  
  res <- c(final_magnitude, final_date)
  return(res)
}

################################## BFAST_3 ###################################
### "BIC" + "Trend" + "h = 0.05"

MNDWI_Characterization_BFAST_3 <- function (x) {
  
  x[x == 0] <- NA
  
  sumNA <- sum(is.na(x))
  if(sumNA > nathresh) {
    return(rep(NA, nbreaks*2))
  }
  
  x <- x / 10000
  
  bts <- bfastts(x, dates, type = "irregular")
  bpp <- bfastpp(bts, order = 1)
  
  breaks <- breakpoints(response ~ trend, data = bpp, breaks = "BIC", h = 0.05)
  
  bpp$segment <- sprintf("segment%s", length(breaks$breakpoints) + 1)
  
  if(is.na(breaks$breakpoints[1])){
    return(rep(NA, nbreaks*2))
  } else {
    i <- length(breaks$breakpoints)
    for(b in rev(breaks$breakpoints)) {
      bpp$segment[c(1:b)] <- sprintf("segment%s", i)
      i <- i - 1
    }
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
  
  
  ##########################################
  getMagnitude <- function(breakdate, trendOnly = FALSE) {
    data_bpp <- new_bpp
    
    if(trendOnly) {
      data_bpp$harmon[] <- 0
      data_bpp$prediction <- predict(model, data_bpp)
    }
    
    beforeBreak_bpp <- subset(data_bpp, time < breakdate)
    afterBreak_bpp <- subset(data_bpp, time > breakdate)
    
    pre_predicted <- beforeBreak_bpp$prediction[nrow(beforeBreak_bpp)]
    post_predicted <- afterBreak_bpp$prediction[1]
    
    return(post_predicted - pre_predicted)
  }
  
  magnitude <- sapply(breakdates, getMagnitude, trendOnly = TRUE)
  
  
  df <- data.frame(break_magntidue = c(magnitude), breakdate = c(breakdates))
  filtered_df <- df[df$break_magntidue >= 0, ]
  sorted_df <- filtered_df[order(filtered_df$break_magntidue, decreasing = TRUE),]
  
  final_magnitude_list <- c(sorted_df$break_magntidue, rep(NA, nbreaks))
  final_date_list <- c(sorted_df$breakdate, rep(NA, nbreaks))
  
  final_magnitude <- final_magnitude_list[1:nbreaks]
  final_date <- final_date_list[1:nbreaks]
  
  res <- c(final_magnitude, final_date)
  return(res)
}

################################## BFAST_4 ###################################
### "BIC" + "Harmonic & Trend" + "h = 0.15"

MNDWI_Characterization_BFAST_4 <- function (x) {
  
  x[x == 0] <- NA
  
  sumNA <- sum(is.na(x))
  if(sumNA > nathresh) {
    return(rep(NA, nbreaks*2))
  }
  
  x <- x / 10000
  
  bts <- bfastts(x, dates, type = "irregular")
  bpp <- bfastpp(bts, order = 1)
  
  breaks <- breakpoints(response ~ harmon+trend, data = bpp, breaks = "BIC", h = 0.15)
  
  bpp$segment <- sprintf("segment%s", length(breaks$breakpoints) + 1)
  
  if(is.na(breaks$breakpoints[1])){
    return(rep(NA, nbreaks*2))
  } else {
    i <- length(breaks$breakpoints)
    for(b in rev(breaks$breakpoints)) {
      bpp$segment[c(1:b)] <- sprintf("segment%s", i)
      i <- i - 1
    }
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
  
  
  ##########################################
  getMagnitude <- function(breakdate, trendOnly = FALSE) {
    data_bpp <- new_bpp
    
    if(trendOnly) {
      data_bpp$harmon[] <- 0
      data_bpp$prediction <- predict(model, data_bpp)
    }
    
    beforeBreak_bpp <- subset(data_bpp, time < breakdate)
    afterBreak_bpp <- subset(data_bpp, time > breakdate)
    
    pre_predicted <- beforeBreak_bpp$prediction[nrow(beforeBreak_bpp)]
    post_predicted <- afterBreak_bpp$prediction[1]
    
    return(post_predicted - pre_predicted)
  }
  
  magnitude <- sapply(breakdates, getMagnitude, trendOnly = TRUE)
  
  
  df <- data.frame(break_magntidue = c(magnitude), breakdate = c(breakdates))
  filtered_df <- df[df$break_magntidue >= 0, ]
  sorted_df <- filtered_df[order(filtered_df$break_magntidue, decreasing = TRUE),]
  
  final_magnitude_list <- c(sorted_df$break_magntidue, rep(NA, nbreaks))
  final_date_list <- c(sorted_df$breakdate, rep(NA, nbreaks))
  
  final_magnitude <- final_magnitude_list[1:nbreaks]
  final_date <- final_date_list[1:nbreaks]
  
  res <- c(final_magnitude, final_date)
  return(res)
}




################################## BEAST_1 ###################################
### "5" + "Deseasonalize" + "1 month"

MNDWI_Characterization_BEAST_1 <- function (x) {
  
  x[x == 0] <- NA
  
  out = try(beast.irreg(x,
                        time = dates,
                        deltat = "1 month",
                        period = "12 months",
                        season = "harmonic",
                        mcmc.seed = 123,
                        print.options = FALSE,
                        print.progress = FALSE,
                        quiet = FALSE,
                        gui = FALSE,
                        tcp.minmax = c(5, 5),
                        tseg.min = 12)
  )
  if(class(out) == "try-error") {
    res <- rep(-999, nbreaks*2)
    return(res)
  }
  
  pos_trend_cpPr = out$trend$pos_cpPr
  pos_trend_cpDate = out$trend$pos_cp
  
  pos_trend_cpPr[is.nan(pos_trend_cpPr)] <- 0
  
  if (length(pos_trend_cpDate) > 0){
    for (i in c(1:nbreaks)){
      if (pos_trend_cpPr[i] < pr_threshold){
        pos_trend_cpPr[i] <- NA
        pos_trend_cpDate[i] <- NA
      }
    }
  } else{
    pos_trend_cpPr <- rep(0, nbreaks)
    pos_trend_cpDate <- rep(0, nbreaks)
  }
  
  res <- c(pos_trend_cpPr, pos_trend_cpDate)
  
  return(res)
}

################################## BEAST_2 ###################################
### "0 - 5" + "Season = None" + "1 month"

MNDWI_Characterization_BEAST_2 <- function (x) { 
  
  x[x == 0] <- NA
  
  out = try(beast.irreg(x,
                        time = dates,
                        deltat = "1 month",
                        season = "none",
                        mcmc.seed = 123,
                        print.options = FALSE,
                        print.progress = FALSE,
                        quiet = TRUE,
                        gui = FALSE,
                        tcp.minmax = c(0, nbreaks),
                        tseg.min = 12)
  )
  if(class(out) == "try-error") {
    res <- rep(-999, nbreaks*2)
    return(res)
  }
  
  pos_trend_cpPr = out$trend$pos_cpPr
  pos_trend_cpDate = out$trend$pos_cp
  
  pos_trend_cpPr[is.nan(pos_trend_cpPr)] <- 0
  
  if (length(pos_trend_cpDate) > 0){
    for (i in c(1:nbreaks)){
      if (pos_trend_cpPr[i] < pr_threshold){
        pos_trend_cpPr[i] <- NA
        pos_trend_cpDate[i] <- NA
      }
    }
  } else{
    pos_trend_cpPr <- rep(0, nbreaks)
    pos_trend_cpDate <- rep(0, nbreaks)
  }
  
  res <- c(pos_trend_cpPr, pos_trend_cpDate)
  
  return(res)
}

################################## BEAST_3 ###################################
### "0 - 5" + "Trend + "1 year"

MNDWI_Characterization_BEAST_3 <- function (x) {
  
  x[x == 0] <- NA
  
  out = try(beast.irreg(x,
                        time = dates,
                        deltat = "1 year",
                        season = "none",
                        mcmc.seed = 123,
                        print.options = FALSE,
                        print.progress = FALSE,
                        quiet = TRUE,
                        gui = FALSE,
                        tcp.minmax = c(0, nbreaks),
                        tseg.min = 1)
  )
  if(class(out) == "try-error") {
    res <- rep(-999, nbreaks*2)
    return(res)
  }
  
  pos_trend_cpPr = out$trend$pos_cpPr
  pos_trend_cpDate = out$trend$pos_cp
  
  pos_trend_cpPr[is.nan(pos_trend_cpPr)] <- 0
  
  if (length(pos_trend_cpDate) > 0){
    for (i in c(1:nbreaks)){
      if (pos_trend_cpPr[i] < pr_threshold){
        pos_trend_cpPr[i] <- NA
        pos_trend_cpDate[i] <- NA
      }
    }
  } else{
    pos_trend_cpPr <- rep(0, nbreaks)
    pos_trend_cpDate <- rep(0, nbreaks)
  }
  
  res <- c(pos_trend_cpPr, pos_trend_cpDate)
  
  return(res)
}

################################## BEAST_4 ###################################
### "0 - 5" + "Deseasonalize" + "1 month"

MNDWI_Characterization_BEAST_4 <- function (x) {
  
  x[x == 0] <- NA
  
  out = try(beast.irreg(x,
                        time = dates,
                        deltat = "1 month",
                        period = "12 months",
                        season = "harmonic",
                        mcmc.seed = 123,
                        print.options = FALSE,
                        print.progress = FALSE,
                        quiet = TRUE,
                        gui = FALSE,
                        tcp.minmax = c(0, nbreaks),
                        tseg.min = 12)
  )
  if(class(out) == "try-error") {
    res <- rep(-999, nbreaks*2)
    return(res)
  }
  
  pos_trend_cpPr = out$trend$pos_cpPr
  pos_trend_cpDate = out$trend$pos_cp
  
  pos_trend_cpPr[is.nan(pos_trend_cpPr)] <- 0
  
  if (length(pos_trend_cpDate) > 0){
    for (i in c(1:nbreaks)){
      if (pos_trend_cpPr[i] < pr_threshold){
        pos_trend_cpPr[i] <- NA
        pos_trend_cpDate[i] <- NA
      }
    }
  } else{
    pos_trend_cpPr <- rep(0, nbreaks)
    pos_trend_cpDate <- rep(0, nbreaks)
  }
  
  res <- c(pos_trend_cpPr, pos_trend_cpDate)
  
  return(res)
}


###############################################################################
########################## For-loop Analysis Starts ###########################
###############################################################################

for (i in c(13:length(text_list))){
  
  input_text <- text_list[i]
  input_raster <- image_list[i]
  
  beast_name <- paste(i, beast_prefix, sep = "_")
  bfast_name <- paste(i, bfast_prefix, sep = "_")
  
  output_bfast_1 <- paste(bfast_output_folder_path_1, bfast_name, sep = "/")
  output_bfast_2 <- paste(bfast_output_folder_path_2, bfast_name, sep = "/")
  output_bfast_3 <- paste(bfast_output_folder_path_3, bfast_name, sep = "/")
  output_bfast_4 <- paste(bfast_output_folder_path_4, bfast_name, sep = "/")
  
  output_beast_1 <- paste(beast_output_folder_path_1, beast_name, sep = "/")
  output_beast_2 <- paste(beast_output_folder_path_2, beast_name, sep = "/")
  output_beast_3 <- paste(beast_output_folder_path_3, beast_name, sep = "/")
  output_beast_4 <- paste(beast_output_folder_path_4, beast_name, sep = "/")
  
  
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
  
  #############################################################################
  ########################### Calculation Begins ##############################
  #############################################################################
  
  print("Currently working on BEAST Model #1:")
  print(i)

  # BEAST_TS_Brick_1 = calc(wetness_TS_stack, MNDWI_Characterization_BEAST_1)
  BEAST_TS_Brick_1 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BEAST_1, mc.cores = UseCores) # , filename = output_beast, overwrite = TRUE
  writeRaster(BEAST_TS_Brick_1, output_beast_1, overwrite = TRUE)
  
  print("BEAST#1: Saved")

  print("Currently working on BEAST Model #2:")
  print(i)

  # BEAST_TS_Brick_2 = calc(wetness_TS_stack, MNDWI_Characterization_BEAST_2)
  BEAST_TS_Brick_2 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BEAST_2, mc.cores = UseCores) # , filename = output_beast, overwrite = TRUE
  writeRaster(BEAST_TS_Brick_2, output_beast_2, overwrite = TRUE)
  
  print("BEAST#2: Saved")

  print("Currently working on BEAST Model #3:")
  print(i)

  # BEAST_TS_Brick_3 = calc(wetness_TS_stack, MNDWI_Characterization_BEAST_3)
  BEAST_TS_Brick_3 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BEAST_3, mc.cores = UseCores) # , filename = output_beast, overwrite = TRUE
  writeRaster(BEAST_TS_Brick_3, output_beast_3, overwrite = TRUE)
  
  print("BEAST#3: Saved")

  print("Currently working on BEAST Model #4:")
  print(i)

  # BEAST_TS_Brick_4 = calc(wetness_TS_stack, MNDWI_Characterization_BEAST_4)
  BEAST_TS_Brick_4 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BEAST_4, mc.cores = UseCores) # , filename = output_beast, overwrite = TRUE
  writeRaster(BEAST_TS_Brick_4, output_beast_4, overwrite = TRUE)
  
  print("BEAST#4: Saved")
  
  print("Currently working on BFAST Model #1:")
  print(i)
  
  # BFAST_TS_Brick_1 = calc(wetness_TS_stack, MNDWI_Characterization_BFAST_1)
  BFAST_TS_Brick_1 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BFAST_1, mc.cores = UseCores) # , filename = output_bfast, overwrite = TRUE
  writeRaster(BFAST_TS_Brick_1, output_bfast_1, overwrite = TRUE)
  
  print("BFAST#1: Saved")
  
  print("Currently working on BFAST Model #2:")
  print(i)
  
  # BFAST_TS_Brick_2 = calc(wetness_TS_stack, MNDWI_Characterization_BFAST_2)
  BFAST_TS_Brick_2 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BFAST_2, mc.cores = UseCores) # , filename = output_bfast, overwrite = TRUE
  writeRaster(BFAST_TS_Brick_2, output_bfast_2, overwrite = TRUE)
  
  print("BFAST#2: Saved")
  
  print("Currently working on BFAST Model #3:")
  print(i)
  
  # BFAST_TS_Brick_3 = calc(wetness_TS_stack, MNDWI_Characterization_BFAST_3)
  BFAST_TS_Brick_3 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BFAST_3, mc.cores = UseCores) # , filename = output_bfast, overwrite = TRUE
  writeRaster(BFAST_TS_Brick_3, output_bfast_3, overwrite = TRUE)
  
  print("BFAST#3: Saved")
  
  print("Currently working on BFAST Model #4:")
  print(i)
  
  # BFAST_TS_Brick_4 = calc(wetness_TS_stack, MNDWI_Characterization_BFAST_4)
  BFAST_TS_Brick_4 <- mc.calc(wetness_TS_stack, fun = MNDWI_Characterization_BFAST_4, mc.cores = UseCores) # , filename = output_bfast, overwrite = TRUE
  writeRaster(BFAST_TS_Brick_4, output_bfast_4, overwrite = TRUE)
  
  print("BFAST#4: Saved")
}

print("All Done!")
