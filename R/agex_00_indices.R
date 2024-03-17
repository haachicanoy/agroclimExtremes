## ------------------------------------------ ##
## Agro-climatic indices
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(caTools))

## ------------------------------------------ ##
## Extreme precipitation
## ------------------------------------------ ##

# Percentile 95 precipitation within growing season
p95 <- function(x){ # x: daily precipitation vector
  if(!all(is.na(x))){
    p <- 1:length(x) # Daily positions within growing season
    p <- p[x > 0]    # Daily positions with non-zero precipitation
    x <- x[x > 0]    # Daily rainfall with non-zero precipitation
    if(length(x) > 0){
      p95 <- stats::quantile(x, probs = 0.95, na.rm = T) # 95th percentile of daily precipitation
      whn <- p[which(x > p95)[1]] # Position when the 95th percentile occur
    } else {
      p95 <- NA
      whn <- NA
    }
  } else {
    p95 <- NA
    whn <- NA
  }
  return(c(p95,whn))
}

# Maximum 5-day running average precipitation
p5d <- function(x){ # x: daily precipitation vector
  if(!all(is.na(x))){
    p <- 1:length(x)
    p5d <- caTools::runmean(x, k = 5, endrule = 'NA')
    whn <- which.max(p5d) # Position when the P5D maximum occur
    p5d <- max(p5d, na.rm = T)
    if(p5d == 0){
      p5d <- NA
      whn <- NA
    }
  } else {
    p5d <- NA
    whn <- NA
  }
  return(c(p5d,whn))
}

## ------------------------------------------ ##
## Drought
## ------------------------------------------ ##

# Maximum number of consecutive days
cdd <- function(x){ # x: daily precipitation vector
  if(!all(is.na(x))){
    runs <- rle(x < 1)
    cdd  <- max(runs$lengths[runs$values==1], na.rm=TRUE)
  } else {
    cdd  <- NA
  }
  return(cdd)
}

# ------------------------------------------ #
## Heat
# ------------------------------------------ #

# Thermal humidity index
thi <- function(tavg, rhum){
  thi <- (1.8 * tavg + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tavg - 26.8))
  return(thi)
}
