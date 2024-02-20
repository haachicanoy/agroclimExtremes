## ------------------------------------------ ##
## Find optimal number of clusters
## By: Harold Achicanoy
## WUR & ABC
## Feb. 2024
## ------------------------------------------ ##

## R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, geodata, Kendall, tidyverse, psych,
                                FactoMineR, factoextra, modifiedmk, NbClust,
                                rnaturalearth, RColorBrewer, MetBrewer,
                                fastcluster, eurostat, giscoR))

## Directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

yrs <- 1979:2022

## List and load files
cdd_roi <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/agex_cdd_25km/one_s1_cdd_25km.tif'))
names(cdd_roi) <- paste0('Y',yrs)

## Apply transformations
cdd_roi_fnl <- terra::app(x = cdd_roi,
                          fun = function(x){
                            if(!all(is.na(x))){
                              # Remove NA and -Inf values
                              x[is.na(x)] <- 0
                              x[is.infinite(x)] <- 0
                              # Scale the time series
                              x_scl <- scale(x)
                              # Mann-Kendall test for evaluating trend significance
                              trd <- Kendall::MannKendall(x_scl)$sl
                              cnd <- ifelse(trd <= 0.05, T, F)
                              if(cnd){ # Remove trend by fitting loess regression
                                dfm <- data.frame(x = yrs, y = x_scl)
                                fit <- loess(y ~ x, data = dfm)
                                x_ntrd <- dfm$y - fit$fitted
                              } else { # Same scaled time series
                                x_ntrd <- x_scl
                              }
                            } else {
                              x_ntrd <- rep(NA, length(x))
                            }
                            return(x_ntrd)
                          })
names(cdd_roi_fnl) <- paste0('Y',yrs); gc(T)

## Transform into a data.frame
cdd_roi_ntrd <- terra::as.data.frame(x = cdd_roi_fnl, xy = T, cell = T, na.rm = T)
# Sample a number of pixels
set.seed(1235)
smp <- sample(x = cdd_roi_ntrd$cell, size = 0.01 * dim(cdd_roi_ntrd)[1], replace = F)
# Sampled locations
x <- t(cdd_roi_ntrd[cdd_roi_ntrd$cell %in% smp,4:ncol(cdd_roi_ntrd)])

## Compute F-madogram distance
get_fmado_dist <- function(x){
  Nnb <- ncol(x)
  Tnb <- nrow(x)
  V <- array(NaN, dim = c(Tnb,Nnb))
  for(p in 1:Nnb){
    x.vec <- as.vector(x[,p])
    if(all(is.na(x.vec))){
      V[,p] <- x.vec
      next
    }
    Femp <- ecdf(x.vec)(x.vec)
    V[,p] <- Femp
  }
  # DD_fmado = dist(t(V), method = "manhattan", diag = TRUE, upper = TRUE)/(2*Tnb)
  DD_fmado <- parallelDist::parDist(x = t(V), method = 'manhattan', diag = T, upper = T)/(2*Tnb)
  return(DD_fmado)
}
# Trim distribution of F-madogram distances greater than 1/6
cap_fmado_dist <- function(DD_fmado){
  DD_fmado_adj <- pmin(DD_fmado, rep(1/6, length(DD_fmado)))
  if(any(DD_fmado_adj > 1/6, na.rm = TRUE) == 1)
    stop("ERROR: cap_fmado_dist didn't work")
  return(DD_fmado_adj)
}

## Obtain distance matrices
fmado_dist <- get_fmado_dist(x); gc(T) # F-madogram distances
# Saving F-madogram distance matrix. Not efficient.
# saveRDS(object = fmado_dist, file = paste0(root,'/agroclimExtremes/agex_results/cdd_extreme_distance.rds'))
fmado_dist <- cap_fmado_dist(fmado_dist); gc(T) # Truncated F-madogram distances

library(cluster)

# Use map_dbl to run many models with varying value of k
sil_width <- purrr::map_dbl(2:200,  function(k){
  model <- cluster::pam(x = fmado_dist, k = k)
  model$silinfo$avg.width
})

# Generate a data frame containing both k and sil_width
sil_df <- data.frame(
  k = 2:200,
  sil_width = sil_width
)

# Plot the relationship between k and sil_width
ggplot(sil_df, aes(x = k, y = sil_width)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 2:200)

pam_cls <- cluster::pam(x = , k = 2, pamonce = 6)
