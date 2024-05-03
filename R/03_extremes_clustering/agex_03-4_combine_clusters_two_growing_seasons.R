## ------------------------------------------ ##
## Combine extreme clusters of places with two growing seasons
## By: Harold Achicanoy
## WUR & ABC
## May 2024
## ------------------------------------------ ##

## R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, geodata, Kendall, tidyverse, psych,
                                FactoMineR, factoextra, modifiedmk, NbClust,
                                rnaturalearth, RColorBrewer, MetBrewer,
                                fastcluster, eurostat, giscoR,parallelDist))

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
yrs    <- 1980:2022
index  <- 'spei-6'
gs     <- 'two'

# Extreme clusters for places with one growing seasons
agex_sgn <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/clusters/agex_global_spei-6_one_s1_fmadogram_clean.tif'))
# Extreme clusters for places with two growing seasons
agex_sgn_s1 <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/clusters/agex_global_spei-6_two_s1_fmadogram_k225.tif'))
agex_sgn_s2 <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/clusters/agex_global_spei-6_two_s2_fmadogram_k222.tif'))

# Lower limit: minimum number of pixels per extreme cluster based on gs1
lwr_lmt <- min(table(as.numeric(terra::values(agex_sgn, na.rm = T))))

# Counting the number of coincident extreme clusters
dfm <- terra::as.data.frame(x = c(agex_sgn_s1, agex_sgn_s2), xy = T, cell = T)
names(dfm)[(ncol(dfm)-1):ncol(dfm)] <- c('s1','s2')
dfm$combined <- paste0(dfm$s1,'-',dfm$s2)
table(dfm$combined)

# Remove clusters with a fewer number of pixels
aux <- base::as.data.frame(sort(table(dfm$combined), decreasing = T))
names(aux) <- c('combined','count')
aux$extreme_cluster <- 1:nrow(aux)
aux$extreme_cluster[aux$count < lwr_lmt] <- NA

dfm <- dplyr::left_join(x = dfm, y = aux, by = 'combined')

agex_sgn_aux <- agex_sgn; rm(agex_sgn)
terra::values(agex_sgn_aux) <- NA
agex_sgn_aux[dfm$cell] <- dfm$extreme_cluster

terra::writeRaster(x = agex_sgn_aux, paste0(root,'/agroclimExtremes/agex_results/clusters/agex_global_spei-6_two_fmadogram_clean.tif'))
