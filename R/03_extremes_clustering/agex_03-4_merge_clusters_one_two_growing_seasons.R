## ------------------------------------------ ##
## Merge extreme clusters from one and two growing seasons
## By: Harold Achicanoy
## WUR & ABC
## May 2024
## ------------------------------------------ ##

## R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, dplyr))

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'

agex_sgn_gs1 <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_one_s1_fmadogram_clean.tif'))
agex_sgn_gs2 <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_two_fmadogram_clean.tif'))

combined <- c(agex_sgn_gs1, agex_sgn_gs2) |>
  terra::as.data.frame(xy = T, cell = T)
names(combined)[(ncol(combined)-1):ncol(combined)] <- c('gs1','gs2')
combined$aux <- paste0(combined$gs1,'-',combined$gs2)

dfm <- data.frame(aux = unique(combined$aux))
dfm$extreme_cluster <- 1:nrow(dfm)

combined <- dplyr::left_join(x = combined, y = dfm, by = 'aux')

agex_sgn <- agex_sgn_gs1
terra::values(agex_sgn) <- NA
agex_sgn[combined$cell] <- combined$extreme_cluster

# Manual cleaning
aux <- agex_sgn
aux[aux != 275] <- NA
aux_dfm <- terra::as.data.frame(aux, xy = T, cell = T)
aux_dfm$extreme_cluster[(nrow(aux_dfm)-1):nrow(aux_dfm)] <- NA
agex_sgn[aux_dfm$cell[is.na(aux_dfm$extreme_cluster)]] <- NA

aux <- agex_sgn
aux[aux != 321] <- NA
aux_dfm <- terra::as.data.frame(aux, xy = T, cell = T)
aux_dfm$extreme_cluster[(nrow(aux_dfm)-3):nrow(aux_dfm)] <- NA
agex_sgn[aux_dfm$cell[is.na(aux_dfm$extreme_cluster)]] <- NA

terra::writeRaster(x = agex_sgn, filename = paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'), overwrite = T)
