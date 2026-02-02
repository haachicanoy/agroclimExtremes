# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Combine extreme clusters of places with two growing seasons
# By: Harold Achicanoy
# WUR & ABC
# Created in May 2024
# Modified in February 2026
# --------------------------------------------------------------- #

# R options and user-defined functions
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, geodata, Kendall, tidyverse, psych,
                                FactoMineR, factoextra, modifiedmk, NbClust,
                                rnaturalearth, RColorBrewer, MetBrewer,
                                fastcluster, eurostat, giscoR,parallelDist,
                                OutliersO3))
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Key arguments
root   <- '//CATALOGUE/AgroclimExtremes'
yrs    <- 1980:2022
index  <- 'spei-6'
gs     <- 'two'

# Extreme clusters for places with one growing seasons
agex_sgn <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_one_s1_fmadogram_clean.tif'))
# Extreme clusters for places with two growing seasons
agex_sgn_s1 <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_two_s1_fmadogram_k225.tif'))
agex_sgn_s2 <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_two_s2_fmadogram_k222.tif'))

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

outfile <- paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_two_fmadogram_clean.tif')

## CLEANING PROCESS
if(!file.exists(outfile)){
  
  # -------------------------------------------------------------------------- #
  # First cleaning phase: Clean clusters from random miss-classification issues
  # -------------------------------------------------------------------------- #
  
  ## Extreme clusters to remove
  
  # Cleaning
  agex_sgn_cln <- agex_sgn_aux
  n_clusters   <- length(unique(as.numeric(terra::values(agex_sgn_aux, na.rm = T))))
  for(cl in 1:n_clusters){
    
    cat(paste0('Processing CLUSTER: ',cl,'\n'))
    aux <- agex_sgn_aux
    aux[aux != cl] <- NA
    
    # Get cluster's coordinates
    aux_dfm <- terra::as.data.frame(x = aux, xy = T, cell = T, na.rm = T)
    
    # Get outliers from coordinates
    set.seed(1235)
    out_pre <- OutliersO3::O3prep(data = aux_dfm[,c('x','y')], method = c('HDo','BAC','adjOut','DDC','MCD')) # 'PCS'
    out_msg <- capture.output(try(OutliersO3::O3plotM(out_pre), silent = T))
    if(length(out_msg) > 1){
      out_res <- OutliersO3::O3plotM(out_pre)
      
      # Number of outliers per method
      out_cns <- out_res$nOut
      print(out_cns)
      n_out <- getmode(out_cns)
      
      cat(paste0('Outliers identified: ',n_out,'\n\n'))
      if(n_out > 0){
        out_cnd <- out_res$outsTable$Method %in% names(which(out_res$nOut == n_out))
        out_css <- unique(out_res$outsTable[out_cnd,'Case'])
        agex_sgn_cln[aux_dfm[out_css,'cell']] <- NA
      }
    }
    
  }
  
  terra::writeRaster(x = agex_sgn_cln, outfile, overwrite = T)
  
} else {
  agex_sgn_cln <- terra::rast(outfile)
}
