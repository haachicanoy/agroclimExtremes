# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Cleaning outliers in places with one growing season
# By: Harold Achicanoy
# WUR & ABC
# Created in May 2024
# Modified in February 2026
# --------------------------------------------------------------- #

# R options and user-defined functions
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,tibble,OutliersO3))
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Key arguments
root   <- '//CATALOGUE/AgroclimExtremes'
index  <- 'spei-6'
gs     <- 'one'
season <- 1

# Load extreme weather clusters
fls <- list.files(path = paste0(root,'/agex_results/clusters'), pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_k[0-9][0-9][0-9].tif$'), full.names = T)
agex_sgn <- terra::rast(fls)
names(agex_sgn) <- 'extreme_cluster'

## CLEANING PROCESS
if(!file.exists(gsub('_k[0-9][0-9][0-9].tif$','_clean.tif',fls))){
  
  # -------------------------------------------------------------------------- #
  # First cleaning phase: Clean clusters from random miss-classification issues
  # -------------------------------------------------------------------------- #
  
  ## Extreme clusters to remove
  # clt_to_remove <- c(1,21,151)
  
  # Cleaning
  agex_sgn_cln <- agex_sgn
  agex_sgn_cln[agex_sgn_cln == 1] <- NA
  agex_sgn_cln[agex_sgn_cln == 21] <- NA
  agex_sgn_cln[agex_sgn_cln == 151] <- NA
  for(cl in c(2:20,22:150,152:299)){
    
    cat(paste0('Processing CLUSTER: ',cl,'\n'))
    # Cluster filtering # cl <- 292 # 188, 209
    aux <- agex_sgn
    aux[aux != cl] <- NA
    
    # Get cluster's coordinates
    aux_dfm <- terra::as.data.frame(x = aux, xy = T, cell = T, na.rm = T)
    
    # Get outliers from coordinates
    out_pre <- OutliersO3::O3prep(data = aux_dfm[,c('x','y')], method = c('HDo','BAC','adjOut','DDC','MCD')) # 'PCS'
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
  
  # -------------------------------------------------------------------------- #
  # Second cleaning phase: Split clusters that have clear structures (manually)
  # -------------------------------------------------------------------------- #
  
  cln_stp <- data.frame(extreme_cluster = c(9,154,131,123,70), tentative_clusters = c(3,2,2,2,2))
  cln_stp <- tibble::as.tibble(cln_stp)
  cln_stp$extreme_values <- list(c(9,300,301),c(154,302),c(131,303),c(123,304),c(70,305))
  for(i in 1:nrow(cln_stp)){
    
    # Filtering to extreme cluster of interest
    aux <- agex_sgn_cln
    aux[aux != cln_stp$extreme_cluster[i]] <- NA
    # Get geographical coordinates
    crds_aux <- terra::as.data.frame(x = aux, xy = T, cell = T, na.rm = T)
    # Fix a random seed and execute k-means clustering
    set.seed(1235)
    km <- stats::kmeans(x = crds_aux[,c('x','y')], centers = cln_stp$tentative_clusters[i])
    replace(x = km$cluster, list = c(cln_stp$extreme_cluster[i],300,301), values = 1:3)
    # Cluster values
    cls_values <- km$cluster
    # Reference extreme cluster values
    ref_values <- cln_stp$extreme_values[[i]]
    # Original k-means values
    kms_values <- unique(cls_values)
    
    upt_values <- c(ref_values, cls_values)[match(cls_values, c(kms_values, cls_values))]
    
    crds_aux$kmeans_cluster <- upt_values
    agex_sgn_cln[crds_aux$cell] <- crds_aux$kmeans_cluster
    
  }
  
  # -------------------------------------------------------------------------- #
  # Third cleaning phase: cluster 9 correction
  # -------------------------------------------------------------------------- #
  
  aux <- agex_sgn_cln
  aux[aux != 9] <- NA
  
  # Get cluster's coordinates
  aux_dfm <- terra::as.data.frame(x = aux, xy = T, cell = T, na.rm = T)
  
  # Get outliers from coordinates
  out_pre <- OutliersO3::O3prep(data = aux_dfm[,c('x','y')], method = c('HDo','BAC','adjOut','DDC','MCD')) # 'PCS'
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
  
  terra::writeRaster(x = agex_sgn_cln, filename = gsub('_k[0-9][0-9][0-9].tif$','_clean.tif',fls), overwrite = T)
  
} else {
  agex_sgn_cln <- terra::rast(gsub('_k[0-9][0-9][0-9].tif$','_clean.tif',fls))
}
