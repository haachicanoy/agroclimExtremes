# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Crop classification counts at 0.1° on agricultural production areas
# By: Harold Achicanoy
# WUR & ABC
# Created in February 2024
# Modified in August 2026
# --------------------------------------------------------------- #

# R options and user-defined functions
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra))
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'pattern')

# Define directories
root <- '//CATALOGUE/AgroclimExtremes'
inp_dir <- paste0(root,'/agex_raw_data/croplands')
if (!file.exists(paste0(root,"/agex_raw_data/croplands/MapSPAM_crops_classification.rds"))) {
  
  # Load MapSPAM readme
  mapspam_info <- readLines(paste0(root,'/agex_raw_data/croplands/spam2010/ReadMe_v2r0_Global.txt'))
  
  # Food crops
  crps <- strsplit(x = mapspam_info[105:132], split = '\t', fixed = T)
  nms <- crps[[1]][!(nchar(crps[[1]]) == 0)]
  crps <- lapply(2:length(crps), function(i){
    dfm <- data.frame(t(crps[[i]][!(nchar(crps[[i]]) == 0)]))
    names(dfm) <- nms
    return(dfm)
  }) |> dplyr::bind_rows()
  
  # Non-food crops
  ncrps <- strsplit(x = mapspam_info[136:151], split = '\t', fixed = T)
  ncrps[[10]] <- c('31','other fibre crops', 'ofib')
  ncrps <- lapply(2:length(ncrps), function(i){
    dfm <- data.frame(t(ncrps[[i]][!(nchar(ncrps[[i]]) == 0)]))
    names(dfm) <- nms
    return(dfm)
  }) |> dplyr::bind_rows()
  
  # Put them together
  agr <- rbind(crps, ncrps)
  agr$`crop #` <- as.numeric(agr$`crop #`)
  agr <- agr |> dplyr::arrange(`crop #`) |> base::as.data.frame()
  
  # Crops classification
  grp <- read.csv('https://raw.githubusercontent.com/wri/MAPSPAM/master/metadata_tables/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv')
  grp <- grp[,c('SPAM.short.name','GROUP')]
  
  # Add classification column
  agr <- dplyr::left_join(x = agr, y = grp, by = c('SPAM name' = 'SPAM.short.name'))
  agr$GROUP[agr$`SPAM name` == 'rest'] <- 'rest of crops'
  rm(crps, ncrps, grp, mapspam_info, nms)
  
  utils::write.csv(x = agr, file = paste0(root,"/agex_raw_data/croplands/MapSPAM_crops_classification.csv"), row.names = F)
  saveRDS(object = agr, file = paste0(root,"/agex_raw_data/croplands/MapSPAM_crops_classification.rds"))
  
} else {
  agr <- readRDS(paste0(root,"/agex_raw_data/croplands/MapSPAM_crops_classification.rds"))
}
# Make the comparison against v2020
if (!file.exists(paste0(root,"/agex_raw_data/croplands/MapSPAM_crops_classification_v2020.rds"))) {
  
  # Load MapSPAM v2020 readme
  mapspam_info <- readLines(paste0(root,'/agex_raw_data/croplands/spam2020/Readme_SPAM2020V2r2.txt'))
  
  # All MapSPAM crops
  crps <- strsplit(x = mapspam_info[116:162], split = '\t', fixed = T)
  nms <- c("ID", "SPAM name 2020", "Crop name 2020")
  crps <- lapply(2:length(crps), function(i){
    dfm <- data.frame(t(crps[[i]][!(nchar(crps[[i]]) == 0)]))
    names(dfm) <- nms
    return(dfm)
  }) |> dplyr::bind_rows()
  crps
  crps$`SPAM name` <- crps$`SPAM name 2020`
  crps$`SPAM name`[5] <- "smil"
  crps$`SPAM name`[32] <- "acof"
  
  agr <- dplyr::left_join(x = agr, y = crps, by = "SPAM name")
  agr$ID <- NULL
  
  saveRDS(object = agr, file = paste0(root,"/agex_raw_data/croplands/MapSPAM_crops_classification_v2020.rds"))
} else {
  agr <- readRDS(paste0(root,"/agex_raw_data/croplands/MapSPAM_crops_classification_v2020.rds"))
}
agr <- agr |> dplyr::select(`SPAM name`, `SPAM name 2020`, name, `Crop name 2020`, GROUP)

# Capital ID
toupper(agr$`SPAM name`)

# Template rasters
tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')

## Cropland areas from MapSPAM 2010
agr_dts <- "spam2020"
crops_dir <- paste0(inp_dir,'/',agr_dts) # Directory

names(agr) <- c("spam2010","spam2020","name2010","name2020","GROUP")

# Food and non-food groups: 10
grps <- unique(agr$GROUP)

# Count the number of crops per crop class
for(i in 1:length(grps)){
  if(!file.exists(paste0(root,'/agex_raw_data/',agr_dts,'_counts/agex_cropclass_',gsub(' ','_',grps[i]),'_count.tif'))){
    # List files per group
    crops_fls <- list.files2(path = crops_dir,
                             pattern = dplyr::case_when(agr_dts == "spam2010" ~ paste0(toupper(agr$spam2010[agr$GROUP == grps[i]]),'_A.tif'),
                                                        agr_dts == "spam2020" ~ paste0(toupper(agr$spam2020[agr$GROUP == grps[i]]),'_A.tif')),
                             full.names = T) |> as.character()
    crops_lnd <- crops_fls |> purrr::map(.f = function(x){
      r <- terra::rast(x); r <- terra::classify(r, rcl = cbind(0, Inf, 1)); return(r)
    }) |> terra::rast() |> sum(na.rm = T)
    terra::writeRaster(x = crops_lnd, filename = paste0(root,'/agex_raw_data/',agr_dts,'_counts/agex_cropclass_',gsub(' ','_',grps[i]),'_count.tif'), overwrite = T)
    if(!file.exists(paste0(root,'/agex_raw_data/',agr_dts,'_counts/agex_cropclass_',gsub(' ','_',grps[i]),'_count_10km.tif'))){
      # Resampling MapSPAM into AgERA5 template resolution
      crops_lnd_10km <- terra::resample(x = crops_lnd, y = tmp_10km, method = 'near')
      terra::writeRaster(x = crops_lnd_10km, filename = paste0(root,'/agex_raw_data/',agr_dts,'_counts/agex_cropclass_',gsub(' ','_',grps[i]),'_count_10km.tif'), overwrite = T)
    }
  }
}
