# ------------------------------------------ #
# Crops classification - croplands
# By: Harold Achicanoy
# WUR & ABC
# Feb. 2024
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'pattern')

# Define directories
root <- '//CATALOGUE/WFP_ClimateRiskPr1'           # Server
inp_dir <- 'D:/Data/Maps'                          # Local
out_dir <- 'D:/OneDrive - CGIAR/PhD/papers/paper1' # local

# Load MapSPAM readme
mapspam_info <- readLines(paste0(root,'/agroclimExtremes/agex_raw_data/croplands/spam2010/ReadMe_v2r0_Global.txt'))

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

# Capital ID
toupper(agr$`SPAM name`)

# Template rasters
tmp_10km <- terra::rast('https://github.com/haachicanoy/agroclimExtremes/raw/main/data/tmp_era5.tif')

## Cropland areas from MapSPAM 2010
crops_dir <- paste0(inp_dir,'/spam2010') # Directory

# Food and non-food groups: 10
grps <- unique(agr$GROUP)

# Count the number of crops per crop class
for(i in 1:length(grps)){
  if(!file.exists(paste0(out_dir,'/data/agex_',gsub(' ','_',grps[i]),'_count.tif'))){
    # List files per group
    crops_fls <- list.files2(path = crops_dir,
                             pattern = paste0(toupper(agr$`SPAM name`[agr$GROUP == grps[i]]),'_A.tif'),
                             full.names = T) |> as.character()
    crops_lnd <- crops_fls |> purrr::map(.f = function(x){
      r <- terra::rast(x); r[r != 0] <- 1; return(r)
    }) |> terra::rast() |> sum()
    terra::writeRaster(x = crops_lnd, filename = paste0(out_dir,'/data/agex_',gsub(' ','_',grps[i]),'_count.tif'), overwrite = T)
    if(!file.exists(paste0(root,'/agroclimExtremes/agex_raw_data/agex_',gsub(' ','_',grps[i]),'_count_10km.tif'))){
      # Resampling MapSPAM into AgERA5 template resolution
      crops_lnd_10km <- terra::resample(x = crops_lnd, y = tmp_10km, method = 'near')
      terra::writeRaster(x = crops_lnd_10km, filename = paste0(root,'/agroclimExtremes/agex_raw_data/agex_',gsub(' ','_',grps[i]),'_count_10km.tif'), overwrite = T)
    }
  }
}
