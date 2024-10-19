# ------------------------------------------ #
# Description of extreme weather clusters
# By: Harold Achicanoy
# WUR & ABC
# May 2024
# ------------------------------------------ #

## R options and packages loading ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,geodata,NbClust,entropy,scales,dplyr,
                                ggplot2,RColorBrewer,landscapemetrics,
                                exactextractr,hrbrthemes,trend,quantreg,
                                OutliersO3,MetBrewer,extrafont))

## Relevant functions ----
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
get_trend <- function(x){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    y <- as.numeric(trend::sens.slope(x)$estimates)
  } else { y <- NA }
  return(y)
}
get_ext_trend <- function(x){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    dfm <- data.frame(time = 1:length(x), ts = x)
    qrf <- quantreg::rq(formula = ts ~ time, tau = .95, data = dfm)
    y <- as.numeric(qrf$coefficients[2])
  } else { y <- NA }
  return(y)
}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Key arguments ----
root   <- '//CATALOGUE/AgroclimExtremes'
index  <- 'spei-6'

## Extreme weather clusters ----
agex_sgn_cln <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_cln) <- 'extreme_cluster'
# Extreme weather clusters by growing seasons
agex_sgn_gs1 <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_one_s1_fmadogram_clean.tif'))
agex_sgn_gs1[!is.na(agex_sgn_gs1)] <- 1
agex_sgn_gs2 <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_two_fmadogram_clean.tif'))
agex_sgn_gs2[!is.na(agex_sgn_gs2)] <- 2
agex_sgn_gs <- terra::merge(agex_sgn_gs1, agex_sgn_gs2); rm(agex_sgn_gs1, agex_sgn_gs2)
names(agex_sgn_gs) <- 'growing_seasons'
# Extreme weather clusters coordinates
crds_sgn_cln <- terra::as.data.frame(x = agex_sgn_cln, xy = T, cell = T, na.rm = T)
aux <- terra::extract(x = agex_sgn_gs, y = crds_sgn_cln[,c('x','y')])
crds_sgn_cln$growing_seasons <- aux$growing_seasons; rm(aux, agex_sgn_gs)

## Global shapefile
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- wrl[,c('adm0_iso','name_en','continent','subregion')]

## Countries collaboration ----
# Countries coverage across extreme drought clusters
out <- paste0(root,'/agex_results/agex_country_collaborations.csv')
if(!file.exists(out)){
  
  # Pixels within countries
  crds_cty <- terra::intersect(x = wrl, y = terra::vect(crds_sgn_cln, c('x','y'), crs = 'EPSG:4326')) |> base::as.data.frame()
  
  # Pixels number per country
  crds_cnt <- crds_cty |>
    dplyr::group_by(adm0_iso, extreme_cluster, growing_seasons) |>
    dplyr::count() |>
    dplyr::arrange(extreme_cluster) |>
    dplyr::ungroup() |>
    base::as.data.frame()
  
  # Unique geographies
  unq_geos <- unique(crds_cty[,c('adm0_iso','name_en','continent','subregion','extreme_cluster','growing_seasons')])
  unq_geos <- base::as.data.frame(dplyr::arrange(.data = unq_geos, extreme_cluster))
  
  # Identify collaboration opportunities among countries
  vls <- sort(unique(unq_geos$extreme_cluster))
  collaborations <- lapply(vls, function(i){
    smm <- data.frame(extreme_cluster = i,
                      growing_seasons = unique(unq_geos$growing_seasons[unq_geos$extreme_cluster == i]),
                      countries_count = length(unq_geos[unq_geos$extreme_cluster == i,'adm0_iso']),
                      isos = paste0(unq_geos[unq_geos$extreme_cluster == i,'adm0_iso'], collapse = ','),
                      countries = paste0(unq_geos[unq_geos$extreme_cluster == i,'name_en'], collapse = ','),
                      continents = paste0(unique(unq_geos[unq_geos$extreme_cluster == i,'continent']), collapse = ',') # ,
                      # regions = paste0(unique(unq_geos[unq_geos$extreme_cluster == i,'subregion']), collapse = ',')
    )
    return(smm)
  }) |> dplyr::bind_rows() |>
    dplyr::arrange(-countries_count) |>
    base::as.data.frame()
  
  utils::write.csv(x = collaborations, file = out, row.names = F)
} else {
  collaborations <- utils::read.csv(out)
}; rm(out)
# Calculate coverage area in km^2
agex_area <- terra::cellSize(x = agex_sgn_cln, unit = 'km') |> terra::mask(mask = agex_sgn_cln)
area_dfm <- terra::zonal(x = agex_area, z = agex_sgn_cln, fun = 'sum') |> base::as.data.frame()
collaborations <- dplyr::left_join(x = collaborations, y = area_dfm, by = 'extreme_cluster'); rm(area_dfm, agex_area)

## Clusters quality ----
# They are computed to provide a quality control of the produced clusters
# Metric: Cohesion index
# Meaning:
# - 0 if patches of class become more isolated
# - 100 if patches of class i become more aggregated
outfile <- paste0(root,'/agex_results/agex_quality_metrics.csv')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  agex_chs <- landscapemetrics::lsm_c_cohesion(landscape = agex_sgn_cln) |> base::as.data.frame()
  agex_chs <- agex_chs[,c('class','value')]
  names(agex_chs) <- c('extreme_cluster','cls_cohesion')
  # Metric: Contiguity index mean
  # Meaning: equals the mean of the contiguity index on class level for all patches
  # - 0 maximally non-contiguous
  # - 1 maximally contiguous
  agex_ctg <- landscapemetrics::lsm_c_contig_mn(landscape = agex_sgn_cln) |> base::as.data.frame()
  agex_ctg <- agex_ctg[,c('class','value')]
  names(agex_ctg) <- c('extreme_cluster','cls_contiguity')
  # Metric: Aggregation index
  # Meaning:
  # - 0 maximally disaggregated
  # - 100 maximally aggregated
  agex_agr <- landscapemetrics::lsm_c_ai(landscape = agex_sgn_cln) |> base::as.data.frame()
  agex_agr <- agex_agr[,c('class','value')]
  names(agex_agr) <- c('extreme_cluster','cls_aggregation')
  
  # Merge quality metrics in one table
  qlt_mtrcs <- dplyr::left_join(x = agex_chs, y = agex_ctg, by = 'extreme_cluster')
  qlt_mtrcs <- dplyr::left_join(x = qlt_mtrcs, y = agex_agr, by = 'extreme_cluster')
  qlt_mtrcs[,-1] <- round(qlt_mtrcs[,-1], 3)
  rm(agex_chs, agex_ctg, agex_agr)
  
  # Run a PCA to produce a quality score per cluster
  agex_qly_pca <- stats::prcomp(x = qlt_mtrcs[,-1], retx = T, center = T, scale. = T)
  qlt_mtrcs$quality_rank <- rank(-1*agex_qly_pca$x[,1])
  utils::write.csv(x = qlt_mtrcs, outfile, row.names = F)
} else {
  qlt_mtrcs <- utils::read.csv(outfile)
}; rm(outfile)

## Crop classes diversity ----
# Metric: Shannon diversity index
# Meaning: High entropy means high variation of food/non-food classes
outfile <- paste0(root,'/agex_results/agex_vulnerability/crop_classes_diversity.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  crp_fls <- list.files(path = paste0(root,'/agex_raw_data'), pattern = '^agex_cropclass_', full.names = T)
  crp_fls <- crp_fls[-grep2(pattern = c('fibres','stimulant','rest_of_crops'), x = crp_fls)] # Categories to exclude
  crp_typ <- terra::rast(crp_fls)
  crp_ntp <- terra::app(x = crp_typ, fun = function(i, ff) ff(i), cores = 20, ff = entropy::entropy)
  names(crp_ntp) <- 'crop_classes_diversity'
  terra::writeRaster(x = crp_ntp, filename = outfile, overwrite = T)
  crp_ntp_25km <- terra::resample(x = crp_ntp, y = agex_sgn_cln, method = 'cubicspline', threads = T)
  crp_ntp_25km <- terra::mask(x = crp_ntp_25km, mask = agex_sgn_cln)
  terra::writeRaster(x = crp_ntp_25km, filename = paste0(root,'/agex_results/agex_vulnerability/crop_classes_diversity_25km.tif'), overwrite = T)
} else {
  crp_ntp <- terra::rast(outfile)
}; rm(outfile)

## Arable lands ----
outfile <- paste0(root,'/agex_results/agex_vulnerability/harvested_area_total.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  crops_fls <- list.files(path = paste0(root,'/agex_raw_data/croplands/spam2010'), pattern = '_A.tif', full.names = T)
  # Crops classification
  grp <- read.csv('https://raw.githubusercontent.com/wri/MAPSPAM/master/metadata_tables/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv')
  grp <- grp[,c('SPAM.short.name','GROUP')]
  foods <- grp[!(grp$GROUP %in% c('fibres','stimulant','')),'SPAM.short.name']
  crops_fls <- crops_fls[grep2(pattern = toupper(foods), x = crops_fls)]
  hrvstd_area <- terra::rast(crops_fls) |> sum(na.rm = T)
  hrvstd_area <- hrvstd_area * 0.01 # in km^2
  names(hrvstd_area) <- 'harvested_area_total'
  terra::writeRaster(x = hrvstd_area, filename = outfile, overwrite = T)
} else {
  hrvstd_area <- terra::rast(outfile)
}; rm(outfile)

## Agricultural economic value ----
# Metric: Value of production
# Meaning: High value of production means high economic value from agriculture
total_vop <- terra::rast(paste0(root,'/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_CROP_A.tif')); names(total_vop) <- 'total_vop'
food_vop  <- terra::rast(paste0(root,'/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_FOOD_A.tif')); names(food_vop) <- 'food_vop'
nonf_vop  <- terra::rast(paste0(root,'/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_NONF_A.tif')); names(nonf_vop) <- 'non-food_vop'

## Livestock diversity ----
# Metric: Livestock equivalent units
# Meaning: grazing equivalent of one adult dairy cow producing 3 000 kg of milk
# annually, without additional concentrated foodstuffs
lvs_dir <- paste0(root,'/agex_raw_data/livestock')
anmls   <- list.dirs(path = lvs_dir, full.names = F, recursive = F)
anmls   <- anmls[-grep('horses',anmls)]
lvstc_fls <- list.files2(path = paste0(lvs_dir,'/',anmls), pattern = '_Da.tif$', full.names = T); rm(lvs_dir, anmls)
names(lvstc_fls) <- NULL
lvstc_cnt <- terra::rast(lvstc_fls)

# Livestock Units ----
# LU: Livestock Unit
# LU = Buffaloes * 1 + Cattle * 1 + Chickens * ((0.007 + 0.014)/2) +
#      Ducks * 0.01 + Goats * 0.1 + Pigs * ((0.5+0.027)/2) + Sheep * 0.1
# Source: https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Livestock_unit_(LSU)
outfile <- paste0(root,'/agex_results/agex_vulnerability/lsu_individuals.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  wghts <- c(1, 1, (0.007 + 0.014)/2, 0.01, 0.1, (0.5+0.027)/2, 0.1)
  lvstc_unt <- lvstc_cnt * wghts
  names(lvstc_unt) <- c('Buffaloes','Cattle','Chickens','Ducks','Goats','Pigs','Sheep')
  terra::writeRaster(lvstc_unt, filename = outfile, overwrite = T)
} else {
  lvstc_unt <- terra::rast(outfile)
}; rm(outfile)

outfile <- paste0(root,'/agex_results/agex_vulnerability/lsu_total.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  lsu <- sum(lvstc_unt)
  rm(lvstc_fls, lvstc_cnt)
  names(lsu) <- 'total_livestock_units'
  terra::writeRaster(x = lsu, filename = outfile, overwrite = T)
  lsu_25km <- terra::resample(x = lsu, y = agex_sgn_cln, method = 'cubicspline', threads = T)
  lsu_25km <- terra::mask(x = lsu_25km, mask = agex_sgn_cln)
  terra::writeRaster(x = lsu_25km, filename = paste0(root,'/agex_results/agex_vulnerability/lsu_total_25km.tif'), overwrite = T)
} else {
  lsu <- terra::rast(outfile)
}; rm(outfile)

outfile <- paste0(root,'/agex_results/agex_vulnerability/lsu_diversity.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  lvs_ntp <- terra::app(x = lvstc_unt, fun = function(i, ff) ff(i), cores = 20, ff = entropy::entropy)
  names(lvs_ntp) <- 'livestock_units_diversity'
  terra::writeRaster(x = lvs_ntp, filename = outfile, overwrite = T)
  lvs_ntp_25km <- terra::resample(x = lvs_ntp, y = agex_sgn_cln, method = 'cubicspline', threads = T)
  lvs_ntp_25km <- terra::mask(x = lvs_ntp_25km, mask = agex_sgn_cln)
  terra::writeRaster(x = lvs_ntp_25km, filename = paste0(root,'/agex_results/agex_vulnerability/lsu_diversity_25km.tif'), overwrite = T)
} else {
  lvs_ntp <- terra::rast(outfile)
}; rm(outfile)


## Exposed population ----
# Metric: Population total
# Meaning: number of people per pixel
pop <- terra::rast(paste0(root,'/agex_raw_data/gpw_v4_population_count_rev11_2020_15_min.tif'))
names(pop) <- 'total_population'

## Index severity ----
# Metric: SPEI-6 trend for 90th percentile. We also computed the number of years when SPEI < -1.5
# For two-growing-seasons pixels we computed the average of the trends
# Meaning: extreme drought trend over time
# Loading SPEI-6 time series
stp <- data.frame(gs = c('one','two','two'), season = c('1','1','2'))
idx_one_s1 <- terra::rast(paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/',stp$gs[1],'_s',stp$season[1],'_',index,'_25km.tif'))
idx_two_s1 <- terra::rast(paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/',stp$gs[2],'_s',stp$season[2],'_',index,'_25km.tif'))
idx_two_s2 <- terra::rast(paste0(root,'/agex_indices/agex_',index,'/agex_',index,'_25km/',stp$gs[3],'_s',stp$season[3],'_',index,'_25km.tif'))

# SPEI-6 average
outfile <- paste0(root,'/agex_results/agex_vulnerability/SPEI-6_average.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  idx_one_s1_avg <- mean(idx_one_s1)
  idx_two_s1_avg <- mean(idx_two_s1)
  idx_two_s2_avg <- mean(idx_two_s2)
  idx_two_avg <- mean(c(idx_two_s1_avg,idx_two_s2_avg)); rm(idx_two_s1_avg, idx_two_s2_avg)
  idx_avg <- terra::merge(idx_one_s1_avg, idx_two_avg)
  names(idx_avg) <- 'SPEI-6_average'
  terra::writeRaster(x = idx_avg, filename = outfile, overwrite = T)
} else {
  idx_avg <- terra::rast(outfile)
}; rm(outfile)

# SPEI-6 range
outfile <- paste0(root,'/agex_results/agex_vulnerability/SPEI-6_range.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  idx_one_s1_range <- range(idx_one_s1)
  idx_two_s1_range <- range(idx_two_s1)
  idx_two_s2_range <- range(idx_two_s2)
  idx_two_range <- c(mean(idx_two_s1_range[[1]],idx_two_s2_range[[1]]),
                     mean(idx_two_s1_range[[2]],idx_two_s2_range[[2]])); rm(idx_two_s1_range, idx_two_s2_range)
  idx_range <- terra::merge(idx_one_s1_range, idx_two_range); rm(idx_one_s1_range, idx_two_range)
  names(idx_range) <- c('SPEI-6_min', 'SPEI-6_max')
  terra::writeRaster(x = idx_range, filename = outfile, overwrite = T)
} else {
  idx_range <- terra::rast(outfile)
}; rm(outfile)

# # Trend
# idx_one_s1_slp <- terra::app(x = idx_one_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_trend); names(idx_one_s1_slp) <- 'SPEI-6_slope'
# idx_two_s1_slp <- terra::app(x = idx_two_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_trend); names(idx_two_s1_slp) <- 'SPEI-6_slope'
# idx_two_s2_slp <- terra::app(x = idx_two_s2, fun = function(i, ff) ff(i), cores = 20, ff = get_trend); names(idx_two_s2_slp) <- 'SPEI-6_slope'

# Extreme trend
outfile <- paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_trend.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  # Calculating extreme trend (95th percentile) for one-growing-season places
  idx_one_s1_sxt <- terra::app(x = idx_one_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_ext_trend); names(idx_one_s1_sxt) <- 'SPEI-6_slope_95th'
  # Calculating extreme trend (95th percentile) for two-growing-season (S1) places
  idx_two_s1_sxt <- terra::app(x = idx_two_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_ext_trend); names(idx_two_s1_sxt) <- 'SPEI-6_slope_95th'
  # Calculating extreme trend (95th percentile) for two-growing-season (S2) places
  idx_two_s2_sxt <- terra::app(x = idx_two_s2, fun = function(i, ff) ff(i), cores = 20, ff = get_ext_trend); names(idx_two_s2_sxt) <- 'SPEI-6_slope_95th'
  # Computing the average of two-growing-season places
  idx_two_sxt <- mean(c(idx_two_s1_sxt, idx_two_s2_sxt))
  # Merge results from one and two growing-season places
  idx_sxt <- terra::merge(idx_one_s1_sxt, idx_two_sxt)
  terra::writeRaster(x = idx_sxt, filename = outfile, overwrite = T)
} else {
  idx_sxt <- terra::rast(outfile)
}; rm(outfile)

# Extreme drought years' count
outfile <- paste0(root,'/agex_results/agex_vulnerability/SPEI-6_extreme_years_count.tif')
if(!file.exists(outfile)){
  dir.create(path = dirname(outfile), F, T)
  # Calculating extreme years count (SPEI-6 < -1.5) for one-growing-season places
  idx_one_s1_cnt <- terra::app(x = idx_one_s1, fun = function(i, ff) ff(i), cores = 20, ff = function(x){sum(x > 1.5)}); names(idx_one_s1_cnt) <- 'SPEI-6_extreme_drought_count'
  # Calculating extreme years count (SPEI-6 < -1.5) for two-growing-season (S1) places
  idx_two_s1_cnt <- terra::app(x = idx_two_s1, fun = function(i, ff) ff(i), cores = 20, ff = function(x){sum(x > 1.5)}); names(idx_two_s1_cnt) <- 'SPEI-6_extreme_drought_count'
  # Calculating extreme years count (SPEI-6 < -1.5) for two-growing-season (S2) places
  idx_two_s2_cnt <- terra::app(x = idx_two_s2, fun = function(i, ff) ff(i), cores = 20, ff = function(x){sum(x > 1.5)}); names(idx_two_s2_cnt) <- 'SPEI-6_extreme_drought_count'
  # Computing the average of two-growing-season places
  idx_two_cnt <- round(mean(c(idx_two_s1_cnt, idx_two_s2_cnt)))
  # Merge results from one and two growing-season places
  idx_cnt <- terra::merge(idx_one_s1_cnt, idx_two_cnt)
  terra::writeRaster(x = idx_cnt, filename = outfile, overwrite = T)
} else {
  idx_cnt <- terra::rast(outfile)
}; rm(outfile)

# Transform extreme drought clusters into shapefile
agex_sgn_poly <- terra::as.polygons(x = agex_sgn_cln)

## Save all metrics ----
rst_mtrcs <- cbind(
  terra::zonal(x = pop, z = agex_sgn_poly, fun = 'sum', na.rm = T),
  terra::zonal(x = c(total_vop, food_vop, nonf_vop), z = agex_sgn_poly, fun = 'sum', na.rm = T),
  terra::zonal(x = hrvstd_area, z = agex_sgn_poly, fun = 'sum', na.rm = T),
  terra::zonal(x = hrvstd_area, z = agex_sgn_poly, fun = 'mean', na.rm = T, exact = T) |> dplyr::rename(harvested_area_average = 'harvested_area_total'),
  terra::zonal(x = hrvstd_area, z = agex_sgn_poly, fun = 'min', na.rm = T, exact = T) |> dplyr::rename(harvested_area_min = 'harvested_area_total'),
  terra::zonal(x = hrvstd_area, z = agex_sgn_poly, fun = 'max', na.rm = T, exact = T) |> dplyr::rename(harvested_area_max = 'harvested_area_total'),
  terra::zonal(x = crp_ntp, z = agex_sgn_poly, fun = 'mean', na.rm = T),
  terra::zonal(x = lsu, z = agex_sgn_poly, fun = 'sum', na.rm = T) |> dplyr::rename(livestock_units_total = 'total_livestock_units'),
  terra::zonal(x = lsu, z = agex_sgn_poly, fun = 'mean', na.rm = T, exact = T) |> dplyr::rename(livestock_units_average = 'total_livestock_units'),
  terra::zonal(x = lsu, z = agex_sgn_poly, fun = 'min', na.rm = T, exact = T) |> dplyr::rename(livestock_units_min = 'total_livestock_units'),
  terra::zonal(x = lsu, z = agex_sgn_poly, fun = 'max', na.rm = T, exact = T) |> dplyr::rename(livestock_units_max = 'total_livestock_units'),
  terra::zonal(x = c(idx_avg, idx_range, idx_sxt), z = agex_sgn_poly, fun = 'mean', na.rm = T)
)
rst_mtrcs$extreme_cluster <- 1:nrow(rst_mtrcs)
rst_mtrcs <- rst_mtrcs[,c('extreme_cluster',names(rst_mtrcs)[-ncol(rst_mtrcs)])]
head(rst_mtrcs,3)
rst_mtrcs[,-1] <- round(rst_mtrcs[,-1], 4)

dfm <- dplyr::left_join(x = collaborations, y = rst_mtrcs, by = 'extreme_cluster')
dfm <- dplyr::left_join(x = dfm, y = qlt_mtrcs, by = 'extreme_cluster')

hist(dfm$`SPEI-6_slope_95th`)
abline(v = quantile(x = dfm$`SPEI-6_slope_95th`, probs = seq(0,1,1/4)), lty = 2, col = 'red')
c('Negative','Low','Medium','High')

hist(dfm$crop_classes_diversity)
abline(v = quantile(x = dfm$crop_classes_diversity, probs = seq(0,1,1/4)), lty = 2, col = 'red')
c('Very low','Low','Medium','High')

hist(dfm$livestock_units_total)
abline(v = quantile(x = dfm$livestock_units_total, probs = seq(0,1,1/4)), lty = 2, col = 'red')
c('Very low','Low','Medium','High')

agex_sgn_poly <- terra::merge(x = agex_sgn_poly, y = dfm)

terra::writeVector(x = agex_sgn_poly, filename = paste0(root,'/agex_results/agex_results_clusters/vct_agex_global_',index,'.gpkg'), overwrite = T)
utils::write.csv(x = dfm, paste0(root,'/agex_results/agex_all_metrics.csv'), row.names = F)
