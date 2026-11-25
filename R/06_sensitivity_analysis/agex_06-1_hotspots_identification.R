# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Sensitivity analysis: hotspot identification for any combination
# of drought index (SPEI-3/6/12) and agricultural dataset
# (MapSPAM 2010/2020)
# By: Harold Achicanoy
# WUR & ABC
# Created in November 2026
# --------------------------------------------------------------- #
# Outputs (per scenario) in agex_results/agex_sensitivity:
# - agex_hotspots_<scenario>.csv            cluster metrics, 4x4 classes, hotspot flags
# - agex_hotspots_thresholds_<scenario>.csv quartile thresholds (answers R3.7)
# - agex_hotspots_<scenario>.tif            pixel layers used in 06-3
# - agex_hotspots_summary.csv               counts of hotspots per scenario
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, dplyr))
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/refs/heads/main/R/06_sensitivity_analysis/agex_06-0_sensitivity_functions.R')
# source('./R/06_sensitivity/agex_06-0_sensitivity_functions.R') # local alternative

# Key arguments ----
root    <- '//CATALOGUE/AgroclimExtremes'
out_dir <- paste0(root,'/agex_results/agex_sensitivity')
indices <- c('spei-3','spei-6','spei-12')
agr_dss <- c('spam2010','spam2020')

# Crop value of production for MapSPAM 2020 (25 km) ----
# The 2010 layer (food_vop.tif) was produced for the main analysis. Build the
# 2020 equivalent with the same steps. CHECK the raw SPAM2020 VoP file name.
outfile <- paste0(root,'/agex_results/agex_vulnerability/food_vop_spam2020.tif')
if(!file.exists(outfile)){
  tmp_25km <- terra::rast(paste0(root,'/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
  vop_2020 <- terra::rast(paste0(root,'/agex_raw_data/agex_spam2020V2r2_global_V_agg_VP_FOOD_A.tif')) # CHECK
  vop_2020_25km <- terra::resample(x = vop_2020, y = tmp_25km, method = 'cubicspline', threads = T)
  vop_2020_25km <- terra::mask(x = vop_2020_25km, mask = tmp_25km)
  names(vop_2020_25km) <- 'food_vop'
  terra::writeRaster(x = vop_2020_25km, filename = outfile, overwrite = T)
  rm(tmp_25km, vop_2020, vop_2020_25km)
}; rm(outfile)

# Scenarios ----
# Baseline = spei-6 x spam2010 (main manuscript)
scn <- expand.grid(index = indices, agr_dts = agr_dss, stringsAsFactors = F)
scn$scenario <- paste0(scn$index,'_',scn$agr_dts,'_glad00')
scn$cls_file <- paste0(root,'/agex_results/agex_results_clusters/agex_global_',scn$index,'_combined_fmadogram_clean.tif')
scn

# Run hotspot identification ----
res <- lapply(1:nrow(scn), function(i){
  cat(paste0('Processing scenario: ',scn$scenario[i],'\n'))
  run_hotspots(root     = root,
               index    = scn$index[i],
               agr_dts  = scn$agr_dts[i],
               cls_file = scn$cls_file[i],
               scenario = scn$scenario[i],
               out_dir  = out_dir)
})
names(res) <- scn$scenario

# Summary of hotspots per scenario ----
smm <- lapply(res, function(r){
  cl <- r$classes
  data.frame(scenario            = unique(cl$scenario),
             n_clusters          = nrow(cl),
             crop_diversity      = sum(cl$crop_diversity_hotspot),
             livestock_diversity = sum(cl$livestock_diversity_hotspot),
             joint_diversity     = sum(cl$crop_diversity_hotspot & cl$livestock_diversity_hotspot),
             food_production     = sum(cl$food_production_hotspot),
             crop_vop            = sum(cl$crop_vop_hotspot),
             livestock_units     = sum(cl$livestock_units_hotspot),
             joint_economic      = sum(cl$crop_vop_hotspot & cl$livestock_units_hotspot),
             economic_value      = sum(cl$economic_value_hotspot))
}) |> dplyr::bind_rows()
utils::write.csv(x = smm, file = paste0(out_dir,'/agex_hotspots_summary.csv'), row.names = F)
smm

# Quartile thresholds of the baseline (Table for R3.7) ----
thr <- utils::read.csv(paste0(out_dir,'/agex_hotspots_thresholds_spei-6_spam2010_glad00.csv'))
thr[,2:6] <- round(thr[,2:6], 4)
utils::write.csv(x = thr, file = paste0(out_dir,'/agex_TableS_quartile_thresholds_baseline.csv'), row.names = F)
thr
