# --------------------------------------------------------------- #
# Global hotspots of co-occurring extreme droughts in agriculture
# Sensitivity analysis: agreement metrics between scenarios
# (Reviewers 1, 2.1, 2.2, 3.2)
# By: Harold Achicanoy
# WUR & ABC
# Created in November 2026
# --------------------------------------------------------------- #
# For every comparison (baseline vs alternative) it reports:
# - Partition similarity: ARI and NMI on common pixels, number of clusters
# - Hotspot agreement: number of hotspot clusters, hotspot area, area-based
#   Jaccard index, % of baseline hotspot area retained, % of baseline hotspot
#   clusters retained (>= 50% of their area still classified as hotspot)
# - Spearman correlation of the continuous scores (trend and exposure), both
#   pixel-level (spearman_px_*) and cluster-level (spearman_cl_*). High rho with low Jaccard = instability from the quartile
#   thresholds, not from the underlying signal.
# - Regions (countries/continents) where hotspot area is gained or lost
# --------------------------------------------------------------- #

# R options and user-defined functions ----
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra, dplyr, tidyr, rnaturalearth, sf))
source('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/refs/heads/main/R/06_sensitivity_analysis/agex_06-0_sensitivity_functions.R')

# Key arguments ----
root     <- '//CATALOGUE/AgroclimExtremes'
out_dir  <- paste0(root,'/agex_results/agex_sensitivity')
baseline <- 'spei-6_spam2010_glad00'

# Comparisons ----
# Group labels distinguish dataset robustness from drought-indicator robustness.
# MapSPAM 2010 vs 2020 is a change of DATASET (model and inputs), not a
# temporal analysis (producers advise against comparing SPAM versions over time).
cmp <- data.frame(
  scenario_a = c(rep(baseline, 5), rep('spei-6_spam2010_gladref', 4), baseline),
  scenario_b = c('spei-3_spam2010_glad00', 'spei-12_spam2010_glad00',
                 'spei-6_spam2020_glad00',
                 'spei-3_spam2020_glad00', 'spei-12_spam2020_glad00',
                 'spei-6_spam2010_glad01', 'spei-6_spam2010_glad05',
                 'spei-6_spam2010_glad10', 'spei-6_spam2010_glad20',
                 'spei-6_spam2010_gladref'),
  group = c('Drought index','Drought index',
            'Agricultural dataset',
            'Index and dataset','Index and dataset',
            rep('Cropland mask (GLAD)', 4),
            'Simplified pipeline (no outlier cleaning)'),
  stringsAsFactors = F)
cmp <- cmp[file.exists(paste0(out_dir,'/agex_hotspots_',cmp$scenario_b,'.tif')),]
cmp

# Countries ----
ctr <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
ctr <- ctr[,c('adm0_iso','name_en','continent')]
ctr_lbl <- base::as.data.frame(ctr)

# Run comparisons ----
res <- lapply(1:nrow(cmp), function(i){
  cat(paste0('Comparing: ',cmp$scenario_a[i],' vs ',cmp$scenario_b[i],'\n'))
  a <- terra::rast(paste0(out_dir,'/agex_hotspots_',cmp$scenario_a[i],'.tif'))
  b <- terra::rast(paste0(out_dir,'/agex_hotspots_',cmp$scenario_b[i],'.tif'))
  r <- compare_scenarios(a, b, cmp$scenario_a[i], cmp$scenario_b[i], ctr = ctr)
  r$group <- cmp$group[i]
  return(r)
})

# Agreement table (Supplementary Table) ----
agr_tbl <- dplyr::bind_rows(res)
agr_tbl <- agr_tbl[,c('group', names(agr_tbl)[names(agr_tbl) != 'group'])]
num <- sapply(agr_tbl, is.numeric); agr_tbl[,num] <- round(agr_tbl[,num], 3)
utils::write.csv(x = agr_tbl, file = paste0(out_dir,'/agex_TableS_sensitivity_agreement.csv'), row.names = F)

# Compact version for the main text / RtR
agr_tbl |>
  dplyr::filter(hotspot %in% c('crop_diversity_hotspot','livestock_diversity_hotspot',
                               'crop_vop_hotspot','livestock_units_hotspot')) |>
  dplyr::select(group, scenario_b, hotspot, ARI, n_hotspots_a, n_hotspots_b,
                jaccard_area, overlap_retained_pct, clusters_retained_pct,
                spearman_cl_extreme_trend, spearman_cl_crop_diversity, spearman_cl_livestock_diversity,
                spearman_cl_crop_vop, spearman_cl_livestock_units) |>
  utils::write.csv(file = paste0(out_dir,'/agex_sensitivity_agreement_compact.csv'), row.names = F)

# Livestock-focused summary for SPEI-12 (R2.1, R3.2) ----
agr_tbl |>
  dplyr::filter(scenario_b == 'spei-12_spam2010_glad00',
                hotspot %in% c('livestock_diversity_hotspot','livestock_units_hotspot')) |>
  utils::write.csv(file = paste0(out_dir,'/agex_sensitivity_spei12_livestock.csv'), row.names = F)

# Regions where hotspot classification changes ----
reg <- lapply(res, function(r) attr(r, 'regions')) |> dplyr::bind_rows()
reg <- dplyr::left_join(reg, ctr_lbl, by = 'adm0_iso')
utils::write.csv(x = reg, file = paste0(out_dir,'/agex_sensitivity_regions_change_country.csv'), row.names = F)
reg_cnt <- reg |>
  dplyr::group_by(scenario_b, hotspot, continent) |>
  dplyr::summarise(lost_km2 = sum(lost_km2), gained_km2 = sum(gained_km2), .groups = 'drop') |>
  dplyr::arrange(scenario_b, hotspot, -(lost_km2 + gained_km2))
utils::write.csv(x = reg_cnt, file = paste0(out_dir,'/agex_sensitivity_regions_change_continent.csv'), row.names = F)

# Top 5 countries with the largest change per comparison and hotspot type
reg |>
  dplyr::group_by(scenario_b, hotspot) |>
  dplyr::slice_max(order_by = lost_km2 + gained_km2, n = 5) |>
  dplyr::ungroup() |>
  base::as.data.frame()

# Pixel-level agreement maps (number of scenarios in which a pixel is hotspot) ----
# Replaces Figure S6 with the same scenarios used in the agreement table.
fls <- paste0(out_dir,'/agex_hotspots_',unique(c(baseline, cmp$scenario_b[cmp$group != 'Simplified pipeline (no outlier cleaning)'])),'.tif')
for(h in c('food_production_hotspot','economic_value_hotspot')){
  stk <- lapply(fls, function(f) terra::rast(f)[[h]]) |> terra::rast()
  stk <- terra::subst(stk, NA, 0)
  agr <- terra::app(stk, mean) # share of scenarios classifying the pixel as hotspot
  agr <- terra::mask(agr, terra::rast(fls[1])[['extreme_cluster']])
  names(agr) <- paste0(h,'_agreement')
  terra::writeRaster(x = agr, filename = paste0(out_dir,'/agex_sensitivity_agreement_',h,'.tif'), overwrite = T)
}; rm(h, stk, agr)
