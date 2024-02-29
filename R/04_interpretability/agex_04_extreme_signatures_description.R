# ------------------------------------------ #
# Descriptives of extreme weather signatures
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,geodata,NbClust,entropy,scales,dplyr,
                                ggplot2,RColorBrewer,landscapemetrics))

# Root directory
root <- '//CATALOGUE/WFP_ClimateRiskPr1'

# Load extreme weather signatures
agex_sgn <- terra::rast('D:/WORLD_fmado_trimmed_100_gs1.tif')
names(agex_sgn) <- 'extreme_signature'
crds_sgn <- terra::as.data.frame(x = agex_sgn, xy = T, cell = T, na.rm = T)

# Load global shapefile
wrl <- geodata::world(resolution = 3, level = 0, path = tempdir())

# Compute diversity of crop types
# High entropy means high variation of food/non-food classes
crp_fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_raw_data'), pattern = '^agex_cropclass_', full.names = T)
crp_typ <- terra::rast(crp_fls)
crp_ntp <- terra::app(x = crp_typ, function(x){entropy::entropy(x, method = 'ML')})
names(crp_ntp) <- 'crop_types_diversity'
plot(crp_ntp)
# terra::writeRaster(x = ent, filename = 'D:/crops_entropy.tif')
crp_ntp_25km <- terra::resample(x = crp_ntp, y = agex_sgn, method = 'bilinear')

# Crop types' diversity
terra::zonal(x = crp_ntp_25km, z = agex_sgn, fun = 'mean', na.rm = T) |> View()

## Cereals analysis
# cereals <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_cropclass_cereals_count_10km.tif'))
# cereals <- terra::resample(x = cereals, y = agex_sgn, method = 'near')
# # Signatures: 2 (few cereals); 45 (rich in cereals)

# Value of production
vop <- terra::rast('C:/Users/haachicanoy/Downloads/spam2010V2r0_global_V_agg_VP_CROP_A.tif')
names(vop) <- 'value_of_production'
vop_25km <- terra::resample(x = vop, y = agex_sgn, method = 'bilinear')
terra::zonal(x = vop_25km, z = agex_sgn, fun = 'mean', na.rm = T) |> View()

tst <- terra::zonal(x = c(crp_ntp_25km, vop_25km), z = agex_sgn, fun = 'mean', na.rm = T)
plot(tst$crop_types_diversity, tst$value_of_production, pch = 20)

tst |>
  ggplot2::ggplot(aes(x = crop_types_diversity, y = value_of_production)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(se = F)

# CDD Severity
cdd <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_cdd/agex_cdd_25km/one_s1_cdd_25km.tif'))
for(i in 1:(terra::nlyr(cdd))){
  cdd[[i]][is.infinite(cdd[[i]])] <- 0
}; rm(i); gc(T)

cdd_avg <- terra::app(x = cdd, mean)
cdd_max <- terra::app(x = cdd, max)
cdd_mdn <- terra::app(x = cdd, median)

tst <- terra::zonal(x = c(crp_ntp_25km, vop_25km, cdd_avg, cdd_max, cdd_mdn), z = agex_sgn, fun = 'mean', na.rm = T)

tst |>
  ggplot2::ggplot(aes(x = median, y = crop_types_diversity)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(se = F)

# Value of production per tertiles
vTert <- stats::quantile(tst$value_of_production, c(0:3/3))
tst$vop_tert = with(tst, cut(value_of_production, vTert, include.lowest = T, labels = c('Low', 'Medium', 'High')))
# CDD severity per tertiles
vTert <- stats::quantile(tst$median, c(0:3/3))
tst$cdd_tert = with(tst, cut(median, vTert, include.lowest = T, labels = c('Low', 'Medium', 'High')))

table(tst$cdd_tert, tst$vop_tert)

View(tst)

tst <- tst |> dplyr::arrange(value_of_production) |> base::as.data.frame()
tst$vop_colors <- NA
tst$vop_colors[tst$vop_tert == 'Low'] <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,'YlGn'))(sum(tst$vop_tert == 'Low'))
tst$vop_colors[tst$vop_tert == 'Medium'] <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,'Oranges'))(sum(tst$vop_tert == 'Medium'))
tst$vop_colors[tst$vop_tert == 'High'] <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,'Reds'))(sum(tst$vop_tert == 'High'))

tst <- tst |> dplyr::arrange(extreme_signature) |> base::as.data.frame()

# Population density
pop <- geodata::population(year = 2020, res = 10, path = tempdir())
pop_25km <- terra::resample(x = pop, y = agex_sgn, method = 'sum')

tst <- terra::zonal(x = c(crp_ntp_25km, vop_25km, cdd_avg, cdd_max, cdd_mdn, pop_25km), z = agex_sgn, fun = 'mean', na.rm = T)

tst[,-1] |> cor(method = 'spearman') |> round(2)

## Graph of global cluster
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- terra::aggregate(wrl)
tmp <- agex_sgn
tmp <- terra::trim(x = tmp)
hist(terra::values(tmp))
png(filename = "D:/global_cluster_3palettes.png", width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(tmp))
plot(tmp, add = T, col = tst$vop_colors, plg = list(cex = 5), cex.main = 7)
dev.off()


lapply(1:length(unique(crds_sgn$extreme_signature)), function(i){
  # Coordinates within countries
  crds_cty <- terra::intersect(x = wrl, terra::vect(crds_sgn[crds_sgn$extreme_signature == i,], c('x','y'), crs = 'EPSG:4326')) |> base::as.data.frame()
  # Get the list of countries
  countries_list <- unique(crds_cty$NAME_0)
  # Count the number of pixels within each country
  cnts_cty <- table(crds_cty$NAME_0) |> base::as.data.frame()
  names(cnts_cty)[1] <- 'Country'
  # Adjacent countries
  ctys_flt <- wrl[wrl$NAME_0 %in% countries_list,]
  adjc_cty <- terra::relate(x = ctys_flt, relation = 'touches')
  adjc_cty <- 1 * adjc_cty
  colnames(adjc_cty) <- rownames(adjc_cty) <- countries_list
  adjc_cty <- adjc_cty[apply(adjc_cty, 2, sd) != 0, apply(adjc_cty, 2, sd) != 0]
  # Best separation
  if(length(countries_list) > 2){
    sprt_cty <- NbClust::NbClust(diss     = as.dist(1-adjc_cty),
                                 distance = NULL,
                                 min.nc   = 2,
                                 max.nc   = nrow(adjc_cty)-1,
                                 method   = 'ward.D2',
                                 index    = c('silhouette'))
    sort(sprt_cty$Best.partition)
    # plot(hclust(as.dist(1-adjc_cty), method = 'ward.D2'), hang = -1)
  } else {
    # The country or two countries
  }
  
})

# Cohesion index: Equals 0 if patches of class i become more isolated.
# Increases if patches of class i become more aggregated
agex_chs <- landscapemetrics::lsm_c_cohesion(landscape = agex_sgn) |> base::as.data.frame()
View(agex_chs)
# Contiguity index mean: equals the mean of the contiguity index on class
# level for all patches
agex_ctg <- landscapemetrics::lsm_c_contig_mn(landscape = agex_sgn) |> base::as.data.frame()
View(agex_ctg)
# Aggregation index: Equals 0 for maximally disaggregated
# and 100 for maximally aggregated classes
landscapemetrics::lsm_c_ai(landscape = agex_sgn) |> base::as.data.frame() |> View()
