# ------------------------------------------ #
# Description of extreme weather clusters
# By: Harold Achicanoy
# WUR & ABC
# May 2024
# ------------------------------------------ #

## R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,geodata,NbClust,entropy,scales,dplyr,
                                ggplot2,RColorBrewer,landscapemetrics,
                                exactextractr,hrbrthemes,trend,quantreg,
                                OutliersO3,MetBrewer,extrafont))

## Relevant functions
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
    qrf <- quantreg::rq(formula = ts ~ time, tau = .90, data = dfm)
    y <- as.numeric(qrf$coefficients[2])
  } else { y <- NA }
  return(y)
}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'

## Extreme weather clusters
agex_sgn_cln <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_cln) <- 'extreme_cluster'
# Extreme weather clusters by growing seasons
agex_sgn_gs1 <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_one_s1_fmadogram_clean.tif'))
agex_sgn_gs1[!is.na(agex_sgn_gs1)] <- 1
agex_sgn_gs2 <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_two_fmadogram_clean.tif'))
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

## Countries collaboration
# Countries coverage across extreme drought clusters
out <- paste0(root,'/agroclimExtremes/agex_results/agex_country_collaborations.csv')
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
collaborations <- dplyr::left_join(x = collaborations, y = area_dfm, by = 'extreme_cluster'); rm(area_dfm)

## Cluster quality
# They are computed to provide a quality control of the produced clusters
# Metric: Cohesion index
# Meaning:
# - 0 if patches of class become more isolated
# - 100 if patches of class i become more aggregated
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

## Crop classes diversity
# Metric: Shannon diversity index
# Meaning: High entropy means high variation of food/non-food classes
crp_fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_raw_data'), pattern = '^agex_cropclass_', full.names = T)
crp_fls <- crp_fls[-grep2(pattern = c('fibres','stimulant','rest_of_crops'), x = crp_fls)] # Categories to exclude
crp_typ <- terra::rast(crp_fls)
crp_ntp <- terra::app(x = crp_typ, fun = function(i, ff) ff(i), cores = 20, ff = entropy::entropy)
names(crp_ntp) <- 'crop_classes_diversity'
plot(crp_ntp)
crp_ntp_25km <- terra::resample(x = crp_ntp, y = agex_sgn_cln, method = 'cubicspline', threads = T)
crp_ntp_25km <- terra::mask(x = crp_ntp_25km, mask = agex_sgn_cln)

## Agricultural economic value
# Metric: Value of production
# Meaning: High value of production means high economic value from agriculture
all_vop <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_CROP_A.tif'))
names(all_vop) <- 'total_vop'
all_vop_25km <- terra::resample(x = all_vop, y = agex_sgn_cln, method = 'cubicspline', threads = T)
all_vop_25km <- terra::mask(x = all_vop_25km, mask = agex_sgn_cln)

food_vop <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_FOOD_A.tif'))
names(food_vop) <- 'food_vop'
food_vop_25km <- terra::resample(x = food_vop, y = agex_sgn_cln, method = 'cubicspline', threads = T)
food_vop_25km <- terra::mask(x = food_vop_25km, mask = agex_sgn_cln)

nonf_vop <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_NONF_A.tif'))
names(nonf_vop) <- 'non-food_vop'
nonf_vop_25km <- terra::resample(x = nonf_vop, y = agex_sgn_cln, method = 'cubicspline', threads = T)
nonf_vop_25km <- terra::mask(x = nonf_vop_25km, mask = agex_sgn_cln)

## Livestock diversity
# Metric: Livestock equivalent units
# Meaning: grazing equivalent of one adult dairy cow producing 3 000 kg of milk
# annually, without additional concentrated foodstuffs
lvs_dir <- paste0(root,'/agroclimExtremes/agex_raw_data/livestock')
anmls   <- list.dirs(path = lvs_dir, full.names = F, recursive = F)
anmls   <- anmls[-grep('horses',anmls)]
lvstc_fls <- list.files2(path = paste0(lvs_dir,'/',anmls), pattern = '_Da.tif$', full.names = T); rm(lvs_dir, anmls)
names(lvstc_fls) <- NULL
lvstc_cnt <- terra::rast(lvstc_fls)

# Livestock Units
# LU: Livestock Unit
# LU = Buffaloes * 1 + Cattle * 1 + Chickens * ((0.007 + 0.014)/2) +
#      Ducks * 0.01 + Goats * 0.1 + Pigs * ((0.5+0.027)/2) + Sheep * 0.1
# Source: https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Livestock_unit_(LSU)
wghts <- c(1, 1, (0.007 + 0.014)/2, 0.01, 0.1, (0.5+0.027)/2, 0.1)
lvstc_unt <- lvstc_cnt * wghts
lsu <- sum(lvstc_unt)
rm(lvstc_fls, lvstc_cnt)
names(lsu) <- 'total_livestock_units'
lsu_25km <- terra::resample(x = lsu, y = agex_sgn_cln, method = 'cubicspline', threads = T)
lsu_25km <- terra::mask(x = lsu_25km, mask = agex_sgn_cln)

## Exposed population
# Metric: Population density
# Meaning: number of people per unit of area
pop <- geodata::population(year = 2020, res = 10, path = tempdir())
pop_25km <- terra::resample(x = pop, y = agex_sgn_cln, method = 'cubicspline', threads = T)
pop_25km <- terra::mask(x = pop_25km, mask = agex_sgn_cln)

## Index severity
# Metric: SPEI-6 trend for 90th percentile. We also computed the number of years when SPEI < -1.5
# For two-growing-seasons pixels we computed the average of the trends
# Meaning: extreme drought trend over time
stp <- data.frame(gs = c('one','two','two'), season = c('1','1','2'))
idx_one_s1 <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/',stp$gs[1],'_s',stp$season[1],'_',index,'_25km.tif'))
idx_two_s1 <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/',stp$gs[2],'_s',stp$season[2],'_',index,'_25km.tif'))
idx_two_s2 <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/',stp$gs[3],'_s',stp$season[3],'_',index,'_25km.tif'))

# # Trend
# idx_one_s1_slp <- terra::app(x = idx_one_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_trend); names(idx_one_s1_slp) <- 'SPEI-6_slope'
# idx_two_s1_slp <- terra::app(x = idx_two_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_trend); names(idx_two_s1_slp) <- 'SPEI-6_slope'
# idx_two_s2_slp <- terra::app(x = idx_two_s2, fun = function(i, ff) ff(i), cores = 20, ff = get_trend); names(idx_two_s2_slp) <- 'SPEI-6_slope'

# Extreme trend
idx_one_s1_sxt <- terra::app(x = idx_one_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_ext_trend); names(idx_one_s1_sxt) <- 'SPEI-6_slope_90th'
idx_two_s1_sxt <- terra::app(x = idx_two_s1, fun = function(i, ff) ff(i), cores = 20, ff = get_ext_trend); names(idx_two_s1_sxt) <- 'SPEI-6_slope_90th'
idx_two_s2_sxt <- terra::app(x = idx_two_s2, fun = function(i, ff) ff(i), cores = 20, ff = get_ext_trend); names(idx_two_s2_sxt) <- 'SPEI-6_slope_90th'
idx_two_sxt <- mean(c(idx_two_s1_sxt, idx_two_s2_sxt))

idx_sxt <- terra::merge(idx_one_s1_sxt, idx_two_sxt)

# Extreme drought years' count
idx_one_s1_cnt <- terra::app(x = idx_one_s1, fun = function(i, ff) ff(i), cores = 20, ff = function(x){sum(x > 1.5)}); names(idx_one_s1_cnt) <- 'SPEI-6_extreme_drought_count'
idx_two_s1_cnt <- terra::app(x = idx_two_s1, fun = function(i, ff) ff(i), cores = 20, ff = function(x){sum(x > 1.5)}); names(idx_two_s1_cnt) <- 'SPEI-6_extreme_drought_count'
idx_two_s2_cnt <- terra::app(x = idx_two_s2, fun = function(i, ff) ff(i), cores = 20, ff = function(x){sum(x > 1.5)}); names(idx_two_s2_cnt) <- 'SPEI-6_extreme_drought_count'
idx_two_cnt <- round(mean(c(idx_two_s1_cnt, idx_two_s2_cnt)))

idx_cnt <- terra::merge(idx_one_s1_cnt, idx_two_cnt)

rst_mtrcs <- terra::zonal(x = c(pop_25km, vop_25km, crp_ntp_25km, lsu_25km, idx_sxt, idx_cnt),
                          z = agex_sgn_cln, fun = 'mean', na.rm = T)
rst_mtrcs[,-1] <- round(rst_mtrcs[,-1], 3)

dfm <- dplyr::left_join(x = collaborations, y = rst_mtrcs, by = 'extreme_cluster')
dfm <- dplyr::left_join(x = dfm, y = qlt_mtrcs, by = 'extreme_cluster')

hist(dfm$`SPEI-6_slope_90th`)
abline(v = quantile(x = dfm$`SPEI-6_slope_90th`, probs = seq(0,1,1/4)), lty = 2, col = 'red')

hist(dfm$crop_classes_diversity)
abline(v = quantile(x = dfm$crop_classes_diversity, probs = seq(0,1,1/4)), lty = 2, col = 'red')

hist(dfm$livestock_units)
abline(v = quantile(x = dfm$livestock_units, probs = seq(0,1,1/4)), lty = 2, col = 'red')

agex_sgn_poly <- terra::as.polygons(x = agex_sgn_cln)
agex_sgn_poly <- terra::merge(x = agex_sgn_poly, y = dfm)

terra::writeVector(x = agex_sgn_poly, filename = paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/vct_agex_global_',index,'.gpkg'), overwrite = T)








# Figure 1
tst <- biscale::bi_class(.data = dfm, x = crop_types_diversity, y = `SPEI-6_slope`, style = 'quantile', dim = 3)
colours <- data.frame(bi_class = sort(unique(tst$bi_class)))
colours$col <- pals::stevens.bluered(n = nrow(colours))
tst <- dplyr::left_join(x = tst, y = colours, by = 'bi_class') |> base::as.data.frame()

tst <- dplyr::arrange(.data = tst, -`SPEI-6_slope`) |> base::as.data.frame()
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'OKeeffe1', n = n_col)
tst$palette <- col_pltt
tst <- dplyr::arrange(.data = tst, extreme_cluster) |> base::as.data.frame()

auxl_plt <- agex_sgn
auxl_plt <- terra::trim(x = auxl_plt)

png(filename = paste0(root,'/agroclimExtremes/agex_global_',index,'_',gs,'_s',season,'_fmadogram_bivariate_crops.png'), width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(auxl_plt))
plot(auxl_plt, add = T, col = tst$col, plg = list(cex = 5), cex.main = 7)
dev.off()

# Figure 2
tst <- biscale::bi_class(.data = dfm, x = livestock_units, y = `SPEI-6_slope`, style = 'quantile', dim = 3)
colours <- data.frame(bi_class = sort(unique(tst$bi_class)))
colours$col <- pals::stevens.bluered(n = nrow(colours))
tst <- dplyr::left_join(x = tst, y = colours, by = 'bi_class') |> base::as.data.frame()

tst <- dplyr::arrange(.data = tst, -`SPEI-6_slope`) |> base::as.data.frame()
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'OKeeffe1', n = n_col)
tst$palette <- col_pltt
tst <- dplyr::arrange(.data = tst, extreme_cluster) |> base::as.data.frame()

auxl_plt <- agex_sgn
auxl_plt <- terra::trim(x = auxl_plt)

png(filename = paste0(root,'/agroclimExtremes/agex_global_',index,'_',gs,'_s',season,'_fmadogram_bivariate_livestock.png'), width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(auxl_plt))
plot(auxl_plt, add = T, col = tst$col, plg = list(cex = 5), cex.main = 7)
dev.off()

# Legend box
tst$x <- as.numeric(substr(x = tst$bi_class, start=1, stop=1))
tst$y <- as.numeric(substr(x = tst$bi_class, start=3, stop=3))
tst %>%
  ggplot2::ggplot(aes(x = x, y = y)) +
  ggplot2::geom_tile(fill = tst$col, alpha = 1) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal() +
  ggplot2::xlab('Livestock diversity') +
  ggplot2::ylab('Drought trend') +
  ggplot2::theme(axis.text       = element_blank(),
                 axis.title      = element_text(size = 35),
                 legend.text     = element_text(size = 17),
                 legend.title    = element_blank(),
                 plot.title      = element_text(size = 25),
                 plot.subtitle   = element_text(size = 17),
                 plot.caption    = element_text(size = 15, hjust = 0),
                 legend.position = "bottom")

# Figure 3
dfm |>
  dplyr::select(extreme_cluster,`SPEI-6_slope`,`SPEI-6_slope_90th`) |>
  tidyr::pivot_longer(cols = -1, names_to = 'variable', values_to = 'value') |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = value, fill = variable)) +
  ggplot2::geom_density(color="#e9ecef", alpha=0.6, position = 'identity') +
  ggplot2::scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ggplot2::geom_vline(xintercept=0, linetype=3, colour='red') +
  ggplot2::theme_bw()


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



tst <- terra::zonal(x = c(crp_ntp_25km, vop_25km, idx_avg, idx_max, idx_mdn, pop_25km), z = agex_sgn, fun = 'mean', na.rm = T)
tst[,2:ncol(tst)] <- round(tst[,2:ncol(tst)], 2)

tst[,-1] |> cor(method = 'kendall') |> round(2)

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

biscale::bi_class(.data = tbl, x = Heat, y = Waterlogging, style = "quantile", dim = 3)
